use std::collections::HashMap;

use crate::cli::ModelType;
use crate::grouping::ReadGroup;
use crate::intern::Symbol;
use crate::matrix::{MatrixCell, MatrixRow, row_get};

/// Result of ambiguity resolution for one read group.
pub struct ResolutionResult {
    /// Annotation → alignment index for confident assignments only.
    pub confident: HashMap<Symbol, usize>,
    /// Annotation → softmax probability for all annotations.
    #[allow(dead_code)]
    pub probabilities: HashMap<Symbol, f64>,
}

/// For each annotation that any alignment in the group overlaps,
/// pick the best representative alignment using tiebreaking order:
///   1. Highest AS (alignment score) tag value
///   2. Proper pair flag (prefer true)
///   3. Highest overlap with this specific annotation
///   4. First in input order (lowest input_order)
///
/// Alignments with no GTF overlap contribute a synthetic `__no_feature__`
/// candidate (overlap = 0), so EM can weigh "originated from no feature"
/// against real annotations.
///
/// Returns a sparse matrix row: a Vec<(Symbol, MatrixCell)> sorted ascending
/// by symbol so downstream lookups can binary-search.
pub fn select_representatives(
    group: &ReadGroup,
    no_feature_symbol: Symbol,
) -> Vec<(Symbol, MatrixCell)> {
    // Collect all unique annotations across all alignments
    #[allow(clippy::type_complexity)]
    let mut annotation_candidates: HashMap<Symbol, Vec<(usize, i64, bool, usize, usize)>> =
        HashMap::new();

    for (idx, entry) in group.alignments.iter().enumerate() {
        if entry.annotations.is_empty() {
            annotation_candidates
                .entry(no_feature_symbol)
                .or_default()
                .push((
                    idx,
                    entry.alignment_score.unwrap_or(i64::MIN),
                    entry.is_proper_pair,
                    0,
                    entry.input_order,
                ));
        } else {
            for ao in &entry.annotations {
                annotation_candidates
                    .entry(ao.annotation)
                    .or_default()
                    .push((
                        idx,
                        entry.alignment_score.unwrap_or(i64::MIN),
                        entry.is_proper_pair,
                        ao.overlap_bp,
                        entry.input_order,
                    ));
            }
        }
    }

    let mut row: Vec<(Symbol, MatrixCell)> = Vec::with_capacity(annotation_candidates.len());

    for (annotation, mut candidates) in annotation_candidates {
        // Sort by: AS desc, proper_pair desc, overlap desc, input_order asc
        candidates.sort_by(|a, b| {
            b.1.cmp(&a.1) // AS desc
                .then_with(|| b.2.cmp(&a.2)) // proper_pair desc (true > false)
                .then_with(|| b.3.cmp(&a.3)) // overlap desc
                .then_with(|| a.4.cmp(&b.4)) // input_order asc
        });

        let winner = candidates[0];
        row.push((
            annotation,
            MatrixCell {
                alignment_idx: winner.0,
                // Store unwrap_or(0) for softmax (i64::MIN was only for tiebreaking)
                alignment_score: group.alignments[winner.0].alignment_score.unwrap_or(0),
            },
        ));
    }

    row.sort_unstable_by_key(|(s, _)| s.idx());
    row
}

/// Per-fragment prior probabilities, in a dense layout keyed by position in
/// `annotations` (so EM can iterate without HashMap hashing costs).
///
/// Invariants:
///   - `annotations.len() == priors.len()`
///   - if non-empty, `priors` sums to 1.0
///   - for a no-feature fragment (empty row), both are empty.
#[derive(Debug, Clone)]
pub struct FragmentPriors {
    pub annotations: Vec<Symbol>,
    pub priors: Vec<f64>,
}

/// Per-fragment EM posteriors. Carries the annotation list alongside the
/// probability vector so the priors Vec can be dropped immediately after EM
/// finishes — saves ~3 GB on the 28M-group dataset.
#[derive(Debug, Clone)]
pub struct FragmentPosteriors {
    pub annotations: Vec<Symbol>,
    pub probs: Vec<f64>,
}

/// Result of an EM run: the per-fragment posteriors plus the number of E-steps
/// actually performed. `iterations` is 0 when EM was skipped (`max_iter == 0`,
/// or — for telescope — no ambiguous fragments to iterate); otherwise it is the
/// iteration at which the max-mass-delta convergence criterion tripped, capped
/// at `max_iter`.
#[derive(Debug, Clone)]
pub struct EmOutcome {
    pub posteriors: Vec<FragmentPosteriors>,
    pub iterations: usize,
}

/// Sum each fragment's probability vector into per-annotation totals, indexed
/// by [`Symbol::idx`]. Each fragment contributes its full distribution (summing
/// to 1, or 0 when empty) unweighted, so the result is the expected fractional
/// fragment count per annotation. Applied to priors it gives the prior
/// assignment; applied to the finalized posteriors it gives ribbonfish θ
/// (`Σ_r γ_{r,a}`) / telescope `Σγ`. The per-iteration form of this same vector
/// (computed by [`accumulate_mass`]) drives the max-delta EM early-stopping
/// criterion.
pub fn fractional_mass<'a>(
    fragments: impl IntoIterator<Item = (&'a [Symbol], &'a [f64])>,
    interner_size: usize,
) -> Vec<f64> {
    let mut mass = vec![0.0; interner_size];
    for (anns, probs) in fragments {
        for (&a, &p) in anns.iter().zip(probs) {
            mass[a.idx()] += p;
        }
    }
    mass
}

/// Compute per-fragment prior probabilities from representative AS scores.
/// This is the model-specific step; the resulting priors are fed into EM.
pub fn compute_priors(
    row: &MatrixRow,
    model: ModelType,
    global_as_min: Option<i64>,
    global_as_max: Option<i64>,
) -> FragmentPriors {
    if row.is_empty() {
        return FragmentPriors {
            annotations: Vec::new(),
            priors: Vec::new(),
        };
    }

    let annotations: Vec<Symbol> = row.iter().map(|(s, _)| *s).collect();
    let scores: Vec<f64> = row.iter().map(|(_, c)| c.alignment_score as f64).collect();

    let priors: Vec<f64> = match model {
        ModelType::Ribbonfish => {
            // Numerically stable softmax with temperature 5
            let max_score = scores.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
            let exps: Vec<f64> = scores
                .iter()
                .map(|s| ((s - max_score) / 5.0).exp())
                .collect();
            let sum: f64 = exps.iter().sum();
            exps.iter().map(|e| e / sum).collect()
        }
        ModelType::Telescope => {
            let n = scores.len();
            match (global_as_min, global_as_max) {
                (Some(lo), Some(hi)) if hi > lo => {
                    let lo_f = lo as f64;
                    let range = (hi - lo) as f64;
                    let weights: Vec<f64> = scores
                        .iter()
                        .map(|s| ((s - lo_f) / range * 100.0).exp_m1())
                        .collect();
                    let sum: f64 = weights.iter().sum();
                    if sum > 0.0 {
                        weights.iter().map(|w| w / sum).collect()
                    } else {
                        vec![1.0 / n as f64; n]
                    }
                }
                _ => vec![1.0 / n as f64; n],
            }
        }
    };

    FragmentPriors {
        annotations,
        priors,
    }
}

/// Per-fragment M-step weight for the telescope model: the *best* (max)
/// alignment score in the row, transformed the telescope way
/// (`expm1((AS − lo)/(hi − lo) · 100)`), taken before normalization to a
/// prior and as the max (so having multiple alignments never reduces it).
///
/// Falls back to `1.0` when global AS bounds are unavailable or degenerate
/// (`hi <= lo`), making the M-step a plain unweighted count. A fragment whose
/// best score sits at the global minimum gets weight `0.0` (consistent with
/// telescope's "global min contributes zero mass").
pub fn fragment_weight(
    row: &MatrixRow,
    global_as_min: Option<i64>,
    global_as_max: Option<i64>,
) -> f64 {
    if row.is_empty() {
        return 1.0;
    }
    match (global_as_min, global_as_max) {
        (Some(lo), Some(hi)) if hi > lo => {
            let max_s = row
                .iter()
                .map(|(_, c)| c.alignment_score)
                .max()
                .unwrap_or(lo);
            let lo_f = lo as f64;
            let range = (hi - lo) as f64;
            ((max_s as f64 - lo_f) / range * 100.0).exp_m1()
        }
        _ => 1.0,
    }
}

/// The largest weight any fragment can be assigned given the global AS bounds:
/// a perfectly-aligned fragment (best score == global AS max) scales to 100, so
/// its weight is `expm1(100)`. In the degenerate case (no usable bounds, or
/// `hi <= lo` — e.g. every read shares one AS value) [`fragment_weight`] falls
/// back to `1.0` for *every* fragment, so the maximum is likewise `1.0`. Used
/// to weight the `theta_prior` regularizing reads as ideal alignments without
/// over-weighting them relative to genuine reads in the degenerate case.
pub fn max_fragment_weight(global_as_min: Option<i64>, global_as_max: Option<i64>) -> f64 {
    match (global_as_min, global_as_max) {
        (Some(lo), Some(hi)) if hi > lo => 100.0_f64.exp_m1(),
        _ => 1.0,
    }
}

/// Run expectation-maximization over per-fragment priors, returning per-fragment
/// posterior γ. θ (read distribution over annotations) is seeded so each non-empty
/// fragment contributes `1/|R_r|` uniformly to each of its annotations.
///
/// Runs at most `max_iter` iterations, stopping early once the largest
/// per-annotation change in assigned fragment mass (`Σ_r γ_{r,a}`) between
/// consecutive iterations drops below `min_mass_delta` (set it to 0 to disable
/// early stopping). The returned [`EmOutcome`] reports the number of E-steps
/// actually performed. When `max_iter == 0`, returns γ == priors unchanged
/// (`iterations == 0`) — callers get the same behavior as running the confidence
/// pick directly on the priors.
///
/// `lengths[a.idx()]` is the per-annotation total length in bp. When
/// `length_correction` is true, the E-step weight becomes `(θ_a / L_a) · π_{r,a}`,
/// biasing posterior mass toward shorter annotations on otherwise-equal evidence.
/// Length 0 (unknown) falls back to no correction for that annotation.
pub fn expectation_maximization(
    priors: Vec<FragmentPriors>,
    interner_size: usize,
    max_iter: usize,
    min_mass_delta: f64,
    lengths: &[u64],
    length_correction: bool,
) -> EmOutcome {
    // γ starts as a copy of the priors so that max_iter=0 yields posteriors == priors.
    let mut gammas: Vec<Vec<f64>> = priors.iter().map(|fp| fp.priors.clone()).collect();

    if max_iter == 0 {
        return EmOutcome {
            posteriors: finalize(priors, gammas),
            iterations: 0,
        };
    }

    // θ_a^(0): each read splits 1.0 uniformly across its own annotations.
    let mut theta = vec![0.0; interner_size];
    for fp in &priors {
        if fp.annotations.is_empty() {
            continue;
        }
        let w = 1.0 / fp.annotations.len() as f64;
        for &a in &fp.annotations {
            theta[a.idx()] += w;
        }
    }

    // Precompute 1/L per annotation if length correction is enabled (empty
    // otherwise). Shared with `telescope_em` so the two E-steps stay in sync.
    let inv_lengths = inverse_length_weights(interner_size, lengths, length_correction);

    // Reused across iterations to avoid per-iteration interner_size-sized allocations.
    let mut theta_next = vec![0.0; interner_size];

    // Convergence tracking: per-annotation assigned mass Σγ, seeded from the
    // initial γ (== priors) so iteration 1 can already converge. Buffers are
    // swapped each iteration to avoid reallocation.
    let mut prev_mass = vec![0.0; interner_size];
    let mut cur_mass = vec![0.0; interner_size];
    accumulate_mass(&mut prev_mass, &priors, &gammas);
    let mut iterations = 0;

    // EM is sequential by design: per-fragment work is a handful of float ops, well below
    // rayon's per-task overhead, and the M-step's full-θ accumulation is bandwidth-bound.
    // Empirically, --threads 1 was ~10× faster than the parallel version on real BAMs.
    for _ in 0..max_iter {
        // E-step: γ_{r,a} ∝ θ_a · π_{r,a}, normalized per fragment.
        // With length correction: γ_{r,a} ∝ (θ_a / L_a) · π_{r,a}.
        for (fp, gamma) in priors.iter().zip(gammas.iter_mut()) {
            if fp.annotations.is_empty() {
                continue;
            }
            let mut sum = 0.0;
            for (i, (&a, &p)) in fp.annotations.iter().zip(fp.priors.iter()).enumerate() {
                let theta_eff = if length_correction {
                    theta[a.idx()] * inv_lengths[a.idx()]
                } else {
                    theta[a.idx()]
                };
                let w = theta_eff * p;
                gamma[i] = w;
                sum += w;
            }
            if sum > 0.0 {
                for g in gamma.iter_mut() {
                    *g /= sum;
                }
            } else {
                // Degenerate: every θ·π product is 0 for this fragment.
                // Fall back to uniform across the fragment's annotations.
                let n = gamma.len() as f64;
                for g in gamma.iter_mut() {
                    *g = 1.0 / n;
                }
            }
        }
        iterations += 1;

        // Convergence check on the just-updated γ. If converged we return these
        // γ directly and skip the (now redundant) M-step.
        accumulate_mass(&mut cur_mass, &priors, &gammas);
        if max_abs_delta(&cur_mass, &prev_mass) < min_mass_delta {
            break;
        }

        // M-step: θ_a = Σ_r γ_{r,a}. Single accumulator buffer reused across iterations.
        // (Length correction lives entirely in the E-step; the M-step accumulates
        // raw γ-mass per annotation, just like the uncorrected EM.)
        theta_next.fill(0.0);
        for (fp, gamma) in priors.iter().zip(gammas.iter()) {
            for (&a, &g) in fp.annotations.iter().zip(gamma.iter()) {
                theta_next[a.idx()] += g;
            }
        }
        std::mem::swap(&mut theta, &mut theta_next);
        std::mem::swap(&mut prev_mass, &mut cur_mass);
    }

    // Move annotation lists out of priors into the returned posteriors; the
    // priors' π Vecs and the FragmentPriors structs themselves are dropped here.
    // This is the savings: the caller no longer needs to retain a separate
    // priors Vec across resolution / output phases.
    EmOutcome {
        posteriors: finalize(priors, gammas),
        iterations,
    }
}

fn finalize(priors: Vec<FragmentPriors>, gammas: Vec<Vec<f64>>) -> Vec<FragmentPosteriors> {
    priors
        .into_iter()
        .zip(gammas)
        .map(|(fp, probs)| FragmentPosteriors {
            annotations: fp.annotations,
            probs,
        })
        .collect()
}

/// Accumulate the per-annotation assigned fragment mass `Σ_r γ_{r,a}` into the
/// caller-owned `buf` (cleared first), index-parallel `priors`/`gammas`. This is
/// the per-iteration form of [`fractional_mass`]; reusing a buffer across
/// iterations avoids an `interner_size` allocation per EM step.
fn accumulate_mass(buf: &mut [f64], priors: &[FragmentPriors], gammas: &[Vec<f64>]) {
    buf.fill(0.0);
    for (fp, gamma) in priors.iter().zip(gammas.iter()) {
        for (&a, &g) in fp.annotations.iter().zip(gamma.iter()) {
            buf[a.idx()] += g;
        }
    }
}

/// Largest absolute element-wise difference between two equal-length mass
/// vectors — the EM convergence statistic `max_a |a[a] − b[a]|`.
fn max_abs_delta(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b.iter())
        .map(|(x, y)| (x - y).abs())
        .fold(0.0_f64, f64::max)
}

/// Normalize the entries at `support` indices so they sum to 1; if the sum is
/// ≤ 0, fall back to uniform over the support.
fn normalize_support(v: &mut [f64], support: &[usize]) {
    let total: f64 = support.iter().map(|&a| v[a]).sum();
    if total > 0.0 {
        for &a in support {
            v[a] /= total;
        }
    } else if !support.is_empty() {
        let u = 1.0 / support.len() as f64;
        for &a in support {
            v[a] = u;
        }
    }
}

/// Precompute the per-annotation inverse-length weights `1/L_a` used by the
/// length-corrected E-step (shared by [`expectation_maximization`] and
/// [`telescope_em`]). Returns an empty Vec when `length_correction` is off so
/// the E-step skips the lookup entirely. Symbols with length 0 (e.g., never
/// seen in the GTF) get weight 1.0 so the E-step degenerates to the
/// no-correction case for those entries.
fn inverse_length_weights(interner_size: usize, lengths: &[u64], length_correction: bool) -> Vec<f64> {
    if length_correction {
        (0..interner_size)
            .map(|i| {
                let len = lengths.get(i).copied().unwrap_or(0);
                if len == 0 { 1.0 } else { 1.0 / len as f64 }
            })
            .collect()
    } else {
        Vec::new()
    }
}

/// Normalize the whole vector so its entries sum to 1; if the sum is ≤ 0 the
/// vector is left unchanged (the E-step's degenerate branch then takes over).
fn normalize_all(v: &mut [f64]) {
    let total: f64 = v.iter().sum();
    if total > 0.0 {
        for x in v.iter_mut() {
            *x /= total;
        }
    }
}

/// Telescope-model EM. Splits read groups into unique (single-annotation rows)
/// and ambiguous (≥2-annotation rows) sets. Unique groups contribute a fixed
/// weighted term to π and are never iterated. Ambiguous groups are reassigned
/// via `γ_{r,a} ∝ prior_{r,a} · π_a · θ_a`; the M-step accumulates weighted
/// γ-mass into π (plus the fixed unique term) and θ (plus the `theta_prior`
/// regularizer on each annotation under contention).
///
/// The `theta_prior` regularizing reads are treated as though they aligned
/// *perfectly*: each is weighted by [`max_fragment_weight`] for the given AS
/// bounds rather than by a raw unit count — `expm1(100)` normally, but `1.0`
/// in the degenerate case so they are weighted equally to genuine fragments
/// (which also fall back to `1.0` there).
///
/// `weights[i]` is the [`fragment_weight`] of fragment `i` and is
/// index-parallel to `priors`. The returned posteriors are index-aligned with
/// `priors` (so `main.rs`'s position-based mapping back to read groups,
/// resolutions and output is preserved).
///
/// When `length_correction` is on, the E-step substitutes `θ_a / L_a` for
/// `θ_a` (RSEM-style effective-length normalization) using `lengths[a]`,
/// exactly as in [`expectation_maximization`]; the M-step and seeds are
/// untouched. Symbols with length 0 are treated as uncorrected.
///
/// Like [`expectation_maximization`], runs at most `max_iter` iterations,
/// stopping early once the largest per-annotation change in assigned fragment
/// mass (`Σ_r γ_{r,a}`, unweighted) between consecutive iterations drops below
/// `min_mass_delta` (0 disables early stopping); the returned [`EmOutcome`]
/// reports the iterations performed. `max_iter == 0` — or no ambiguous fragments
/// at all — returns γ == priors with `iterations == 0`.
#[allow(clippy::too_many_arguments)]
pub fn telescope_em(
    priors: Vec<FragmentPriors>,
    weights: &[f64],
    interner_size: usize,
    max_iter: usize,
    min_mass_delta: f64,
    theta_prior: f64,
    global_as_min: Option<i64>,
    global_as_max: Option<i64>,
    lengths: &[u64],
    length_correction: bool,
) -> EmOutcome {
    // γ starts as a copy of priors so max_iter=0 (or all-unique input) yields
    // posteriors == priors. Unique fragments keep this γ ([1.0]); empty stay [].
    let mut gammas: Vec<Vec<f64>> = priors.iter().map(|fp| fp.priors.clone()).collect();

    // Ambiguous = ≥2 candidate annotations. Unique (len == 1) and empty rows
    // are never iterated.
    let ambiguous: Vec<usize> = (0..priors.len())
        .filter(|&i| priors[i].annotations.len() >= 2)
        .collect();

    if max_iter == 0 || ambiguous.is_empty() {
        return EmOutcome {
            posteriors: finalize(priors, gammas),
            iterations: 0,
        };
    }

    // Precompute 1/L per annotation if length correction is enabled (empty
    // otherwise). Shared with `expectation_maximization` so the two E-steps
    // stay in sync; only the E-step consults it.
    let inv_lengths = inverse_length_weights(interner_size, lengths, length_correction);

    // π_fixed: weighted unique-read mass per annotation, constant across
    // iterations. A unique fragment maps to exactly one annotation.
    let mut pi_fixed = vec![0.0_f64; interner_size];
    for (i, fp) in priors.iter().enumerate() {
        if fp.annotations.len() == 1 {
            pi_fixed[fp.annotations[0].idx()] += weights[i];
        }
    }

    // θ/π support = annotation indices appearing in ≥1 ambiguous fragment.
    // The theta_prior regularizer is applied to exactly these annotations.
    let mut support: Vec<usize> = Vec::new();
    {
        let mut seen = vec![false; interner_size];
        for &r in &ambiguous {
            for &a in &priors[r].annotations {
                if !seen[a.idx()] {
                    seen[a.idx()] = true;
                    support.push(a.idx());
                }
            }
        }
    }

    // Weighted-uniform ambiguous seed: Σ_{r∈Amb} w_r/|R_r| onto each a ∈ R_r.
    let mut amb_seed = vec![0.0_f64; interner_size];
    for &r in &ambiguous {
        let fp = &priors[r];
        let share = weights[r] / fp.annotations.len() as f64;
        for &a in &fp.annotations {
            amb_seed[a.idx()] += share;
        }
    }

    // θ⁰ = normalize(amb_seed) over the support.
    let mut theta = amb_seed.clone();
    normalize_support(&mut theta, &support);

    // π⁰ = normalize(π_fixed + amb_seed) over all annotations. Any global
    // scaling of π/θ cancels in the per-fragment E-step ratio, and π is fully
    // recomputed each M-step, so this seed only shapes iteration 1.
    let mut pi = amb_seed;
    for (a, &pf) in pi_fixed.iter().enumerate() {
        pi[a] += pf;
    }
    normalize_all(&mut pi);

    // Each theta_prior regularizing read is treated as a perfectly-aligned
    // fragment: weighted by the maximum possible fragment weight given the AS
    // bounds — expm1(100) normally, or 1.0 in the degenerate case so it is
    // weighted equally to genuine ambiguous reads (which are all 1.0 there).
    let theta_reg = theta_prior * max_fragment_weight(global_as_min, global_as_max);

    let mut theta_next = vec![0.0_f64; interner_size];
    let mut pi_next = vec![0.0_f64; interner_size];

    // Convergence tracking on unweighted assigned mass Σγ (matching
    // `fractional_mass` / the run summary), seeded from the initial γ. Unique and
    // empty fragments contribute a constant term that cancels in the delta.
    let mut prev_mass = vec![0.0; interner_size];
    let mut cur_mass = vec![0.0; interner_size];
    accumulate_mass(&mut prev_mass, &priors, &gammas);
    let mut iterations = 0;

    for _ in 0..max_iter {
        // E-step (ambiguous fragments only): γ_{r,a} ∝ prior · π_a · θ_a.
        // With length correction: γ_{r,a} ∝ prior · π_a · (θ_a / L_a).
        for &r in &ambiguous {
            let fp = &priors[r];
            let gamma = &mut gammas[r];
            let mut sum = 0.0;
            for (i, (&a, &p)) in fp.annotations.iter().zip(fp.priors.iter()).enumerate() {
                let theta_eff = if length_correction {
                    theta[a.idx()] * inv_lengths[a.idx()]
                } else {
                    theta[a.idx()]
                };
                let w = p * pi[a.idx()] * theta_eff;
                gamma[i] = w;
                sum += w;
            }
            if sum > 0.0 {
                for g in gamma.iter_mut() {
                    *g /= sum;
                }
            } else {
                // Degenerate: every prior·π·θ product is 0 for this fragment.
                let n = gamma.len() as f64;
                for g in gamma.iter_mut() {
                    *g = 1.0 / n;
                }
            }
        }
        iterations += 1;

        // Convergence check on the just-updated γ. If converged we return these
        // γ directly and skip the (now redundant) M-step.
        accumulate_mass(&mut cur_mass, &priors, &gammas);
        if max_abs_delta(&cur_mass, &prev_mass) < min_mass_delta {
            break;
        }

        // M-step. raw_θ[a] = theta_prior·expm1(100) + Σ_{r∈Amb} w_r·γ_{r,a};
        //         raw_π[a] = π_fixed[a]            + Σ_{r∈Amb} w_r·γ_{r,a}.
        theta_next.fill(0.0);
        pi_next.fill(0.0);
        for &r in &ambiguous {
            let fp = &priors[r];
            let w = weights[r];
            for (&a, &g) in fp.annotations.iter().zip(gammas[r].iter()) {
                let mass = w * g;
                theta_next[a.idx()] += mass;
                pi_next[a.idx()] += mass;
            }
        }
        for &a in &support {
            theta_next[a] += theta_reg;
        }
        for (a, &pf) in pi_fixed.iter().enumerate() {
            pi_next[a] += pf;
        }
        normalize_support(&mut theta_next, &support);
        normalize_all(&mut pi_next);
        std::mem::swap(&mut theta, &mut theta_next);
        std::mem::swap(&mut pi, &mut pi_next);
        std::mem::swap(&mut prev_mass, &mut cur_mass);
    }

    EmOutcome {
        posteriors: finalize(priors, gammas),
        iterations,
    }
}

/// Apply the confidence threshold to EM posteriors and pick the single best
/// annotation per fragment (tiebreakers: AS desc, input_order asc).
pub fn select_confident(
    group: &ReadGroup,
    row: &MatrixRow,
    posteriors: &FragmentPosteriors,
    confidence: f64,
) -> ResolutionResult {
    let mut probabilities = HashMap::new();
    for (&ann, &prob) in posteriors.annotations.iter().zip(posteriors.probs.iter()) {
        probabilities.insert(ann, prob);
    }

    let mut confident = HashMap::new();
    if posteriors.annotations.is_empty() {
        return ResolutionResult {
            confident,
            probabilities,
        };
    }

    let best = posteriors
        .annotations
        .iter()
        .zip(posteriors.probs.iter())
        .max_by(|(a_ann, a_prob), (b_ann, b_prob)| {
            a_prob.partial_cmp(b_prob).unwrap().then_with(|| {
                let a_entry = &group.alignments[row_get(row, **a_ann).unwrap().alignment_idx];
                let b_entry = &group.alignments[row_get(row, **b_ann).unwrap().alignment_idx];
                // Tiebreak: AS desc, input_order asc
                a_entry
                    .alignment_score
                    .unwrap_or(0)
                    .cmp(&b_entry.alignment_score.unwrap_or(0))
                    .then_with(|| b_entry.input_order.cmp(&a_entry.input_order))
            })
        });

    if let Some((&ann, &prob)) = best
        && prob >= confidence
    {
        confident.insert(ann, row_get(row, ann).unwrap().alignment_idx);
    }

    ResolutionResult {
        confident,
        probabilities,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::grouping::{AlignmentEntry, ReadGroup};
    use crate::intern::Interner;
    use crate::overlap::AnnotationOverlap;

    fn make_entry(
        input_order: usize,
        as_score: Option<i64>,
        proper_pair: bool,
        annotations: Vec<(Symbol, usize)>,
    ) -> AlignmentEntry {
        AlignmentEntry {
            qname: Some("test_read".to_string()),
            input_order,
            mate_input_order: None,
            annotations: annotations
                .into_iter()
                .map(|(sym, bp)| AnnotationOverlap {
                    annotation: sym,
                    overlap_bp: bp,
                })
                .collect(),
            alignment_score: as_score,
            is_proper_pair: proper_pair,
        }
    }

    /// Test helper: compute priors + confidence pick without EM.
    /// Matches the behavior of the pre-EM `resolve_ambiguity` and keeps the
    /// ribbonfish/telescope tests readable.
    fn resolve_no_em(
        group: &ReadGroup,
        row: &MatrixRow,
        confidence: f64,
        model: ModelType,
        global_as_min: Option<i64>,
        global_as_max: Option<i64>,
    ) -> ResolutionResult {
        let fp = compute_priors(row, model, global_as_min, global_as_max);
        let fposts = FragmentPosteriors {
            annotations: fp.annotations,
            probs: fp.priors,
        };
        select_confident(group, row, &fposts, confidence)
    }

    fn make_group(alignments: Vec<AlignmentEntry>) -> ReadGroup {
        ReadGroup {
            qname: "test_read".to_string(),
            alignments,
        }
    }

    fn no_feat(int: &mut Interner) -> Symbol {
        int.intern(crate::overlap::NO_FEATURE)
    }

    #[test]
    fn test_as_score_wins() {
        let mut int = Interner::new();
        let gene_a = int.intern("GENE_A");
        let group = make_group(vec![
            make_entry(0, Some(100), false, vec![(gene_a, 50)]),
            make_entry(1, Some(200), false, vec![(gene_a, 50)]),
        ]);
        let row = select_representatives(&group, no_feat(&mut int));
        assert_eq!(row_get(&row, gene_a).unwrap().alignment_idx, 1); // higher AS
        assert_eq!(row_get(&row, gene_a).unwrap().alignment_score, 200);
    }

    #[test]
    fn test_proper_pair_breaks_as_tie() {
        let mut int = Interner::new();
        let gene_a = int.intern("GENE_A");
        let group = make_group(vec![
            make_entry(0, Some(100), false, vec![(gene_a, 50)]),
            make_entry(1, Some(100), true, vec![(gene_a, 50)]),
        ]);
        let row = select_representatives(&group, no_feat(&mut int));
        assert_eq!(row_get(&row, gene_a).unwrap().alignment_idx, 1); // proper pair
    }

    #[test]
    fn test_overlap_breaks_tie() {
        let mut int = Interner::new();
        let gene_a = int.intern("GENE_A");
        let group = make_group(vec![
            make_entry(0, Some(100), true, vec![(gene_a, 30)]),
            make_entry(1, Some(100), true, vec![(gene_a, 80)]),
        ]);
        let row = select_representatives(&group, no_feat(&mut int));
        assert_eq!(row_get(&row, gene_a).unwrap().alignment_idx, 1); // more overlap
    }

    #[test]
    fn test_input_order_breaks_tie() {
        let mut int = Interner::new();
        let gene_a = int.intern("GENE_A");
        let group = make_group(vec![
            make_entry(5, Some(100), true, vec![(gene_a, 50)]),
            make_entry(3, Some(100), true, vec![(gene_a, 50)]),
        ]);
        let row = select_representatives(&group, no_feat(&mut int));
        assert_eq!(row_get(&row, gene_a).unwrap().alignment_idx, 1); // input_order 3 < 5
    }

    #[test]
    fn test_multi_annotation_independent() {
        let mut int = Interner::new();
        let gene_a = int.intern("GENE_A");
        let gene_b = int.intern("GENE_B");
        let group = make_group(vec![
            make_entry(0, Some(200), false, vec![(gene_a, 80), (gene_b, 20)]),
            make_entry(1, Some(100), false, vec![(gene_a, 30), (gene_b, 90)]),
        ]);
        let row = select_representatives(&group, no_feat(&mut int));
        // GENE_A: entry 0 wins (higher AS)
        assert_eq!(row_get(&row, gene_a).unwrap().alignment_idx, 0);
        // GENE_B: entry 0 wins (higher AS, even though less overlap)
        assert_eq!(row_get(&row, gene_b).unwrap().alignment_idx, 0);
    }

    #[test]
    fn test_no_as_tag_loses() {
        let mut int = Interner::new();
        let gene_a = int.intern("GENE_A");
        let group = make_group(vec![
            make_entry(0, None, false, vec![(gene_a, 50)]),
            make_entry(1, Some(1), false, vec![(gene_a, 50)]),
        ]);
        let row = select_representatives(&group, no_feat(&mut int));
        assert_eq!(row_get(&row, gene_a).unwrap().alignment_idx, 1); // None → i64::MIN
    }

    // ========== priors + select_confident tests (via resolve_no_em) ==========

    #[test]
    fn test_resolve_single_annotation_always_confident() {
        let mut int = Interner::new();
        let gene_a = int.intern("GENE_A");
        let group = make_group(vec![make_entry(0, Some(100), false, vec![(gene_a, 50)])]);
        let row = select_representatives(&group, no_feat(&mut int));
        let res = resolve_no_em(&group, &row, 0.5, ModelType::Ribbonfish, None, None);
        assert_eq!(res.probabilities[&gene_a], 1.0);
        assert!(res.confident.contains_key(&gene_a));

        // Even at threshold 1.0
        let res = resolve_no_em(&group, &row, 1.0, ModelType::Ribbonfish, None, None);
        assert!(res.confident.contains_key(&gene_a));
    }

    #[test]
    fn test_resolve_equal_scores_threshold_half() {
        let mut int = Interner::new();
        let gene_a = int.intern("GENE_A");
        let gene_b = int.intern("GENE_B");
        // Two annotations with equal AS → each gets probability 0.5
        // Only ONE should be confident (single-winner with tiebreak)
        let group = make_group(vec![
            make_entry(0, Some(100), false, vec![(gene_a, 50)]),
            make_entry(1, Some(100), false, vec![(gene_b, 50)]),
        ]);
        let row = select_representatives(&group, no_feat(&mut int));
        let res = resolve_no_em(&group, &row, 0.5, ModelType::Ribbonfish, None, None);

        assert!((res.probabilities[&gene_a] - 0.5).abs() < 1e-10);
        assert!((res.probabilities[&gene_b] - 0.5).abs() < 1e-10);
        // Only one wins via tiebreak (input_order 0 < 1)
        assert_eq!(res.confident.len(), 1);
    }

    #[test]
    fn test_resolve_equal_scores_threshold_high() {
        let mut int = Interner::new();
        let gene_a = int.intern("GENE_A");
        let gene_b = int.intern("GENE_B");
        let group = make_group(vec![
            make_entry(0, Some(100), false, vec![(gene_a, 50)]),
            make_entry(1, Some(100), false, vec![(gene_b, 50)]),
        ]);
        let row = select_representatives(&group, no_feat(&mut int));
        // At threshold 0.6, neither 0.5 passes
        let res = resolve_no_em(&group, &row, 0.6, ModelType::Ribbonfish, None, None);
        assert!(res.confident.is_empty());
    }

    #[test]
    fn test_resolve_different_scores() {
        let mut int = Interner::new();
        let gene_a = int.intern("GENE_A");
        let gene_b = int.intern("GENE_B");
        // GENE_A has much higher AS → high probability
        let group = make_group(vec![
            make_entry(0, Some(200), false, vec![(gene_a, 50)]),
            make_entry(1, Some(0), false, vec![(gene_b, 50)]),
        ]);
        let row = select_representatives(&group, no_feat(&mut int));
        let res = resolve_no_em(&group, &row, 0.9, ModelType::Ribbonfish, None, None);

        assert!(res.probabilities[&gene_a] > 0.9);
        assert!(res.probabilities[&gene_b] < 0.1);
        assert!(res.confident.contains_key(&gene_a));
        assert!(!res.confident.contains_key(&gene_b));
    }

    #[test]
    fn test_resolve_unique_mode() {
        let mut int = Interner::new();
        let gene_a = int.intern("GENE_A");
        let gene_b = int.intern("GENE_B");
        // threshold=1.0: only single-annotation groups pass
        // Use close AS scores so neither probability rounds to exactly 1.0
        let group = make_group(vec![
            make_entry(0, Some(10), false, vec![(gene_a, 50)]),
            make_entry(1, Some(0), false, vec![(gene_b, 50)]),
        ]);
        let row = select_representatives(&group, no_feat(&mut int));
        let res = resolve_no_em(&group, &row, 1.0, ModelType::Ribbonfish, None, None);
        assert!(res.confident.is_empty()); // multi-annotation → nothing passes

        // Single annotation should pass at threshold 1.0
        let group2 = make_group(vec![make_entry(0, Some(10), false, vec![(gene_a, 50)])]);
        let row2 = select_representatives(&group2, no_feat(&mut int));
        let res2 = resolve_no_em(&group2, &row2, 1.0, ModelType::Ribbonfish, None, None);
        assert!(res2.confident.contains_key(&gene_a));
    }

    #[test]
    fn test_resolve_no_as_tag_uses_zero() {
        let mut int = Interner::new();
        let gene_a = int.intern("GENE_A");
        let gene_b = int.intern("GENE_B");
        // Both have no AS → both get score 0 → equal probability 0.5
        // Only one wins via tiebreak
        let group = make_group(vec![
            make_entry(0, None, false, vec![(gene_a, 50)]),
            make_entry(1, None, false, vec![(gene_b, 50)]),
        ]);
        let row = select_representatives(&group, no_feat(&mut int));
        let res = resolve_no_em(&group, &row, 0.5, ModelType::Ribbonfish, None, None);
        assert!((res.probabilities[&gene_a] - 0.5).abs() < 1e-10);
        assert!((res.probabilities[&gene_b] - 0.5).abs() < 1e-10);
        assert_eq!(res.confident.len(), 1);
    }

    #[test]
    fn test_resolve_no_overlap_yields_no_feature() {
        // Every alignment misses → row contains only __no_feature__ →
        // confidently picked at any threshold (single-annotation row).
        let mut int = Interner::new();
        let no_feature = no_feat(&mut int);
        let group = make_group(vec![make_entry(0, None, false, vec![])]);
        let row = select_representatives(&group, no_feature);
        let res = resolve_no_em(&group, &row, 0.5, ModelType::Ribbonfish, None, None);
        assert_eq!(row.len(), 1);
        assert!(row_get(&row, no_feature).is_some());
        assert!((res.probabilities[&no_feature] - 1.0).abs() < 1e-10);
        assert!(res.confident.contains_key(&no_feature));
    }

    #[test]
    fn test_resolve_telescope_linear_rescale() {
        let mut int = Interner::new();
        let gene_a = int.intern("GENE_A");
        let gene_b = int.intern("GENE_B");
        // AS 80 and 100; global bounds 50..=150 → scaled 30 and 50.
        let group = make_group(vec![
            make_entry(0, Some(80), false, vec![(gene_a, 50)]),
            make_entry(1, Some(100), false, vec![(gene_b, 50)]),
        ]);
        let row = select_representatives(&group, no_feat(&mut int));
        let res = resolve_no_em(&group, &row, 0.5, ModelType::Telescope, Some(50), Some(150));

        let pa = res.probabilities[&gene_a];
        let pb = res.probabilities[&gene_b];
        assert!(pb > pa, "higher AS should get higher probability");
        assert!((pa + pb - 1.0).abs() < 1e-10);

        // Expected: weights = expm1(30), expm1(50); pa = expm1(30) / (expm1(30)+expm1(50))
        let w_a = 30f64.exp_m1();
        let w_b = 50f64.exp_m1();
        let expected_pa = w_a / (w_a + w_b);
        let expected_pb = w_b / (w_a + w_b);
        assert!((pa - expected_pa).abs() < 1e-10);
        assert!((pb - expected_pb).abs() < 1e-10);
    }

    #[test]
    fn test_resolve_telescope_global_min_contributes_zero() {
        let mut int = Interner::new();
        let gene_a = int.intern("GENE_A");
        let gene_b = int.intern("GENE_B");
        // gene_a is at the global min → its expm1(0) = 0 → prob 0.
        let group = make_group(vec![
            make_entry(0, Some(50), false, vec![(gene_a, 50)]),
            make_entry(1, Some(150), false, vec![(gene_b, 50)]),
        ]);
        let row = select_representatives(&group, no_feat(&mut int));
        let res = resolve_no_em(&group, &row, 0.9, ModelType::Telescope, Some(50), Some(150));
        assert_eq!(res.probabilities[&gene_a], 0.0);
        assert!((res.probabilities[&gene_b] - 1.0).abs() < 1e-10);
        assert!(res.confident.contains_key(&gene_b));
    }

    #[test]
    fn test_resolve_telescope_degenerate_uniform() {
        let mut int = Interner::new();
        let gene_a = int.intern("GENE_A");
        let gene_b = int.intern("GENE_B");
        let group = make_group(vec![
            make_entry(0, Some(100), false, vec![(gene_a, 50)]),
            make_entry(1, Some(100), false, vec![(gene_b, 50)]),
        ]);
        let row = select_representatives(&group, no_feat(&mut int));
        // global_min == global_max: no range, fall back to uniform.
        let res = resolve_no_em(
            &group,
            &row,
            0.5,
            ModelType::Telescope,
            Some(100),
            Some(100),
        );
        assert!((res.probabilities[&gene_a] - 0.5).abs() < 1e-10);
        assert!((res.probabilities[&gene_b] - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_resolve_telescope_no_global_bounds_uniform() {
        let mut int = Interner::new();
        let gene_a = int.intern("GENE_A");
        let gene_b = int.intern("GENE_B");
        let group = make_group(vec![
            make_entry(0, Some(80), false, vec![(gene_a, 50)]),
            make_entry(1, Some(100), false, vec![(gene_b, 50)]),
        ]);
        let row = select_representatives(&group, no_feat(&mut int));
        // No global bounds available → uniform fallback.
        let res = resolve_no_em(&group, &row, 0.5, ModelType::Telescope, None, None);
        assert!((res.probabilities[&gene_a] - 0.5).abs() < 1e-10);
        assert!((res.probabilities[&gene_b] - 0.5).abs() < 1e-10);
    }

    // ========== expectation_maximization tests ==========

    #[test]
    fn test_em_zero_iterations_returns_priors() {
        let fp = FragmentPriors {
            annotations: vec![Symbol_from(0), Symbol_from(1)],
            priors: vec![0.7, 0.3],
        };
        let priors_copy = fp.priors.clone();
        let out = expectation_maximization(vec![fp], 2, 0, 0.0, &[], false).posteriors;
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].probs, priors_copy);
    }

    #[test]
    fn test_em_single_annotation_stays_one() {
        // A fragment with a single annotation has prior [1.0]; γ must remain [1.0]
        // regardless of how many EM iterations run.
        let fp = FragmentPriors {
            annotations: vec![Symbol_from(0)],
            priors: vec![1.0],
        };
        let out = expectation_maximization(vec![fp], 1, 25, 0.0, &[], false).posteriors;
        assert_eq!(out.len(), 1);
        assert!((out[0].probs[0] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_em_two_fragments_reinforce() {
        // Both fragments slightly favor X over Y. EM should amplify X's share.
        let frag_a = FragmentPriors {
            annotations: vec![Symbol_from(0), Symbol_from(1)],
            priors: vec![0.6, 0.4],
        };
        let frag_b = FragmentPriors {
            annotations: vec![Symbol_from(0), Symbol_from(1)],
            priors: vec![0.55, 0.45],
        };
        let frag_a_priors = frag_a.priors.clone();
        let frag_b_priors = frag_b.priors.clone();
        let out = expectation_maximization(vec![frag_a, frag_b], 2, 10, 0.0, &[], false).posteriors;
        assert!(
            out[0].probs[0] > frag_a_priors[0],
            "fragment A's γ on X ({}) should exceed its prior ({})",
            out[0].probs[0],
            frag_a_priors[0]
        );
        assert!(
            out[1].probs[0] > frag_b_priors[0],
            "fragment B's γ on X ({}) should exceed its prior ({})",
            out[1].probs[0],
            frag_b_priors[0]
        );
        assert!((out[0].probs.iter().sum::<f64>() - 1.0).abs() < 1e-10);
        assert!((out[1].probs.iter().sum::<f64>() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_partial_overlap_yields_real_and_no_feature() {
        // One alignment overlaps GENE_A, one alignment misses → row contains both.
        let mut int = Interner::new();
        let no_feature = no_feat(&mut int);
        let gene_a = int.intern("GENE_A");
        let group = make_group(vec![
            make_entry(0, Some(100), false, vec![(gene_a, 50)]),
            make_entry(1, Some(150), false, vec![]),
        ]);
        let row = select_representatives(&group, no_feature);
        assert_eq!(row.len(), 2);
        assert!(row_get(&row, gene_a).is_some());
        assert!(row_get(&row, no_feature).is_some());
        // The synthetic candidate's representative is the higher-AS no-overlap alignment.
        assert_eq!(row_get(&row, no_feature).unwrap().alignment_idx, 1);
        assert_eq!(row_get(&row, no_feature).unwrap().alignment_score, 150);
    }

    #[test]
    fn test_em_no_feature_competes_with_real_annotations() {
        // Two fragments each have a weak GENE_X prior and a stronger __no_feature__ prior.
        // EM should reinforce __no_feature__ (shared across both fragments) over GENE_X.
        let mut int = Interner::new();
        let no_feature = no_feat(&mut int);
        let gene_x = int.intern("GENE_X");
        let frag_a = FragmentPriors {
            annotations: vec![gene_x, no_feature],
            priors: vec![0.4, 0.6],
        };
        let frag_b = FragmentPriors {
            annotations: vec![gene_x, no_feature],
            priors: vec![0.3, 0.7],
        };
        let frag_a_priors = frag_a.priors.clone();
        let frag_b_priors = frag_b.priors.clone();
        let out = expectation_maximization(vec![frag_a, frag_b], int.len(), 10, 0.0, &[], false).posteriors;
        assert!(
            out[0].probs[1] > frag_a_priors[1],
            "fragment A's γ on __no_feature__ ({}) should exceed its prior ({})",
            out[0].probs[1],
            frag_a_priors[1]
        );
        assert!(
            out[1].probs[1] > frag_b_priors[1],
            "fragment B's γ on __no_feature__ ({}) should exceed its prior ({})",
            out[1].probs[1],
            frag_b_priors[1]
        );
    }

    #[test]
    fn test_em_empty_fragment_unchanged() {
        let frag_a = FragmentPriors {
            annotations: vec![Symbol_from(0), Symbol_from(1)],
            priors: vec![0.5, 0.5],
        };
        let frag_empty = FragmentPriors {
            annotations: vec![],
            priors: vec![],
        };
        let out = expectation_maximization(vec![frag_a, frag_empty], 2, 5, 0.0, &[], false).posteriors;
        assert_eq!(out[1].probs.len(), 0); // empty stays empty, no panic
        assert!((out[0].probs.iter().sum::<f64>() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_em_length_correction_favors_shorter() {
        // Equal priors, equal seed θ. After exactly one EM iteration with
        // length correction on, γ should be in inverse proportion to length:
        // L_A = 100, L_B = 200 → γ_A / γ_B = 2.
        let mut int = Interner::new();
        let a = int.intern("A");
        let b = int.intern("B");
        let lengths = vec![100u64, 200u64];

        let fp = FragmentPriors {
            annotations: vec![a, b],
            priors: vec![0.5, 0.5],
        };
        let out = expectation_maximization(vec![fp], int.len(), 1, 0.0, &lengths, true).posteriors;
        let probs = &out[0].probs;
        let ratio = probs[0] / probs[1];
        assert!(
            (ratio - 2.0).abs() < 1e-6,
            "expected γ_A / γ_B == 2 after 1 iter, got ratio {ratio} (probs {probs:?})",
        );
        // probs sum to 1.
        assert!((probs.iter().sum::<f64>() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_em_length_correction_off_matches_legacy() {
        // Same fixture, length_correction = false → uniform priors stay uniform
        // through any number of iterations (no length skew).
        let mut int = Interner::new();
        let a = int.intern("A");
        let b = int.intern("B");
        let lengths = vec![100u64, 200u64];
        let fp = FragmentPriors {
            annotations: vec![a, b],
            priors: vec![0.5, 0.5],
        };
        let out = expectation_maximization(vec![fp], int.len(), 10, 0.0, &lengths, false).posteriors;
        let probs = &out[0].probs;
        assert!((probs[0] - 0.5).abs() < 1e-10);
        assert!((probs[1] - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_em_length_correction_zero_length_falls_back() {
        // Defensive: a length-0 entry must not divide-by-zero. With both
        // lengths zero, length_correction is effectively off.
        let mut int = Interner::new();
        let a = int.intern("A");
        let b = int.intern("B");
        let lengths = vec![0u64, 0u64];
        let fp = FragmentPriors {
            annotations: vec![a, b],
            priors: vec![0.6, 0.4],
        };
        let out = expectation_maximization(vec![fp], int.len(), 5, 0.0, &lengths, true).posteriors;
        // With inv_len=1.0 for both, this collapses to the legacy update.
        assert!(out[0].probs.iter().all(|p| p.is_finite()));
        assert!((out[0].probs.iter().sum::<f64>() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_em_early_stop_reports_iterations() {
        // Two fragments → total assigned mass is only 2, so no per-annotation
        // mass can change by a whole fragment in a single step. With the default
        // threshold of 1.0 EM converges on the very first iteration and reports
        // it, well short of the 100 cap.
        let frag_a = FragmentPriors {
            annotations: vec![Symbol_from(0), Symbol_from(1)],
            priors: vec![0.6, 0.4],
        };
        let frag_b = FragmentPriors {
            annotations: vec![Symbol_from(0), Symbol_from(1)],
            priors: vec![0.55, 0.45],
        };
        let outcome = expectation_maximization(vec![frag_a, frag_b], 2, 100, 1.0, &[], false);
        assert_eq!(outcome.iterations, 1);
        assert!(outcome.iterations < 100);
    }

    #[test]
    fn test_em_min_mass_delta_zero_runs_full_max_iter() {
        // threshold 0 → delta < 0 never trips → all max_iter iterations run.
        let frag_a = FragmentPriors {
            annotations: vec![Symbol_from(0), Symbol_from(1)],
            priors: vec![0.6, 0.4],
        };
        let frag_b = FragmentPriors {
            annotations: vec![Symbol_from(0), Symbol_from(1)],
            priors: vec![0.55, 0.45],
        };
        let outcome = expectation_maximization(vec![frag_a, frag_b], 2, 7, 0.0, &[], false);
        assert_eq!(outcome.iterations, 7);
    }

    #[test]
    fn test_em_zero_max_iter_reports_zero_iterations() {
        let fp = FragmentPriors {
            annotations: vec![Symbol_from(0), Symbol_from(1)],
            priors: vec![0.7, 0.3],
        };
        let outcome = expectation_maximization(vec![fp], 2, 0, 1.0, &[], false);
        assert_eq!(outcome.iterations, 0);
        assert_eq!(outcome.posteriors[0].probs, vec![0.7, 0.3]);
    }

    // ========== fragment_weight tests ==========

    fn cell(score: i64) -> MatrixCell {
        MatrixCell {
            alignment_idx: 0,
            alignment_score: score,
        }
    }

    fn weight_row(scores: &[i64]) -> Vec<(Symbol, MatrixCell)> {
        scores
            .iter()
            .enumerate()
            .map(|(i, &s)| (Symbol_from(i as u32), cell(s)))
            .collect()
    }

    #[test]
    fn test_fragment_weight_max_of_transformed() {
        let row = weight_row(&[80, 100]);
        let w = fragment_weight(&row, Some(50), Some(150));
        let expected = ((100.0_f64 - 50.0) / 100.0 * 100.0).exp_m1();
        assert_eq!(w, expected);
    }

    #[test]
    fn test_fragment_weight_uses_max_not_sum() {
        // Multiple alignments must not reduce or inflate the weight: a row with
        // [80,100,90] has the same weight as one with just [100].
        let multi = fragment_weight(&weight_row(&[80, 100, 90]), Some(50), Some(150));
        let single = fragment_weight(&weight_row(&[100]), Some(50), Some(150));
        assert_eq!(multi, single);
    }

    #[test]
    fn test_fragment_weight_no_bounds_is_one() {
        assert_eq!(fragment_weight(&weight_row(&[80, 100]), None, None), 1.0);
    }

    #[test]
    fn test_fragment_weight_degenerate_bounds_is_one() {
        // hi == lo → no range → unweighted (1.0).
        assert_eq!(
            fragment_weight(&weight_row(&[100, 100]), Some(100), Some(100)),
            1.0
        );
    }

    #[test]
    fn test_fragment_weight_global_min_is_zero() {
        // Best score sits at the global minimum → expm1(0) == 0.
        assert_eq!(
            fragment_weight(&weight_row(&[50]), Some(50), Some(150)),
            0.0
        );
    }

    // ========== telescope_em tests ==========

    #[test]
    fn test_telescope_em_zero_iterations_returns_priors() {
        let amb = FragmentPriors {
            annotations: vec![Symbol_from(0), Symbol_from(1)],
            priors: vec![0.6, 0.4],
        };
        let uniq = FragmentPriors {
            annotations: vec![Symbol_from(0)],
            priors: vec![1.0],
        };
        let out = telescope_em(vec![amb, uniq], &[1.0, 1.0], 2, 0, 0.0, 0.0, Some(0), Some(100), &[], false).posteriors;
        assert_eq!(out.len(), 2);
        assert_eq!(out[0].probs, vec![0.6, 0.4]);
        assert_eq!(out[1].probs, vec![1.0]);
    }

    #[test]
    fn test_telescope_em_unique_fragment_stays_one() {
        // A unique fragment's γ must stay [1.0] even while ambiguous EM runs.
        let amb = FragmentPriors {
            annotations: vec![Symbol_from(0), Symbol_from(1)],
            priors: vec![0.5, 0.5],
        };
        let uniq = FragmentPriors {
            annotations: vec![Symbol_from(0)],
            priors: vec![1.0],
        };
        let out = telescope_em(vec![amb, uniq], &[1.0, 1.0], 2, 25, 0.0, 0.0, Some(0), Some(100), &[], false).posteriors;
        assert_eq!(out[1].probs, vec![1.0]);
    }

    #[test]
    fn test_telescope_em_uniques_boost_pi() {
        // Heavy unique mass on A lifts π_A, so an equal-prior ambiguous {A,B}
        // fragment is pulled onto A.
        let amb = FragmentPriors {
            annotations: vec![Symbol_from(0), Symbol_from(1)],
            priors: vec![0.5, 0.5],
        };
        let uniq = FragmentPriors {
            annotations: vec![Symbol_from(0)],
            priors: vec![1.0],
        };
        let out = telescope_em(
            vec![amb, uniq],
            &[1.0, 10.0],
            2,
            20,
            0.0,
            0.0,
            Some(0),
            Some(100),
            &[],
            false,
        )
        .posteriors;
        assert!(
            out[0].probs[0] > out[0].probs[1],
            "ambiguous γ should favor A: {:?}",
            out[0].probs
        );
        assert!(out[0].probs[0] > 0.5);
    }

    #[test]
    fn test_telescope_em_weight_scales_mstep() {
        // f1 (weight 10) favors A, f2 (weight 1) favors B. The heavy fragment
        // dominates the shared θ, dragging the light fragment's γ toward A
        // (above its own prior of 0.4 on A).
        let f1 = FragmentPriors {
            annotations: vec![Symbol_from(0), Symbol_from(1)],
            priors: vec![0.6, 0.4],
        };
        let f2 = FragmentPriors {
            annotations: vec![Symbol_from(0), Symbol_from(1)],
            priors: vec![0.4, 0.6],
        };
        let out = telescope_em(vec![f1, f2], &[10.0, 1.0], 2, 20, 0.0, 0.0, Some(0), Some(100), &[], false).posteriors;
        assert!(
            out[1].probs[0] > 0.4,
            "light B-favoring fragment should be dragged toward A: {:?}",
            out[1].probs
        );
        assert!(out[0].probs[0] > 0.6);
    }

    #[test]
    fn test_telescope_em_theta_prior_regularizes() {
        // Two fragments both favoring A. A huge theta_prior flattens θ so only
        // π (unregularized) reinforces → less collapse than theta_prior == 0.
        let mk = || {
            vec![
                FragmentPriors {
                    annotations: vec![Symbol_from(0), Symbol_from(1)],
                    priors: vec![0.7, 0.3],
                },
                FragmentPriors {
                    annotations: vec![Symbol_from(0), Symbol_from(1)],
                    priors: vec![0.7, 0.3],
                },
            ]
        };
        let no_reg = telescope_em(mk(), &[1.0, 1.0], 2, 15, 0.0, 0.0, Some(0), Some(100), &[], false).posteriors;
        let reg = telescope_em(mk(), &[1.0, 1.0], 2, 15, 0.0, 1e7, Some(0), Some(100), &[], false).posteriors;
        assert!(
            reg[0].probs[0] + 1e-6 < no_reg[0].probs[0],
            "regularized γ_A ({}) should be below unregularized ({})",
            reg[0].probs[0],
            no_reg[0].probs[0],
        );
    }

    #[test]
    fn test_telescope_em_e_step_multiplies_pi_and_theta() {
        // Closed form for one iteration: ambiguous {A,B} priors [0.5,0.5],
        // one unique on A (weight 1) → π_fixed=[1,0]. amb_seed=[0.5,0.5];
        // θ⁰=[0.5,0.5]; π⁰=normalize([1.5,0.5])=[0.75,0.25].
        // γ_A ∝ 0.5·0.75·0.5, γ_B ∝ 0.5·0.25·0.5 → γ=[0.75,0.25].
        let amb = FragmentPriors {
            annotations: vec![Symbol_from(0), Symbol_from(1)],
            priors: vec![0.5, 0.5],
        };
        let uniq = FragmentPriors {
            annotations: vec![Symbol_from(0)],
            priors: vec![1.0],
        };
        let out = telescope_em(vec![amb, uniq], &[1.0, 1.0], 2, 1, 0.0, 0.0, Some(0), Some(100), &[], false).posteriors;
        assert!((out[0].probs[0] - 0.75).abs() < 1e-12, "{:?}", out[0].probs);
        assert!((out[0].probs[1] - 0.25).abs() < 1e-12, "{:?}", out[0].probs);
        assert_eq!(out[1].probs, vec![1.0]);
    }

    #[test]
    fn test_telescope_em_length_correction_favors_shorter() {
        // Same setup as the E-step test above, but with length correction on
        // and A twice as long as B (lengths [2,1]). The E-step substitutes
        // θ_a/L_a, so for one iteration: θ⁰=[0.5,0.5], π⁰=[0.75,0.25].
        // γ_A ∝ 0.5·0.75·(0.5/2) = 0.09375, γ_B ∝ 0.5·0.25·(0.5/1) = 0.0625;
        // sum 0.15625 → γ=[0.6,0.4]: mass shifts toward the shorter B versus
        // the uncorrected [0.75,0.25].
        let amb = FragmentPriors {
            annotations: vec![Symbol_from(0), Symbol_from(1)],
            priors: vec![0.5, 0.5],
        };
        let uniq = FragmentPriors {
            annotations: vec![Symbol_from(0)],
            priors: vec![1.0],
        };
        let out = telescope_em(
            vec![amb, uniq],
            &[1.0, 1.0],
            2,
            1,
            0.0,
            0.0,
            Some(0),
            Some(100),
            &[2, 1],
            true,
        )
        .posteriors;
        assert!((out[0].probs[0] - 0.6).abs() < 1e-12, "{:?}", out[0].probs);
        assert!((out[0].probs[1] - 0.4).abs() < 1e-12, "{:?}", out[0].probs);
        assert_eq!(out[1].probs, vec![1.0]);
    }

    #[test]
    fn test_telescope_em_degenerate_uniform() {
        // Zero weight, no uniques → π collapses to all-zero → every prior·π·θ
        // is 0 → uniform fallback regardless of the prior.
        let amb = FragmentPriors {
            annotations: vec![Symbol_from(0), Symbol_from(1)],
            priors: vec![0.9, 0.1],
        };
        let out = telescope_em(vec![amb], &[0.0], 2, 5, 0.0, 0.0, Some(0), Some(100), &[], false).posteriors;
        assert!((out[0].probs[0] - 0.5).abs() < 1e-12, "{:?}", out[0].probs);
        assert!((out[0].probs[1] - 0.5).abs() < 1e-12, "{:?}", out[0].probs);
    }

    #[test]
    fn test_telescope_em_empty_fragment_unchanged() {
        let amb = FragmentPriors {
            annotations: vec![Symbol_from(0), Symbol_from(1)],
            priors: vec![0.5, 0.5],
        };
        let empty = FragmentPriors {
            annotations: vec![],
            priors: vec![],
        };
        let out = telescope_em(vec![amb, empty], &[1.0, 1.0], 2, 5, 0.0, 0.0, Some(0), Some(100), &[], false).posteriors;
        assert_eq!(out[1].probs.len(), 0);
        assert!((out[0].probs.iter().sum::<f64>() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_telescope_em_all_unique_short_circuits() {
        // No ambiguous fragments → priors returned verbatim, no panic.
        let u0 = FragmentPriors {
            annotations: vec![Symbol_from(0)],
            priors: vec![1.0],
        };
        let u1 = FragmentPriors {
            annotations: vec![Symbol_from(1)],
            priors: vec![1.0],
        };
        let out = telescope_em(
            vec![u0, u1],
            &[1.0, 1.0],
            2,
            10,
            0.0,
            200_000.0,
            Some(0),
            Some(100),
            &[],
            false,
        )
        .posteriors;
        assert_eq!(out.len(), 2);
        assert_eq!(out[0].probs, vec![1.0]);
        assert_eq!(out[1].probs, vec![1.0]);
    }

    #[test]
    fn test_telescope_em_output_index_aligned() {
        // Order: [ambiguous, unique, empty, ambiguous]. Output must stay
        // index-aligned with input (main.rs maps posteriors[i] -> groups[i]).
        let priors = vec![
            FragmentPriors {
                annotations: vec![Symbol_from(0), Symbol_from(1)],
                priors: vec![0.6, 0.4],
            },
            FragmentPriors {
                annotations: vec![Symbol_from(0)],
                priors: vec![1.0],
            },
            FragmentPriors {
                annotations: vec![],
                priors: vec![],
            },
            FragmentPriors {
                annotations: vec![Symbol_from(0), Symbol_from(1)],
                priors: vec![0.3, 0.7],
            },
        ];
        let out = telescope_em(priors, &[1.0, 1.0, 1.0, 1.0], 2, 3, 0.0, 0.0, Some(0), Some(100), &[], false).posteriors;
        assert_eq!(out.len(), 4);
        assert_eq!(out[0].annotations.len(), 2);
        assert_eq!(out[1].annotations.len(), 1);
        assert_eq!(out[1].probs, vec![1.0]);
        assert_eq!(out[2].annotations.len(), 0);
        assert_eq!(out[3].annotations.len(), 2);
        assert!((out[0].probs.iter().sum::<f64>() - 1.0).abs() < 1e-10);
        assert!((out[3].probs.iter().sum::<f64>() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_telescope_em_early_stop_reports_iterations() {
        // Small ambiguous + unique set: total assigned mass is a few fragments,
        // so the per-annotation delta drops below 1 within the first iteration.
        let amb = FragmentPriors {
            annotations: vec![Symbol_from(0), Symbol_from(1)],
            priors: vec![0.6, 0.4],
        };
        let uniq = FragmentPriors {
            annotations: vec![Symbol_from(0)],
            priors: vec![1.0],
        };
        let outcome = telescope_em(
            vec![amb, uniq],
            &[1.0, 1.0],
            2,
            100,
            1.0,
            0.0,
            Some(0),
            Some(100),
            &[],
            false,
        );
        assert!(outcome.iterations >= 1 && outcome.iterations < 100);

        // No ambiguous fragments → EM never iterates → 0 even with max_iter > 0.
        let only_unique = FragmentPriors {
            annotations: vec![Symbol_from(0)],
            priors: vec![1.0],
        };
        let none = telescope_em(
            vec![only_unique],
            &[1.0],
            1,
            100,
            1.0,
            0.0,
            Some(0),
            Some(100),
            &[],
            false,
        );
        assert_eq!(none.iterations, 0);
    }

    #[test]
    fn test_max_fragment_weight() {
        // Usable bounds → a perfectly-aligned fragment scales to 100.
        assert_eq!(max_fragment_weight(Some(0), Some(100)), 100.0_f64.exp_m1());
        // Degenerate: hi == lo, hi < lo, or missing bounds → 1.0, matching
        // fragment_weight's own fallback so genuine reads are all 1.0 too.
        assert_eq!(max_fragment_weight(Some(50), Some(50)), 1.0);
        assert_eq!(max_fragment_weight(Some(150), Some(50)), 1.0);
        assert_eq!(max_fragment_weight(None, None), 1.0);
        assert_eq!(max_fragment_weight(Some(0), None), 1.0);
    }

    #[test]
    fn test_telescope_em_degenerate_regularizer_weighted_like_genuine_reads() {
        // In the degenerate case (max AS == min AS) every genuine fragment's
        // weight falls back to 1.0, so each theta_prior regularizing read must
        // also be weighted 1.0 — NOT expm1(100). Pinned via a closed-form
        // 2-iteration computation: single ambiguous {A,B}, no uniques, prior
        // [0.6,0.4], weight w=2, theta_prior T=1, degenerate bounds.
        //
        //   θ⁰=π⁰=[0.5,0.5]; iter1 leaves γ=prior; M-step (reg = T·1):
        //     θ = [(wp+T),(w(1-p)+T)] / (w+2T) = [2.2,1.8]/4 = [0.55,0.45]
        //     π = [wp, w(1-p)] / w               = [0.6,0.4]
        //   iter2 E-step: γ_A ∝ p²·θ_A = 0.36·0.55, γ_B ∝ 0.4²·θ_B = 0.16·0.45
        //     → γ_A = 0.198 / 0.270 = 0.733333…
        //
        // If the regularizer used T·expm1(100) instead, θ would be ≈[0.5,0.5]
        // and γ_A would be ≈0.6923 — this exact assertion would then fail.
        let amb = FragmentPriors {
            annotations: vec![Symbol_from(0), Symbol_from(1)],
            priors: vec![0.6, 0.4],
        };
        let out = telescope_em(vec![amb], &[2.0], 2, 2, 0.0, 1.0, None, None, &[], false).posteriors;
        assert!(
            (out[0].probs[0] - 0.7333333333333333).abs() < 1e-12,
            "degenerate regularizer must use weight 1.0, got γ_A={:?}",
            out[0].probs,
        );
        assert!((out[0].probs[1] - 0.26666666666666666).abs() < 1e-12);
    }

    /// Test-only: build a `Symbol` from a raw u32 index without going through
    /// an `Interner` (the EM tests don't care about string resolution).
    #[allow(non_snake_case)]
    fn Symbol_from(i: u32) -> Symbol {
        let mut int = Interner::new();
        // Intern `i+1` placeholders so the returned symbol has index `i`.
        for k in 0..=i {
            int.intern(&format!("__test_sym_{k}"));
        }
        int.intern(&format!("__test_sym_{i}"))
    }

    #[test]
    fn test_fractional_mass_sums_unweighted() {
        let a = Symbol_from(0);
        let b = Symbol_from(1);
        // Fragment 1 splits 0.5/0.5 across {A,B}; fragment 2 is unique on A;
        // an empty fragment contributes nothing.
        let frags: Vec<(Vec<Symbol>, Vec<f64>)> = vec![
            (vec![a, b], vec![0.5, 0.5]),
            (vec![a], vec![1.0]),
            (vec![], vec![]),
        ];
        let mass = fractional_mass(
            frags.iter().map(|(an, p)| (an.as_slice(), p.as_slice())),
            2,
        );
        assert!((mass[0] - 1.5).abs() < 1e-12, "{mass:?}");
        assert!((mass[1] - 0.5).abs() < 1e-12, "{mass:?}");
        // Total mass equals the number of non-empty fragments.
        assert!((mass.iter().sum::<f64>() - 2.0).abs() < 1e-12, "{mass:?}");
    }
}
