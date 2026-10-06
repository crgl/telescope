//! Reassignment models for ambiguous fragments.
//!
//! Everything downstream (reassignment modes, the report) only needs a
//! [`ModelFit`], so a new model is one more [`ReassignmentModel`] impl.
//! [`TelescopeEm`] is the reference: Telescope's `TelescopeLikelihood`, with
//! each sum taken in the order scipy takes it so results match bit for bit.

use super::numpy::{np_sum, recip0, reduceat_sum};

/// Fragment x feature score matrix in CSR form (Telescope's `raw_scores`).
/// Columns are sorted within each row; scores are the rescaled
/// `alignment score + aligned length` Telescope stores as uint16.
pub struct ScoreMatrix {
    pub n_rows: usize,
    pub n_cols: usize,
    pub indptr: Vec<usize>,
    pub indices: Vec<u32>,
    pub data: Vec<u16>,
}

impl ScoreMatrix {
    pub fn row(&self, i: usize) -> std::ops::Range<usize> {
        self.indptr[i]..self.indptr[i + 1]
    }
    /// Telescope's `Y`: does fragment `i` have more than one candidate?
    pub fn is_ambiguous(&self, i: usize) -> bool {
        self.indptr[i + 1] - self.indptr[i] > 1
    }
}

/// What a fitted model hands to reassignment and reporting. `z` and `z_init`
/// are parallel to the matrix's stored entries.
pub struct ModelFit {
    /// Membership weight of each fragment for each of its candidates.
    pub z: Vec<f64>,
    /// Entries the model dropped from `z` altogether (parallel to `z`; empty
    /// when none were). Telescope's sparse arithmetic discards a candidate
    /// whose weight underflows to exactly zero, which is not the same as
    /// keeping it with weight zero: the row is summed and ranked without it.
    pub absent: Vec<bool>,
    /// Weights before any fitting (the `init_*` report columns).
    pub z_init: Vec<f64>,
    /// Final and first-iteration feature proportions.
    pub pi: Vec<f64>,
    pub pi_init: Vec<f64>,
    pub iterations: u32,
    pub converged: bool,
    pub log_likelihood: f64,
}

pub trait ReassignmentModel {
    /// Fit the model. `progress` receives one line per iteration.
    fn fit(&self, m: &ScoreMatrix, progress: &mut dyn FnMut(&str)) -> ModelFit;
}

/// Telescope's EM (`TelescopeLikelihood`).
pub struct TelescopeEm {
    pub pi_prior: i64,
    pub theta_prior: i64,
    pub epsilon: f64,
    pub max_iter: u32,
    pub use_likelihood: bool,
}

const SCALE_FACTOR: f64 = 100.0;

struct Precomputed {
    /// `Q`: expm1 of the scaled scores
    q: Vec<f64>,
    /// row maxima of `Q`
    weights: Vec<f64>,
    total_wt: f64,
    ambig_wt: f64,
    pi_prior_wt: f64,
    theta_prior_wt: f64,
    /// column sums of `Q` over unambiguous fragments
    pisum0: Vec<f64>,
}

impl TelescopeEm {
    fn precompute(&self, m: &ScoreMatrix) -> Precomputed {
        // Q = expm1(raw * (1/max) * 100). At most 65536 distinct scores, so
        // tabulate instead of calling expm1 per entry.
        let max_score = m.data.iter().copied().max().unwrap_or(0);
        let recip = 1.0 / max_score as f64;
        let table: Vec<f64> = (0..=max_score as u32)
            .map(|s| (s as f64 * recip * SCALE_FACTOR).exp_m1())
            .collect();
        let q: Vec<f64> = m.data.iter().map(|&s| table[s as usize]).collect();

        let mut weights = vec![0.0; m.n_rows];
        let mut ambig = vec![0.0; m.n_rows];
        let mut pisum0 = vec![0.0; m.n_cols];
        for i in 0..m.n_rows {
            let r = m.row(i);
            let mut w = q[r.clone()].iter().copied().fold(f64::NEG_INFINITY, f64::max);
            if r.len() < m.n_cols {
                w = w.max(0.0); // implicit zeros take part in a sparse max
            }
            weights[i] = w;
            if m.is_ambiguous(i) {
                ambig[i] = w;
            } else {
                for k in r {
                    pisum0[m.indices[k] as usize] += q[k];
                }
            }
        }
        let wmax = weights.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        Precomputed {
            q,
            total_wt: np_sum(&weights),
            ambig_wt: np_sum(&ambig),
            weights,
            pi_prior_wt: self.pi_prior as f64 * wmax,
            theta_prior_wt: self.theta_prior as f64 * wmax,
            pisum0,
        }
    }

    /// Unnormalised E-step numerators: `Q * pi * theta` for ambiguous
    /// fragments, `Q * pi` otherwise.
    fn numerators(m: &ScoreMatrix, q: &[f64], pi: &[f64], theta: &[f64], out: &mut [f64]) {
        for i in 0..m.n_rows {
            if m.is_ambiguous(i) {
                for k in m.row(i) {
                    let c = m.indices[k] as usize;
                    out[k] = q[k] * (pi[c] * theta[c]);
                }
            } else {
                for k in m.row(i) {
                    out[k] = q[k] * pi[m.indices[k] as usize];
                }
            }
        }
    }

    fn lnl(m: &ScoreMatrix, q: &[f64], z: &[f64], pi: &[f64], theta: &[f64], buf: &mut [f64]) -> f64 {
        Self::numerators(m, q, pi, theta, buf);
        let rows: Vec<f64> = (0..m.n_rows)
            .map(|i| {
                let mut s = 0.0;
                for k in m.row(i) {
                    s += z[k] * buf[k].ln_1p();
                }
                s
            })
            .collect();
        np_sum(&rows)
    }
}

/// `csr_matrix_plus.norm(1)`: scale each row by the reciprocal of its sum.
pub fn normalize_rows(m: &ScoreMatrix, v: &mut [f64]) {
    for i in 0..m.n_rows {
        let r = m.row(i);
        let recip = recip0(reduceat_sum(&v[r.clone()]));
        for x in &mut v[r] {
            *x *= recip;
        }
    }
}

impl ReassignmentModel for TelescopeEm {
    fn fit(&self, m: &ScoreMatrix, progress: &mut dyn FnMut(&str)) -> ModelFit {
        let pre = self.precompute(m);
        let k = m.n_cols;
        let kf = k as f64;
        let mut pi = vec![1.0 / kf; k];
        let mut theta = vec![1.0 / kf; k];
        let mut pi_init = Vec::new();
        let mut z = vec![0.0; m.data.len()];
        let mut z_prev = vec![0.0; m.data.len()];
        let mut scratch = vec![0.0; if self.use_likelihood { m.data.len() } else { 0 }];
        let theta_denom = pre.ambig_wt + pre.theta_prior_wt * kf;
        let pi_denom = pre.total_wt + pre.pi_prior_wt * kf;

        let mut lnl = f64::INFINITY;
        let mut inum = 0u32;
        let mut converged = false;
        let mut reached_max = false;
        let mut thetasum = vec![0.0; k];
        let mut diffs = vec![0.0; k];
        let mut absent: Vec<bool> = Vec::new();
        let mut kept: Vec<f64> = Vec::new();
        while !(converged || reached_max) {
            // E-step into z_prev's storage, then swap.
            let z_new = &mut z_prev;
            Self::numerators(m, &pre.q, &pi, &theta, z_new);
            // Row-normalise. scipy's sparse add drops numerators that are
            // exactly zero before the row is summed, which changes how the
            // remaining values pair up in the sum.
            absent.iter_mut().for_each(|a| *a = false);
            for i in 0..m.n_rows {
                let r = m.row(i);
                let row = &mut z_new[r.clone()];
                let recip = if row.contains(&0.0) {
                    if absent.is_empty() {
                        absent = vec![false; m.data.len()];
                    }
                    kept.clear();
                    for (off, &v) in row.iter().enumerate() {
                        if v == 0.0 {
                            absent[r.start + off] = true;
                        } else {
                            kept.push(v);
                        }
                    }
                    recip0(reduceat_sum(&kept))
                } else {
                    recip0(reduceat_sum(row))
                };
                for x in row {
                    *x *= recip;
                }
            }

            // M-step: column sums run down the rows in order.
            thetasum.iter_mut().for_each(|x| *x = 0.0);
            for i in 0..m.n_rows {
                if m.is_ambiguous(i) {
                    let w = pre.weights[i];
                    for kk in m.row(i) {
                        thetasum[m.indices[kk] as usize] += z_new[kk] * w;
                    }
                }
            }
            let theta_new: Vec<f64> =
                thetasum.iter().map(|&t| (t + pre.theta_prior_wt) / theta_denom).collect();
            let pi_new: Vec<f64> = (0..k)
                .map(|j| ((pre.pisum0[j] + thetasum[j]) + pre.pi_prior_wt) / pi_denom)
                .collect();

            inum += 1;
            if inum == 1 {
                pi_init = pi_new.clone();
            }
            for j in 0..k {
                diffs[j] = (pi_new[j] - pi[j]).abs();
            }
            let diff_est = np_sum(&diffs);

            if self.use_likelihood {
                let cur = Self::lnl(m, &pre.q, z_new, &pi_new, &theta_new, &mut scratch);
                progress(&format!("Iteration {inum}, lnl= {cur:.5e}, diff={diff_est:.5e}"));
                converged = (cur - lnl).abs() < self.epsilon;
                lnl = cur;
            } else {
                progress(&format!("Iteration {inum}, diff={diff_est:.5e}"));
                converged = diff_est < self.epsilon;
            }
            reached_max = inum >= self.max_iter;
            std::mem::swap(&mut z, &mut z_prev);
            pi = pi_new;
            theta = theta_new;
        }
        drop(z_prev);
        if !self.use_likelihood {
            let mut buf = vec![0.0; m.data.len()];
            lnl = Self::lnl(m, &pre.q, &z, &pi, &theta, &mut buf);
        }

        let mut z_init = pre.q;
        normalize_rows(m, &mut z_init);
        if !absent.contains(&true) {
            absent = Vec::new();
        }
        ModelFit { z, absent, z_init, pi, pi_init, iterations: inum, converged, log_likelihood: lnl }
    }
}
