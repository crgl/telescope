use std::collections::{HashMap, HashSet};
use std::fmt;
use std::fs;
use std::io::{self, Write};
use std::path::Path;

use rayon::prelude::*;
use rayon::slice::ParallelSliceMut;

use noodles::{
    bam,
    sam::{
        Header,
        alignment::{
            io::Write as AlignmentWrite, record::data::field::Tag, record_buf::data::field::Value,
        },
    },
};

use crate::grouping::ReadGroup;
use crate::gtf::NoFeatureClass;
use crate::intern::{Interner, Symbol};
use crate::logging;
use crate::matrix::ReadAnnotationMatrix;
use crate::selection::ResolutionResult;

/// Per-annotation summary counts. Computed once during the pipeline and held
/// in `SummaryCounts::per_annotation`. Adding a new column to the summary TSV
/// is just a new field here + one increment in `compute_summary_counts`.
#[derive(Default, Clone)]
pub struct AnnotationCounts {
    /// Number of read groups in which this annotation appeared as a representative.
    pub initial_count: usize,
    /// Number of read groups for which this annotation was the confident assignment.
    pub final_count: usize,
    /// Number of read groups for which this annotation was the *sole* row entry —
    /// equivalently, the assignment that would survive a confidence threshold of 1.0.
    pub unique_count: usize,
}

/// Single source of truth for all summary statistics produced by the pipeline.
/// Built once via `compute_summary_counts` and threaded through to `write_summary`,
/// `RunSummary`, and the verbose log lines — replacing the four separate ad-hoc
/// counting passes that used to live in `write_summary` and `main.rs`.
#[derive(Default)]
pub struct SummaryCounts {
    pub per_annotation: HashMap<Symbol, AnnotationCounts>,
    pub total_groups: usize,
    pub num_confident: usize,
    pub num_ambiguous: usize,
    pub num_no_feature: usize,
}

impl SummaryCounts {
    /// Combine two partial `SummaryCounts` produced by parallel fold chunks.
    fn merge(mut self, other: Self) -> Self {
        self.total_groups += other.total_groups;
        self.num_confident += other.num_confident;
        self.num_ambiguous += other.num_ambiguous;
        self.num_no_feature += other.num_no_feature;
        for (sym, counts) in other.per_annotation {
            let entry = self.per_annotation.entry(sym).or_default();
            entry.initial_count += counts.initial_count;
            entry.final_count += counts.final_count;
            entry.unique_count += counts.unique_count;
        }
        self
    }
}

/// Compute all summary counts in a single parallel fold-reduce over the read
/// groups. Reads `matrix.row(i)` and `resolutions[i]` together so each group is
/// touched exactly once. Replaces the sequential loops in the old `write_summary`
/// plus the scattered `num_confident`/`num_no_feature`/`num_ambiguous`/
/// `unique_annotations` counting that used to live in `main.rs`.
pub fn compute_summary_counts(
    matrix: &ReadAnnotationMatrix,
    resolutions: &[ResolutionResult],
    no_feature_class: NoFeatureClass,
) -> SummaryCounts {
    (0..matrix.num_groups())
        .into_par_iter()
        .fold(SummaryCounts::default, |mut acc, i| {
            let row = matrix.row(i);
            let res = &resolutions[i];
            acc.total_groups += 1;

            for (ann, _) in row {
                acc.per_annotation.entry(*ann).or_default().initial_count += 1;
            }
            for &ann in res.confident.keys() {
                acc.per_annotation.entry(ann).or_default().final_count += 1;
            }
            // unique == "would be confident at threshold 1.0" == row of size 1.
            // Per the softmax in selection.rs, prob == 1.0 only when n == 1.
            if row.len() == 1 {
                let only = row[0].0;
                acc.per_annotation.entry(only).or_default().unique_count += 1;
            }

            // Group categorization. A confident pick of __no_feature__ (or any
            // cytoband) counts toward num_no_feature (not num_confident), so the
            // run-summary aggregates across bands; no confident pick is
            // ambiguous; everything else is a confident real-annotation pick.
            if res.confident.is_empty() {
                acc.num_ambiguous += 1;
            } else if res
                .confident
                .keys()
                .any(|&k| no_feature_class.is_no_feature(k))
            {
                acc.num_no_feature += 1;
            } else {
                acc.num_confident += 1;
            }
            acc
        })
        .reduce(SummaryCounts::default, SummaryCounts::merge)
}

/// Tags to apply to a BAM record during the output pass.
///
/// Annotation strings are kept as interned `Symbol`s here and only resolved to
/// bytes inside `write_bam_streaming`. With ~28M directives in flight this
/// saves ~80 B/directive vs. the previous `String`-carrying layout — roughly
/// 2-3 GB on a real dataset, with no behavioral change in the output BAM.
#[derive(Clone, Copy)]
pub struct OutputDirective {
    pub zb: Symbol,
    pub zf: Option<Symbol>,
    /// `None` → no ZR tag. `Some(true)` → "1" (representative).
    /// `Some(false)` → "0". Only set when `--all-alignments` is in effect.
    pub zr: Option<bool>,
}

/// Build a sorted list of `(input_order, OutputDirective)` for all records that
/// should be output. The list is sorted ascending by `input_order` so that
/// `write_bam_streaming` can walk it with a cursor in lockstep with the input
/// BAM stream — no hashing needed.
///
/// Implementation note: a previous version returned `HashMap<usize, OutputDirective>`
/// built by parallel-mapping each group to a small HashMap and then sequentially
/// `extend`-ing 34M of them into one 68M-entry HashMap. That sequential merge
/// dominated the post-grouping phase (~280s on full_JEG3_tester.bam). Switching
/// to a flat `Vec<(usize, OutputDirective)>` flattened in parallel via
/// `flat_map_iter` and then `par_sort_unstable_by_key` eliminates both the
/// per-group HashMap allocations and the merge step.
#[allow(clippy::too_many_arguments)]
pub fn build_output_directives(
    groups: &[ReadGroup],
    matrix: &ReadAnnotationMatrix,
    resolutions: &[ResolutionResult],
    all_alignments: bool,
    include_no_feature: bool,
    no_feature_class: NoFeatureClass,
) -> Vec<(usize, OutputDirective)> {
    // The generic __no_feature__ symbol is the ZB/ZF fallback; the class also
    // recognises every cytoband so band reads are gated out of the BAM unless
    // --include-no-feature, exactly like __no_feature__.
    let no_feature_symbol = no_feature_class.no_feature_symbol();
    let mut directives: Vec<(usize, OutputDirective)> = (0..groups.len())
        .into_par_iter()
        .flat_map_iter(|group_id| -> Vec<(usize, OutputDirective)> {
            let group = &groups[group_id];
            let row = matrix.row(group_id);
            let resolution = &resolutions[group_id];
            let mut local: Vec<(usize, OutputDirective)> = Vec::new();

            let rep_indices: HashSet<usize> = row.iter().map(|(_, c)| c.alignment_idx).collect();
            let confident_indices: HashSet<usize> =
                resolution.confident.values().copied().collect();

            // Build map: alignment index → annotation symbols it represents (initial)
            let mut idx_to_annotations: HashMap<usize, Vec<Symbol>> = HashMap::new();
            for (annotation, cell) in row {
                idx_to_annotations
                    .entry(cell.alignment_idx)
                    .or_default()
                    .push(*annotation);
            }

            // Build map: alignment index → confident annotation symbol
            let mut idx_to_confident: HashMap<usize, Vec<Symbol>> = HashMap::new();
            for (&annotation, &idx) in &resolution.confident {
                idx_to_confident.entry(idx).or_default().push(annotation);
            }

            // No confident annotation → __no_feature__ (default mode only)
            if resolution.confident.is_empty() && !all_alignments {
                if include_no_feature && let Some(entry) = group.alignments.first() {
                    let directive = OutputDirective {
                        zb: no_feature_symbol,
                        zf: Some(no_feature_symbol),
                        zr: None,
                    };
                    if let Some(mate_io) = entry.mate_input_order {
                        local.push((mate_io, directive));
                    }
                    local.push((entry.input_order, directive));
                }
                return local;
            }

            for (idx, entry) in group.alignments.iter().enumerate() {
                let is_rep = rep_indices.contains(&idx);
                let is_confident = confident_indices.contains(&idx);

                if !all_alignments && !is_confident {
                    continue;
                }

                // Build overlap lookup for this entry (Symbol → overlap_bp)
                let overlap_map: HashMap<Symbol, usize> = entry
                    .annotations
                    .iter()
                    .map(|ao| (ao.annotation, ao.overlap_bp))
                    .collect();

                // Determine initial annotation (ZB)
                let zb = if is_rep {
                    let rep_annotations = idx_to_annotations.get(&idx);
                    match rep_annotations {
                        Some(anns) if !anns.is_empty() => *anns
                            .iter()
                            .max_by_key(|&&ann| overlap_map.get(&ann).copied().unwrap_or(0))
                            .unwrap(),
                        _ => no_feature_symbol,
                    }
                } else {
                    entry
                        .annotations
                        .iter()
                        .max_by_key(|ao| ao.overlap_bp)
                        .map(|ao| ao.annotation)
                        .unwrap_or(no_feature_symbol)
                };

                // Determine confident annotation (ZF)
                let zf = if is_confident {
                    Some(
                        idx_to_confident
                            .get(&idx)
                            .and_then(|anns| {
                                anns.iter()
                                    .max_by_key(|&&ann| overlap_map.get(&ann).copied().unwrap_or(0))
                            })
                            .copied()
                            .unwrap_or(no_feature_symbol),
                    )
                } else {
                    None
                };

                // Default mode: suppress records whose initial OR confident
                // annotation is __no_feature__ or a cytoband. With
                // --include-no-feature this gate is disabled.
                let zf_is_no_feature = zf.is_some_and(|z| no_feature_class.is_no_feature(z));
                if (no_feature_class.is_no_feature(zb) || zf_is_no_feature) && !include_no_feature {
                    continue;
                }

                // ZR tag (only when --all-alignments)
                let zr = if all_alignments { Some(is_rep) } else { None };

                let directive = OutputDirective { zb, zf, zr };
                if let Some(mate_io) = entry.mate_input_order {
                    local.push((mate_io, directive));
                }
                local.push((entry.input_order, directive));
            }

            local
        })
        .collect();

    // input_order values are unique across groups, so a stable sort is unnecessary.
    directives.par_sort_unstable_by_key(|(io, _)| *io);
    directives
}

/// Write the output BAM by streaming through the input BAM and applying directives.
/// Mutates record_buf in place to avoid cloning entire BAM records.
///
/// `directives` must be sorted ascending by `input_order` (as produced by
/// `build_output_directives`). The function walks both the BAM stream and the
/// directive slice in lockstep with a cursor — O(n) with no hashing.
pub fn write_bam_streaming(
    input_bam: &str,
    output_dir: &Path,
    input_stem: &str,
    header: &Header,
    directives: &[(usize, OutputDirective)],
    interner: &Interner,
) -> io::Result<()> {
    let bam_path = output_dir.join(format!("{}_annotated.bam", input_stem));
    let file = fs::File::create(&bam_path)?;
    let mut writer = bam::io::Writer::new(file);
    writer.write_header(header)?;

    let zb_tag = Tag::new(b'Z', b'B');
    let zf_tag = Tag::new(b'Z', b'F');
    let zr_tag = Tag::new(b'Z', b'R');

    let mut reader = bam::io::reader::Builder.build_from_path(input_bam)?;
    let _ = reader.read_header()?;

    let mut record_buf = noodles::sam::alignment::RecordBuf::default();
    let mut input_order: usize = 0;
    let mut cursor: usize = 0;

    loop {
        match reader.read_record_buf(header, &mut record_buf) {
            Ok(0) => break,
            Ok(_) => {
                if cursor < directives.len() && directives[cursor].0 == input_order {
                    let directive = &directives[cursor].1;
                    // Resolve interned Symbols → &str at write time. Avoids
                    // ~28M heap String allocations during build_output_directives.
                    let zb_str = interner.resolve(directive.zb);
                    record_buf
                        .data_mut()
                        .insert(zb_tag, Value::String(zb_str.into()));
                    if let Some(zf_sym) = directive.zf {
                        let zf_str = interner.resolve(zf_sym);
                        record_buf
                            .data_mut()
                            .insert(zf_tag, Value::String(zf_str.into()));
                    }
                    if let Some(is_rep) = directive.zr {
                        let zr_str = if is_rep { "1" } else { "0" };
                        record_buf
                            .data_mut()
                            .insert(zr_tag, Value::String(zr_str.into()));
                    }
                    writer.write_alignment_record(header, &record_buf)?;
                    cursor += 1;
                }
                input_order += 1;
            }
            Err(e) => return Err(e),
        }
    }

    Ok(())
}

/// Write all outputs to the output directory. Expects `directives` to be
/// precomputed (in `main`, so the cost of `build_output_directives` can be
/// measured independently). When `debug` is true, `write_bam_streaming` is
/// skipped; jaccard and summary TSVs are still produced.
///
/// `skip_jaccard` suppresses the `{stem}_jaccard.tsv` file entirely (no header,
/// no rows). When `false` the file is always created — even if no pairs were
/// found, a header-only file is written.
#[allow(clippy::too_many_arguments)]
pub fn write_outputs(
    input_bam: &str,
    output_dir: &Path,
    input_stem: &str,
    header: &Header,
    directives: &[(usize, OutputDirective)],
    summary_counts: &SummaryCounts,
    jaccard_tuples: &[(String, String, f64, f64)],
    interner: &Interner,
    lengths: &[u64],
    prior_mass: &[f64],
    posterior_mass: &[f64],
    debug: bool,
    skip_jaccard: bool,
) -> io::Result<()> {
    fs::create_dir_all(output_dir)?;

    if !debug {
        write_bam_streaming(
            input_bam, output_dir, input_stem, header, directives, interner,
        )?;
    }
    if !skip_jaccard {
        write_jaccard(output_dir, input_stem, jaccard_tuples)?;
    }
    write_summary(
        output_dir,
        input_stem,
        summary_counts,
        interner,
        lengths,
        prior_mass,
        posterior_mass,
    )?;

    Ok(())
}

fn write_jaccard(
    output_dir: &Path,
    input_stem: &str,
    tuples: &[(String, String, f64, f64)],
) -> io::Result<()> {
    let tsv_path = output_dir.join(format!("{}_jaccard.tsv", input_stem));
    let mut file = fs::File::create(&tsv_path)?;

    writeln!(file, "annotation_a\tannotation_b\tjaccard\toverlap")?;
    for (a, b, j, o) in tuples {
        writeln!(file, "{}\t{}\t{:.2}\t{:.2}", a, b, j, o)?;
    }

    Ok(())
}

/// Stream Jaccard pairs straight from a `SparseMatrix` to the `_jaccard.tsv`
/// file, with no in-memory tuples buffer and no sort. Used by `--low-memory`.
/// Output is unsorted; the user can `sort` it post-hoc if needed.
pub fn write_jaccard_streaming(
    output_dir: &Path,
    input_stem: &str,
    sparse_matrix: &crate::matrix::SparseMatrix,
    interner: &Interner,
    similarity_threshold: f64,
    min_reads: usize,
) -> io::Result<()> {
    let tsv_path = output_dir.join(format!("{}_jaccard.tsv", input_stem));
    let file = fs::File::create(&tsv_path)?;
    let mut writer = std::io::BufWriter::new(file);
    sparse_matrix.write_jaccard_streaming(
        interner,
        similarity_threshold,
        min_reads,
        &mut writer,
    )?;
    writer.flush()?;
    Ok(())
}

fn write_summary(
    output_dir: &Path,
    input_stem: &str,
    counts: &SummaryCounts,
    interner: &Interner,
    lengths: &[u64],
    prior_mass: &[f64],
    posterior_mass: &[f64],
) -> io::Result<()> {
    let tsv_path = output_dir.join(format!("{}_summary.tsv", input_stem));
    let mut file = fs::File::create(&tsv_path)?;

    // Resolve symbols to strings and sort by (final_count DESC, initial_count DESC).
    let mut rows: Vec<(String, Symbol, &AnnotationCounts)> = counts
        .per_annotation
        .iter()
        .map(|(sym, c)| (interner.resolve(*sym).to_string(), *sym, c))
        .collect();
    rows.sort_by(|a, b| {
        b.2.final_count
            .cmp(&a.2.final_count)
            .then_with(|| b.2.initial_count.cmp(&a.2.initial_count))
    });

    writeln!(
        file,
        "annotation\tfinal_count\tinitial_count\tunique_count\tlength\tprior_mass\tposterior_mass"
    )?;
    for (ann, sym, c) in &rows {
        let length = lengths.get(sym.idx()).copied().unwrap_or(0);
        let prior = prior_mass.get(sym.idx()).copied().unwrap_or(0.0);
        let posterior = posterior_mass.get(sym.idx()).copied().unwrap_or(0.0);
        writeln!(
            file,
            "{}\t{}\t{}\t{}\t{}\t{:.6}\t{:.6}",
            ann, c.final_count, c.initial_count, c.unique_count, length, prior, posterior
        )?;
    }

    Ok(())
}

/// Aggregate statistics for a pipeline run.
pub struct RunSummary {
    pub total_records: usize,
    pub dropped_unmapped: usize,
    pub num_groups: usize,
    pub num_dropped_no_feature: usize,
    pub num_confident: usize,
    pub num_ambiguous: usize,
    pub num_no_feature: usize,
    pub validation_failures: usize,
    pub confidence: f64,
    pub em_iterations: usize,
    pub em_stop_reason: EmStopReason,
    pub elapsed_secs: f64,
    pub peak_rss_bytes: Option<u64>,
    pub threads_used: usize,
}

/// Why EM stopped, reported alongside the iteration count in the run summary.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EmStopReason {
    /// Max per-annotation mass delta fell below `--min-mass-delta`.
    Converged,
    /// Hit the `--max-iter` cap without converging.
    MaxIter,
    /// EM did not run (`--max-iter 0`, or no ambiguous reads to iterate).
    Skipped,
}

impl fmt::Display for EmStopReason {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s = match self {
            EmStopReason::Converged => "converged",
            EmStopReason::MaxIter => "max-iter reached",
            EmStopReason::Skipped => "skipped",
        };
        f.write_str(s)
    }
}

impl fmt::Display for RunSummary {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "rusty_telescope run summary")?;
        writeln!(f, "==========================")?;
        writeln!(f, "Total BAM records:    {:>10}", self.total_records)?;
        if self.dropped_unmapped > 0 {
            writeln!(f, "Dropped (unmapped):   {:>10}", self.dropped_unmapped)?;
        }
        writeln!(f, "Read groups:          {:>10}", self.num_groups)?;
        if self.num_dropped_no_feature > 0 {
            writeln!(
                f,
                "Dropped (no overlap): {:>10}",
                self.num_dropped_no_feature,
            )?;
        }
        writeln!(f, "Confident assignments:{:>10}", self.num_confident)?;
        writeln!(f, "Ambiguous (dropped):  {:>10}", self.num_ambiguous)?;
        writeln!(f, "No feature:           {:>10}", self.num_no_feature)?;
        if self.validation_failures > 0 {
            writeln!(f, "Pair validation fails:{:>10}", self.validation_failures)?;
        }
        writeln!(f, "Confidence threshold: {:>10.2}", self.confidence)?;
        writeln!(
            f,
            "EM iterations:        {:>10}  ({})",
            self.em_iterations, self.em_stop_reason,
        )?;
        writeln!(f, "Threads used:         {:>10}", self.threads_used)?;
        writeln!(f, "Elapsed time:         {:>8.2}s", self.elapsed_secs)?;
        if let Some(b) = self.peak_rss_bytes {
            writeln!(f, "Peak RAM:             {:>10}", logging::format_bytes(b))?;
        }
        Ok(())
    }
}

impl RunSummary {
    pub fn print_stderr(&self) {
        eprint!("{}", self);
    }
}

/// Write run summary to a text file in the output directory.
pub fn write_run_summary(
    output_dir: &Path,
    input_stem: &str,
    summary: &RunSummary,
) -> io::Result<()> {
    let path = output_dir.join(format!("{}_run_summary.txt", input_stem));
    let mut file = fs::File::create(&path)?;
    write!(file, "{}", summary)?;
    Ok(())
}
