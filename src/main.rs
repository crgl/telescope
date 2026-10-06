mod cli;
mod cytoband;
mod detect_strand;
mod grouping;
mod gtf;
mod intern;
mod legacy;
mod logging;
mod matrix;
mod output;
mod overlap;
mod selection;

use std::io;
use std::path::Path;

use clap::Parser;
use noodles::{
    bam,
    sam::alignment::{RecordBuf, record::data::field::Tag},
};
use rayon::prelude::*;

use cli::{
    AnnotateArgs, BandNoFeature, Cli, LengthCorrection, ModelType, StrandedMode, Subcommand,
};
use grouping::AlignmentEntry;
use gtf::NoFeatureClass;
use logging::{Logger, Verbosity};
use overlap::{AnnotationOverlap, StrandFilter};

/// Lightweight BAM record fields extracted during sequential read,
/// before parallel annotation. Borrows chrom name from the BAM header.
pub struct RawRecord<'a> {
    pub qname: Option<String>,
    pub input_order: usize,
    pub alignment_score: Option<i64>,
    pub is_proper_pair: bool,
    pub span: Option<(&'a str, usize, usize)>,
    // Mate pairing fields
    pub is_paired: bool,
    pub is_first_segment: bool,
    pub is_last_segment: bool,
    pub is_reverse: bool,
    pub is_mate_reverse: bool,
    pub is_mate_unmapped: bool,
    pub mate_ref_id: Option<usize>,
    pub mate_alignment_start: Option<usize>,
    pub ref_id: Option<usize>,
    pub alignment_start_0based: Option<usize>,
}

/// Bridge between `--stranded` library mode + the alignment-unit's strand
/// (read1 strand, or the reverse of read2 for read2-only fragments) into a
/// `StrandFilter` consumed by `annotate_from_span(s)`.
///
/// Returns `None` when no filtering should occur — either `--stranded` was not
/// passed, or the alignment unit has no determinable strand (e.g., unmapped).
pub fn build_strand_filter(
    mode: Option<StrandedMode>,
    unit_strand: Option<char>,
) -> Option<StrandFilter> {
    let mode = mode?;
    let s = unit_strand?;
    Some(match mode {
        StrandedMode::Fr | StrandedMode::F => StrandFilter::Same(s),
        StrandedMode::Rf | StrandedMode::R => StrandFilter::Opposite(s),
    })
}

/// Resolve `--length-correction` (auto/on/off) into a concrete bool given the
/// selected prior-probability `model`. `auto` enables length correction for
/// `ribbonfish` and disables it for `telescope`.
fn resolve_length_correction(setting: LengthCorrection, model: ModelType) -> bool {
    match setting {
        LengthCorrection::On => true,
        LengthCorrection::Off => false,
        LengthCorrection::Auto => matches!(model, ModelType::Ribbonfish),
    }
}

/// Resolve `--band-no-feature` (auto/on/off) into a concrete bool given the
/// selected prior-probability `model`. `auto` enables banding for `ribbonfish`
/// and disables it for `telescope`.
fn resolve_band_no_feature(setting: BandNoFeature, model: ModelType) -> bool {
    match setting {
        BandNoFeature::On => true,
        BandNoFeature::Off => false,
        BandNoFeature::Auto => matches!(model, ModelType::Ribbonfish),
    }
}

/// Whether an alignment unit overlaps a *real* GTF feature — i.e. a non-band,
/// non-`__no_feature__` annotation. Band assignment writes a band symbol into
/// `annotations`, so this is what decides "processed vs dropped": banding must
/// only relabel surviving no-feature reads, never make a bandless group eligible.
fn overlaps_real_feature(entry: &AlignmentEntry, class: &NoFeatureClass) -> bool {
    entry
        .annotations
        .iter()
        .any(|ao| !class.is_no_feature(ao.annotation))
}

fn main() -> io::Result<()> {
    let cli = Cli::parse();
    match cli.command {
        Subcommand::Annotate(args) => run_annotate(args),
        Subcommand::DetectStrand(args) => detect_strand::run_detect_strand(args),
        Subcommand::Assign(args) => legacy::run_assign(args),
    }
}

fn run_annotate(args: AnnotateArgs) -> io::Result<()> {
    if args.low_memory && args.all_alignments {
        eprintln!("error: --low-memory is incompatible with --all-alignments");
        std::process::exit(2);
    }
    // --low-memory streams Jaccard pairs straight to disk (no in-memory buffer,
    // no sort). --skip-jaccard short-circuits both paths.
    let stream_jaccard = args.low_memory && !args.skip_jaccard;
    // --low-memory defaults to single-threaded; explicit --threads N overrides.
    let effective_threads = if args.low_memory && args.threads == 0 {
        1
    } else {
        args.threads
    };

    let verbosity = if args.quiet {
        Verbosity::Quiet
    } else if args.verbose {
        Verbosity::Verbose
    } else {
        Verbosity::Normal
    };
    let log = Logger::new(verbosity);

    // Step 1: Configure thread pool
    if effective_threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(effective_threads)
            .build_global()
            .expect("failed to initialize thread pool");
    }

    // Step 2: Load GTF (streaming, line-by-line)
    let mut index = gtf::FeatureIndex::from_file(&args.gtf, &args.field)?;
    log.stage(&format!(
        "Loaded GTF: {} features across {} chromosomes",
        index.num_features(),
        index.num_chromosomes(),
    ));

    // Step 2b: Optionally load cytobands to partition __no_feature__ reads into
    // genomic bands. Resolves on for ribbonfish / off for telescope under
    // `auto`. A missing file surfaces as an io::Error, consistent with --gtf.
    let band_no_feature = resolve_band_no_feature(args.band_no_feature, args.model);
    if band_no_feature {
        index.attach_bands(&args.cytoband)?;
        log.stage(&format!("Loaded {} cytobands", index.num_bands()));
    }

    // === Pass 1: Analyze ===
    // Phase A: Stream BAM sequentially, extract lightweight fields.
    // Phase B: Annotate all records in parallel using rayon.

    let mut reader = bam::io::reader::Builder.build_from_path(&args.bam)?;
    let header = reader.read_header()?;
    let as_tag = Tag::new(b'A', b'S');

    let mut raw_records: Vec<RawRecord<'_>> = Vec::new();
    let mut record_buf = RecordBuf::default();
    let mut input_order: usize = 0;
    let mut dropped_unmapped: usize = 0;
    let mut global_as_min: Option<i64> = None;
    let mut global_as_max: Option<i64> = None;

    loop {
        match reader.read_record_buf(&header, &mut record_buf) {
            Ok(0) => break,
            Ok(_) => {
                // Drop unmapped reads (SAM 0x4 flag) as they're read: count them,
                // otherwise ignore (no pairing, annotation, or output). This
                // covers both RNAME='*' reads and "placed unmapped" reads (CIGAR
                // '*' carrying a mate's RNAME/POS). `input_order` still advances so
                // it stays the true BAM record index, keeping the output-write
                // pass — which re-reads every record — in lockstep.
                if record_buf.flags().is_unmapped() {
                    dropped_unmapped += 1;
                    input_order += 1;
                    continue;
                }
                let qname = record_buf.name().map(|n| n.to_string());
                let span = overlap::reference_span(&record_buf, &header);
                let alignment_score = record_buf.data().get(&as_tag).and_then(|v| v.as_int());
                if let Some(s) = alignment_score {
                    global_as_min = Some(global_as_min.map_or(s, |m| m.min(s)));
                    global_as_max = Some(global_as_max.map_or(s, |m| m.max(s)));
                }
                let flags = record_buf.flags();

                let mate_alignment_start = record_buf
                    .mate_alignment_start()
                    .map(|pos| usize::from(pos) - 1);
                let alignment_start_0based =
                    record_buf.alignment_start().map(|pos| usize::from(pos) - 1);

                raw_records.push(RawRecord {
                    qname,
                    input_order,
                    alignment_score,
                    is_proper_pair: flags.is_properly_segmented(),
                    span,
                    is_paired: flags.is_segmented(),
                    is_first_segment: flags.is_first_segment(),
                    is_last_segment: flags.is_last_segment(),
                    is_reverse: flags.is_reverse_complemented(),
                    is_mate_reverse: flags.is_mate_reverse_complemented(),
                    is_mate_unmapped: flags.is_mate_unmapped(),
                    mate_ref_id: record_buf.mate_reference_sequence_id(),
                    mate_alignment_start,
                    ref_id: record_buf.reference_sequence_id(),
                    alignment_start_0based,
                });
                input_order += 1;
            }
            Err(e) => return Err(e),
        }
    }

    log.stage(&format!("Read {} BAM records", input_order));
    if dropped_unmapped > 0 {
        log.stage(&format!("Dropped {} unmapped alignments", dropped_unmapped));
    }

    // Drop the reader to release BAM I/O resources
    drop(reader);
    drop(record_buf);

    // Phase A5: Pair mates (PE mode) or pass through (SE mode)
    let pairing = grouping::pair_mates(raw_records, args.single_end);
    let mut alignment_units = pairing.units;
    let pairing_stats = pairing.stats;

    if pairing_stats.total_failures() > 0 {
        log.stage(&format!(
            "Mate pairing: {} validation failures (treated as unpaired)",
            pairing_stats.total_failures(),
        ));
        if verbosity == Verbosity::Verbose {
            log.detail(&format!(
                "Failure breakdown: {} flag, {} ref, {} position",
                pairing_stats.flag_failures, pairing_stats.ref_failures, pairing_stats.pos_failures,
            ));
        }
    }

    // Filter discordant pairs if requested
    if args.exclude_discordant {
        let before = alignment_units.len();
        alignment_units.retain(|u| {
            if u.spans().len() < 2 {
                return true; // single-read units are never discordant
            }
            // Concordant if both spans on same chromosome
            u.spans()[0].0 == u.spans()[1].0
        });
        let discarded = before - alignment_units.len();
        if discarded > 0 {
            log.stage(&format!("Excluded {} discordant pairs", discarded));
        }
    }

    log.stage(&format!(
        "Paired into {} alignment units ({} mate pairs, {} unpaired)",
        alignment_units.len(),
        alignment_units
            .iter()
            .filter(|u| u.mate_input_order.is_some())
            .count(),
        alignment_units
            .iter()
            .filter(|u| u.mate_input_order.is_none())
            .count(),
    ));

    // Verbose diagnostic: where do alignment units that can't be placed in any
    // cytoband land? Surfaces reference-naming mismatches (e.g. unplaced contigs
    // named with GenBank accessions) vs genuinely intergenic reads. Read-only
    // over alignment_units; only runs under --verbose with band mode on.
    if band_no_feature && verbosity == Verbosity::Verbose {
        use std::collections::HashMap;
        let unbanded: HashMap<&str, u64> = alignment_units
            .par_iter()
            .filter_map(|u| match u.spans().first() {
                None => Some("(unmapped)"),
                Some(_) if index.assign_band(u.spans()).is_none() => Some(u.spans()[0].0),
                Some(_) => None,
            })
            .fold(HashMap::new, |mut m, r| {
                *m.entry(r).or_insert(0) += 1;
                m
            })
            .reduce(HashMap::new, |mut a, b| {
                for (k, v) in b {
                    *a.entry(k).or_insert(0) += v;
                }
                a
            });
        let total: u64 = unbanded.values().sum();
        let unmapped = unbanded.get("(unmapped)").copied().unwrap_or(0);
        log.detail(&format!(
            "Out-of-cytoband alignment units: {} (unmapped: {})",
            total, unmapped
        ));
        let mut by_ref: Vec<(&str, u64)> = unbanded.into_iter().collect();
        by_ref.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(b.0)));
        for (rname, count) in by_ref.into_iter().take(15) {
            log.detail(&format!("no band: {rname}\t{count}"));
        }
    }

    // Phase B: Annotate in parallel (using merged spans for mate pairs)
    let stranded = args.stranded;
    let entries: Vec<AlignmentEntry> = alignment_units
        .into_par_iter()
        .map(|unit| {
            let strand_filter = build_strand_filter(stranded, unit.strand);
            let mut annotations = overlap::annotate_from_spans(
                unit.spans(),
                &index,
                &args.field,
                args.normalize_chr,
                strand_filter,
                args.min_overlap,
            );
            // Band mode: a read overlapping no GTF feature is partitioned into
            // its best-overlapping cytoband instead of the generic
            // __no_feature__. `assign_band` returns None when band mode is off
            // or no span hits a band, leaving the generic fallback path intact.
            if annotations.is_empty()
                && let Some((band, overlap_bp)) = index.assign_band(unit.spans())
            {
                annotations.push(AnnotationOverlap {
                    annotation: band,
                    overlap_bp,
                });
            }
            AlignmentEntry {
                qname: unit.qname,
                input_order: unit.input_order,
                mate_input_order: unit.mate_input_order,
                annotations,
                alignment_score: unit.alignment_score,
                is_proper_pair: unit.is_proper_pair,
            }
        })
        .collect();

    log.stage(&format!(
        "Annotated {} alignment units in parallel",
        entries.len()
    ));
    // Classifier recognising __no_feature__ + every band; used to decide which
    // groups are processed (real-feature overlap only) so banding never changes
    // the processed/dropped set.
    let no_feature_class = index.no_feature_class();
    if verbosity == Verbosity::Verbose {
        let annotated = entries
            .iter()
            .filter(|e| overlaps_real_feature(e, &no_feature_class))
            .count();
        log.detail(&format!("{} units with ≥1 annotation", annotated));
    }

    // Step 4: Group by QNAME, then drop groups with no alignment overlapping a
    // *real* GTF feature (not included in biological model; would only inflate
    // __no_feature__'s θ mass during EM and bloat downstream outputs). Band
    // symbols do NOT count here: banding only relabels the no-feature reads in
    // groups that already overlap a real feature, so the processed/dropped set
    // is invariant to --band-no-feature and --model.
    let mut groups = grouping::group_by_qname(entries);
    let pre_drop_count = groups.len();
    groups.retain(|g| {
        g.alignments
            .iter()
            .any(|a| overlaps_real_feature(a, &no_feature_class))
    });
    let num_dropped_no_feature = pre_drop_count - groups.len();
    if num_dropped_no_feature > 0 {
        log.stage(&format!(
            "Dropped {} read groups with no annotation overlap",
            num_dropped_no_feature,
        ));
    }
    log.stage(&format!("Grouped into {} read groups", groups.len()));

    // Step 5: Build read-annotation matrix (parallel)
    let no_feature_symbol = index.no_feature_symbol;
    let read_ann_matrix = matrix::ReadAnnotationMatrix::new(
        groups
            .par_iter()
            .map(|g| selection::select_representatives(g, no_feature_symbol))
            .collect(),
    );
    log.stage("Built read-annotation matrix");

    // Step 5.5a: Compute per-fragment priors (model-specific, parallel)
    let priors: Vec<selection::FragmentPriors> = (0..groups.len())
        .into_par_iter()
        .map(|i| {
            selection::compute_priors(
                read_ann_matrix.row(i),
                args.model,
                global_as_min,
                global_as_max,
            )
        })
        .collect();
    log.stage("Computed per-fragment priors");

    // Fractional prior assignment per annotation (Σ prior over fragments,
    // unweighted). Computed here because the EM call below moves `priors`.
    let prior_mass = selection::fractional_mass(
        priors
            .iter()
            .map(|fp| (fp.annotations.as_slice(), fp.priors.as_slice())),
        index.interner.len(),
    );

    // Step 5.5b: EM reassignment (sequential outer loop, parallel per iteration).
    // EM consumes `priors` and returns owned `FragmentPosteriors` carrying the
    // annotation list, so the priors Vec is freed here rather than living
    // alongside posteriors through resolution and output.
    let length_correction = resolve_length_correction(args.length_correction, args.model);
    let selection::EmOutcome {
        posteriors,
        iterations: em_iterations,
    } = match args.model {
        ModelType::Ribbonfish => selection::expectation_maximization(
            priors,
            index.interner.len(),
            args.max_iter,
            args.min_mass_delta,
            &index.lengths,
            length_correction,
        ),
        ModelType::Telescope => {
            // Per-fragment weights = max transformed AS, computed in the same
            // parallel pattern as the priors above (index-parallel to groups).
            let weights: Vec<f64> = (0..groups.len())
                .into_par_iter()
                .map(|i| {
                    selection::fragment_weight(read_ann_matrix.row(i), global_as_min, global_as_max)
                })
                .collect();
            selection::telescope_em(
                priors,
                &weights,
                index.interner.len(),
                args.max_iter,
                args.min_mass_delta,
                args.theta_prior as f64,
                global_as_min,
                global_as_max,
                &index.lengths,
                length_correction,
            )
        }
    };
    // EM stops at the first of --max-iter or max-mass-delta convergence.
    // 0 iterations means EM never ran: --max-iter 0, or (telescope) no ambiguous
    // reads to iterate. Fewer than max_iter but >0 means it converged early.
    let em_stop_reason = if em_iterations == 0 {
        output::EmStopReason::Skipped
    } else if em_iterations < args.max_iter {
        output::EmStopReason::Converged
    } else {
        output::EmStopReason::MaxIter
    };
    match args.model {
        ModelType::Ribbonfish => log.stage(&format!(
            "Ran EM for {} iterations ({}, length-correction {})",
            em_iterations,
            em_stop_reason,
            if length_correction { "on" } else { "off" },
        )),
        ModelType::Telescope => log.stage(&format!(
            "Ran telescope EM for {} iterations ({}, theta-prior {}, length-correction {})",
            em_iterations,
            em_stop_reason,
            args.theta_prior,
            if length_correction { "on" } else { "off" },
        )),
    }

    // Fractional posterior assignment per annotation (Σ γ over fragments,
    // unweighted). For ribbonfish this equals the final θ; for telescope it is
    // Σγ (expected fragment count).
    let posterior_mass = selection::fractional_mass(
        posteriors
            .iter()
            .map(|fp| (fp.annotations.as_slice(), fp.probs.as_slice())),
        index.interner.len(),
    );

    // Step 5.5c: Confidence pick + threshold on final γ (parallel)
    let resolutions: Vec<selection::ResolutionResult> = (0..groups.len())
        .into_par_iter()
        .map(|i| {
            selection::select_confident(
                &groups[i],
                read_ann_matrix.row(i),
                &posteriors[i],
                args.confidence,
            )
        })
        .collect();
    log.stage("Resolved ambiguity");

    // Step 5.6: Compute all summary counts in one parallel fold-reduce.
    // Replaces four separate ad-hoc passes (write_summary's two loops + the
    // num_confident/num_no_feature filters + the unique_annotations HashSet).
    let summary_counts =
        output::compute_summary_counts(&read_ann_matrix, &resolutions, no_feature_class);
    log.stage(&format!(
        "Computed summary counts ({} confident, {} ambiguous, {} no-feature)",
        summary_counts.num_confident, summary_counts.num_ambiguous, summary_counts.num_no_feature,
    ));
    if verbosity == Verbosity::Verbose {
        log.detail(&format!(
            "{} unique annotations",
            summary_counts.per_annotation.len()
        ));
    }

    // Step 6+7: Build Jaccard co-occurrence matrix and either collect+sort
    // pairs (default mode) or stream them straight to the output TSV
    // (--low-memory). --skip-jaccard short-circuits both, leaving the
    // SparseMatrix unbuilt.
    let input_stem = Path::new(&args.bam)
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("output");
    let output_dir = Path::new(&args.output_dir);

    let jaccard_tuples: Vec<(String, String, f64, f64)> = if args.skip_jaccard {
        log.stage("Skipped Jaccard computation");
        Vec::new()
    } else {
        let mut sparse_matrix = matrix::SparseMatrix::new();
        for group_id in 0..read_ann_matrix.num_groups() {
            sparse_matrix.insert_group(
                group_id,
                read_ann_matrix.row(group_id).iter().map(|(s, _)| *s),
            );
        }
        if stream_jaccard {
            // Streaming path: write directly to file, no tuples buffer.
            std::fs::create_dir_all(output_dir)?;
            output::write_jaccard_streaming(
                output_dir,
                input_stem,
                &sparse_matrix,
                &index.interner,
                args.similarity_threshold,
                args.min_reads,
            )?;
            log.stage("Streamed Jaccard pairs (unsorted)");
            Vec::new()
        } else {
            let tuples = sparse_matrix.jaccard_tuples(
                &index.interner,
                args.similarity_threshold,
                args.min_reads,
            );
            log.stage(&format!(
                "Computed {} Jaccard pairs (threshold {:.2}, min_reads {})",
                tuples.len(),
                args.similarity_threshold,
                args.min_reads,
            ));
            tuples
        }
    };

    // Step 8: Build output directives (timed separately from streaming write)
    let directives = output::build_output_directives(
        &groups,
        &read_ann_matrix,
        &resolutions,
        args.all_alignments,
        args.include_no_feature,
        no_feature_class,
    );
    log.stage(&format!("Built {} output directives", directives.len()));

    // Step 9: Write remaining TSV outputs + stream BAM output (Pass 2). The
    // jaccard.tsv was already produced above (either by the streaming path or
    // by the earlier write phase, depending on mode); skip it here in either
    // of those cases.
    let skip_jaccard_write = args.skip_jaccard || stream_jaccard;
    output::write_outputs(
        &args.bam,
        output_dir,
        input_stem,
        &header,
        &directives,
        &summary_counts,
        &jaccard_tuples,
        &index.interner,
        &index.lengths,
        &prior_mass,
        &posterior_mass,
        args.debug,
        skip_jaccard_write,
    )?;

    // Write run summary
    let peak_rss_bytes = logging::peak_rss_bytes();
    // rayon::current_num_threads reflects whatever the global pool resolved to,
    // including cases where we set it explicitly above and cases where rayon
    // auto-detected (effective_threads == 0).
    let threads_used = rayon::current_num_threads();
    let summary = output::RunSummary {
        total_records: input_order,
        dropped_unmapped,
        num_groups: groups.len(),
        num_dropped_no_feature,
        num_confident: summary_counts.num_confident,
        num_ambiguous: summary_counts.num_ambiguous,
        num_no_feature: summary_counts.num_no_feature,
        validation_failures: pairing_stats.total_failures(),
        confidence: args.confidence,
        em_iterations,
        em_stop_reason,
        elapsed_secs: log.elapsed_secs(),
        peak_rss_bytes,
        threads_used,
    };
    output::write_run_summary(output_dir, input_stem, &summary)?;

    log.stage(&format!("Wrote output to {}/", args.output_dir));
    if let Some(bytes) = peak_rss_bytes {
        log.stage(&format!("Peak RAM: {}", logging::format_bytes(bytes)));
    }
    log.stage(&format!("Done in {:.2}s", log.elapsed_secs()));

    if verbosity != Verbosity::Quiet {
        summary.print_stderr();
    }

    Ok(())
}
