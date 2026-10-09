//! `assign`: a drop-in for `telescope assign` that reproduces the Python
//! reference's report exactly, in a fraction of the memory and time.
//!
//! The pieces are deliberately separable:
//!   * [`annotation`] - GTF loading and overlap rules (`--overlap_compat`);
//!   * [`loader`]     - BAM streaming into a compact score matrix;
//!   * [`model`]      - the reassignment model behind a trait;
//!   * [`report`]     - reassignment modes and the report writer;
//!   * [`updated_sam`] - the annotated BAM of `--updated_sam`;
//!   * [`numpy`]      - numeric primitives matched to numpy/scipy;
//!   * `pyset`/`pytree` - Python set and `intervaltree` ordering, which decide
//!     ties between equally-overlapped loci in Telescope-compatible mode.

pub mod annotation;
mod coverage;
pub mod gtf_ties;
pub mod loader;
pub mod model;
pub mod numpy;
mod pyset;
mod pytree;
mod rawbam;
pub mod report;
pub mod updated_sam;

use std::io;
use std::path::Path;

use clap::{Parser, ValueEnum};

use crate::logging::{Logger, Verbosity, format_bytes, peak_rss_bytes};
use annotation::{Annotation, OverlapCompat, OverlapRules};
use coverage::Coverage;
use loader::{LoadOptions, SamOutputs, Stranded};
use rawbam::{CompressionLevel, RawReader, RawWriter};
use numpy::Mt19937;
use updated_sam::Content;
use model::{ReassignmentModel, TelescopeEm};
use report::{ReassignMode, ReportInputs};

/// Telescope release whose behaviour (and report header) this mode reproduces.
const REFERENCE_VERSION: &str = "1.0.4.1";

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum ReassignModeArg {
    Exclude,
    Choose,
    Average,
    Conf,
    Unique,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum StrandedArg {
    #[value(name = "None")]
    None,
    #[value(name = "RF")]
    Rf,
    #[value(name = "R")]
    R,
    #[value(name = "FR")]
    Fr,
    #[value(name = "F")]
    F,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum OverlapCompatArg {
    /// Reproduce Telescope's overlap arithmetic, quirks included
    Telescope,
    /// True coordinates; no abort on bridging rows
    Corrected,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum UpdatedSamContentArg {
    /// Every alignment of every fragment that overlaps the annotation, plus
    /// `other.bam` for the rest (what Telescope writes)
    All,
    /// Only the alignment each fragment was assigned to
    Assigned,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum TieHashArg {
    /// Python 3.8 and newer
    Python38,
    /// Python 3.7 and older
    Python37,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum MathArg {
    Portable,
    System,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum FloatSumsArg {
    Numpy,
    Sequential,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum ModelArg {
    /// Telescope's EM, bit-for-bit
    Telescope,
    /// Telescope's EM without theta, the extra per-feature weight applied to
    /// ambiguous fragments: plain EM on the proportions. Loading, overlap
    /// rules, reassignment modes and the report are unchanged;
    /// --theta_prior has no effect
    #[value(name = "no-theta")]
    NoTheta,
}

/// Option names and defaults follow `telescope assign`.
#[derive(Parser, Debug)]
pub struct AssignArgs {
    /// Alignment file (SAM or BAM), collated so a fragment's alignments are adjacent
    pub samfile: String,

    /// Annotation file (GTF)
    pub gtffile: String,

    /// GTF attribute that defines a locus
    #[arg(long, default_value = "locus")]
    pub attribute: String,

    /// Name used for alignments that overlap no feature
    #[arg(long = "no_feature_key", default_value = "__no_feature")]
    pub no_feature_key: String,

    /// Output directory
    #[arg(long, default_value = ".")]
    pub outdir: String,

    /// Experiment tag (output file prefix)
    #[arg(long = "exp_tag", default_value = "telescope")]
    pub exp_tag: String,

    /// Reassignment mode used for the final_count column (and for the updated
    /// BAM and bigWig). `choose` draws one candidate with probability equal
    /// to its fitted weight and `average` gives every candidate its fitted
    /// weight; with --legacy both act as in Telescope, only on exact ties
    /// between best hits (uniform draw, even split)
    #[arg(long = "reassign_mode", value_enum, default_value = "exclude")]
    pub reassign_mode: ReassignModeArg,

    /// Minimum probability for a high-confidence assignment; must be above
    /// 0.5 unless --legacy is given
    #[arg(long = "conf_prob", default_value_t = 0.9)]
    pub conf_prob: f64,

    /// Fraction of a fragment that must lie within a feature
    #[arg(long = "overlap_threshold", default_value_t = 0.2)]
    pub overlap_threshold: f64,

    /// Library orientation for strand-aware assignment
    #[arg(long = "stranded_mode", value_enum, default_value = "None")]
    pub stranded_mode: StrandedArg,

    /// Prior on pi (equivalent to adding n unique reads)
    #[arg(long = "pi_prior", default_value_t = 0)]
    pub pi_prior: i64,

    /// Prior on theta (equivalent to adding n non-unique reads)
    #[arg(long = "theta_prior", default_value_t = 200000)]
    pub theta_prior: i64,

    /// EM convergence cutoff
    #[arg(long = "em_epsilon", default_value_t = 1e-7)]
    pub em_epsilon: f64,

    /// EM maximum iterations
    #[arg(long = "max_iter", default_value_t = 100)]
    pub max_iter: u32,

    /// Use the change in log-likelihood as the convergence criterion
    #[arg(long = "use_likelihood")]
    pub use_likelihood: bool,

    /// Reassignment model for ambiguous fragments
    #[arg(long, value_enum, default_value = "telescope")]
    pub model: ModelArg,

    /// How feature overlaps are computed
    #[arg(long = "overlap_compat", value_enum, default_value = "telescope")]
    pub overlap_compat: OverlapCompatArg,

    /// Overlap coordinates only: Telescope's (shifted one base, abort on a
    /// row bridging two intervals of its locus) or corrected
    /// [default: follows --overlap_compat]
    #[arg(long = "overlap_coords", value_enum)]
    pub overlap_coords: Option<OverlapCompatArg>,

    /// Ties between equally-overlapped loci only: Telescope's (Python set
    /// order) or corrected (first locus in the GTF)
    /// [default: follows --overlap_compat]
    #[arg(long = "overlap_ties", value_enum)]
    pub overlap_ties: Option<OverlapCompatArg>,

    /// Which Python's hashing decides Telescope-compatible ties. Python 3.8
    /// changed how tuples hash, so Telescope under 3.7 or older breaks ties
    /// differently from Telescope under 3.8 or newer
    #[arg(long = "tie_hash", value_enum, default_value = "python38")]
    pub tie_hash: TieHashArg,

    /// Where expm1, log1p and log10 come from. `portable`: a correctly
    /// rounded implementation, so output is identical on every machine.
    /// `system`: the platform's math library, which is what numpy uses, so
    /// output matches Python Telescope run on the same machine but can
    /// differ in the last digit between machines
    #[arg(long, value_enum, default_value = "portable")]
    pub math: MathArg,

    /// How floating-point sums are accumulated: as numpy does (needed for
    /// bit-identical results) or plainly left to right
    #[arg(long = "float_sums", value_enum, default_value = "numpy")]
    pub float_sums: FloatSumsArg,

    /// Generate an updated alignment file (<exp_tag>-updated.bam)
    #[arg(long = "updated_sam")]
    pub updated_sam: bool,

    /// With --updated_sam, write exactly the files Telescope writes: every
    /// alignment of every overlapping fragment in <exp_tag>-updated.bam,
    /// <exp_tag>-other.bam for the rest, the intermediate
    /// <exp_tag>-tmp_tele.bam left in place, and compression level 6.
    /// Without it the updated file holds only the alignment each fragment
    /// was assigned to, nothing else is kept, and compression is level 1.
    /// It also selects Telescope's own form of --reassign_mode choose and
    /// average and lifts the --conf_prob limit; the report is otherwise the
    /// same either way
    #[arg(long)]
    pub legacy: bool,

    /// Which alignments the updated file holds [default: assigned, or all
    /// with --legacy]
    #[arg(long = "updated_sam_content", value_enum)]
    pub updated_sam_content: Option<UpdatedSamContentArg>,

    /// Compression level for the output BAMs, 0 (none) to 9. Lower is faster
    /// and larger; the alignments inside are the same [default: 1, or 6
    /// with --legacy]
    #[arg(long = "compression_level", value_parser = clap::value_parser!(u8).range(0..=9))]
    pub compression_level: Option<u8>,

    /// Write coverage of assigned fragments as bigWig: <exp_tag>-coverage.bw,
    /// or -coverage.plus.bw and -coverage.minus.bw when --stranded_mode is
    /// set. Each assigned alignment adds its share of the fragment (1, or the
    /// fitted weight under --reassign_mode average) over the bases it
    /// covers. Independent of --updated_sam; on its own it leaves no BAM
    /// behind
    #[arg(long)]
    pub bigwig: bool,

    /// Do not write <exp_tag>-other.bam (the fragments that are unmapped or
    /// overlap no feature) even when all alignments are requested
    #[arg(long = "no-other", alias = "no_other")]
    pub no_other: bool,

    /// Silence progress output
    #[arg(long)]
    pub quiet: bool,

    /// Print per-iteration detail
    #[arg(long)]
    pub debug: bool,
}

pub fn run_assign(args: AssignArgs) -> io::Result<()> {
    let verbosity = if args.quiet {
        Verbosity::Quiet
    } else if args.debug {
        Verbosity::Verbose
    } else {
        Verbosity::Normal
    };
    let log = Logger::new(verbosity);

    let stranded = match args.stranded_mode {
        StrandedArg::None => Stranded::None,
        StrandedArg::Rf => Stranded::Rf,
        StrandedArg::R => Stranded::R,
        StrandedArg::Fr => Stranded::Fr,
        StrandedArg::F => Stranded::F,
    };
    let legacy = args.legacy;
    if !legacy && args.conf_prob <= 0.5 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "--conf_prob must be above 0.5 (got {}): at 0.5 or below one fragment can pass for \
                 several features at once. Telescope accepts such values; use --legacy to do the same.",
                args.conf_prob
            ),
        ));
    }
    let compat = match args.overlap_compat {
        OverlapCompatArg::Telescope => OverlapCompat::Telescope,
        OverlapCompatArg::Corrected => OverlapCompat::Corrected,
    };
    let pick = |arg: Option<OverlapCompatArg>| match arg {
        Some(OverlapCompatArg::Telescope) => OverlapCompat::Telescope,
        Some(OverlapCompatArg::Corrected) => OverlapCompat::Corrected,
        None => compat,
    };
    let rules = OverlapRules { coords: pick(args.overlap_coords), ties: pick(args.overlap_ties) };
    pyset::set_pre38_tuple_hash(args.tie_hash == TieHashArg::Python37);
    numpy::set_system_math(args.math == MathArg::System);
    numpy::set_sequential_sums(args.float_sums == FloatSumsArg::Sequential);
    let content = match args.updated_sam_content {
        Some(UpdatedSamContentArg::All) => Content::All,
        Some(UpdatedSamContentArg::Assigned) => Content::Assigned,
        None if legacy => Content::All,
        None => Content::Assigned,
    };
    let outfile = |suffix: &str| -> String {
        Path::new(&args.outdir).join(format!("{}-{suffix}", args.exp_tag)).to_string_lossy().into_owned()
    };

    let annot = Annotation::from_gtf(&args.gtffile, &args.attribute, stranded != Stranded::None, rules)?;
    log.stage(&format!("Loaded {} features.", annot.loci.len()));
    log.detail(&format!("Annotation loaded; peak memory so far {}", peak_rss_bytes().map(format_bytes).unwrap_or_default()));

    // Like Telescope, park overlapping fragments in <tag>-tmp_tele.bam while
    // loading; Telescope leaves that file (and other.bam) behind, so the
    // "all" content does too. The "assigned" content keeps only updated.bam.
    let tagged_path = outfile("tmp_tele.bam");
    let level = CompressionLevel::new(args.compression_level.unwrap_or(if legacy { 6 } else { 1 }))
        .expect("level range checked by the parser");
    let reader = RawReader::open(&args.samfile)?;
    let sam_out = if args.updated_sam || args.bigwig {
        Some(SamOutputs {
            // The tagged file is an intermediate, read back once: compress it lightly.
            tagged: RawWriter::create(&tagged_path, &reader.header, CompressionLevel::FAST)?,
            other: match content {
                Content::All if args.updated_sam && !args.no_other => {
                    Some(RawWriter::create(&outfile("other.bam"), &reader.header, level)?)
                }
                Content::All => None,
                Content::Assigned => None,
            },
        })
    } else {
        None
    };
    drop(reader);
    // The intermediate file exists from here on. Telescope leaves it behind
    // however the run ends, and so does --legacy; otherwise it is removed on
    // every exit path, including the early ones below.
    let wrote_intermediate = sam_out.is_some();
    let keep_intermediate = args.updated_sam && legacy;
    let discard_intermediate = || {
        if wrote_intermediate && !keep_intermediate {
            std::fs::remove_file(&tagged_path).ok();
        }
    };
    let loaded = loader::load(
        &args.samfile,
        &annot,
        &LoadOptions {
            overlap_threshold: args.overlap_threshold,
            stranded,
            no_feature_key: args.no_feature_key.clone(),
        },
        sam_out,
    )
    .inspect_err(|_| discard_intermediate())?;
    drop(annot);
    let i = &loaded.info;
    log.stage(&format!(
        "Alignment summary: {} total fragments ({} pairs, {} mixed, {} single, {} unmapped); \
         {} unique, {} multi-mapped; {} overlap annotation ({} one locus, {} several)",
        i.total_fragments,
        i.pair_mapped,
        i.pair_mixed,
        i.single_mapped,
        i.unmapped,
        i.unique,
        i.ambig,
        i.overlap_unique + i.overlap_ambig,
        i.overlap_unique,
        i.overlap_ambig
    ));
    if loaded.tied_alignments > 0 {
        log.detail(&format!(
            "{} alignments overlapped two or more loci equally ({})",
            loaded.tied_alignments,
            match rules.ties {
                OverlapCompat::Telescope => "resolved in Telescope's order",
                OverlapCompat::Corrected => "resolved to the first locus in the GTF",
            }
        ));
    }
    if i.overlap_unique + i.overlap_ambig == 0 {
        log.stage("No alignments overlapping annotation");
        discard_intermediate();
        return Ok(());
    }

    let m = &loaded.matrix;
    // Telescope.get_random_seed()
    let seed = ((i.total_fragments % m.n_rows as u64) * m.n_cols as u64 % 4294967295) as u32;

    let model: Box<dyn ReassignmentModel> = match args.model {
        ModelArg::Telescope | ModelArg::NoTheta => Box::new(TelescopeEm {
            use_theta: args.model == ModelArg::Telescope,
            pi_prior: args.pi_prior,
            theta_prior: args.theta_prior,
            epsilon: args.em_epsilon,
            max_iter: args.max_iter,
            use_likelihood: args.use_likelihood,
        }),
    };
    let em_start = log.elapsed_secs();
    let fit = model.fit(m, &mut |line| log.detail(line));
    log.stage(&format!(
        "EM {} after {} iterations ({:.2}s). Final log-likelihood: {:.6}.",
        if fit.converged { "converged" } else { "terminated" },
        fit.iterations,
        log.elapsed_secs() - em_start,
        fit.log_likelihood
    ));

    // What a dataset exercises: the numpy-style sums only differ from plain
    // ones on long rows (9+ candidates) or long vectors (8192+), and the
    // zero-weight rule only when a candidate underflows.
    let row_len = |r: usize| m.indptr[r + 1] - m.indptr[r];
    log.detail(&format!(
        "Matrix: {} fragments x {} features, {} entries; longest row {}, rows with 9+ candidates {}; \
         candidates dropped at zero weight in the last iteration {}",
        m.n_rows,
        m.n_cols,
        m.data.len(),
        (0..m.n_rows).map(row_len).max().unwrap_or(0),
        (0..m.n_rows).filter(|&r| row_len(r) >= 9).count(),
        fit.absent.iter().filter(|&&a| a).count()
    ));

    // Developer aid for chasing last-bit differences against the reference:
    // TELESCOPE_RS_DUMP_PI=<path> writes feature, pi and pi_init as exact hex floats.
    if let Ok(dump) = std::env::var("TELESCOPE_RS_DUMP_PI") {
        use std::io::Write;
        let mut f = io::BufWriter::new(std::fs::File::create(dump)?);
        for (j, name) in loaded.feat_names.iter().enumerate() {
            writeln!(f, "{name}\t{:016x}\t{:016x}", fit.pi[j].to_bits(), fit.pi_init[j].to_bits())?;
        }
    }

    let mode = match args.reassign_mode {
        ReassignModeArg::Exclude => ReassignMode::Exclude,
        // Telescope's choose and average only act on exact ties between best
        // hits and ignore the fitted weights otherwise; that form is kept
        // for --legacy.
        ReassignModeArg::Choose if legacy => ReassignMode::Choose,
        ReassignModeArg::Choose => ReassignMode::ChooseWeighted,
        ReassignModeArg::Average if legacy => ReassignMode::Average,
        ReassignModeArg::Average => ReassignMode::Fractional,
        ReassignModeArg::Conf => ReassignMode::Conf,
        ReassignModeArg::Unique => ReassignMode::Unique,
    };
    let path = outfile("telescope_report.tsv");
    let mut rng = Mt19937::new(seed);
    let weights = report::Weights { z: &fit.z, absent: &fit.absent };
    // The two modes that are not Telescope's are computed once, from their own
    // generator, and that one result feeds the report, the BAM and the bigWig.
    let own_mode = matches!(mode, ReassignMode::ChooseWeighted | ReassignMode::Fractional);
    let precomputed =
        own_mode.then(|| report::reassign_entries(m, weights, mode, args.conf_prob, &mut Mt19937::new(seed)));
    report::write_report(
        &path,
        &ReportInputs {
            matrix: m,
            fit: &fit,
            feat_names: &loaded.feat_names,
            feat_lengths: &loaded.feat_lengths,
            info: i,
            version: REFERENCE_VERSION,
            mode,
            conf_prob: args.conf_prob,
            final_entries: precomputed.as_deref(),
        },
        &mut rng,
    )?;
    if args.updated_sam || args.bigwig {
        // Telescope recomputes the reassignment here, drawing from the same
        // generator again in `choose` mode.
        let assigned = match precomputed {
            Some(entries) => entries,
            None => report::reassign_entries(m, weights, mode, args.conf_prob, &mut rng),
        };
        let updated = outfile("updated.bam");
        let n_refs = loaded.header.reference_sequences().len();
        let mut coverage = args.bigwig.then(|| Coverage::new(if stranded == Stranded::None { 1 } else { 2 }, n_refs));
        updated_sam::write(
            &tagged_path,
            &loaded,
            &fit,
            &assigned,
            updated_sam::Outputs {
                bam: args.updated_sam.then_some((updated.as_str(), content, level)),
                coverage: coverage.as_mut().map(|c| (c, stranded)),
            },
        )?;
        // Telescope leaves its intermediate file behind; only --legacy does too.
        if !(args.updated_sam && legacy) {
            std::fs::remove_file(&tagged_path)?;
        }
        if args.updated_sam {
            log.stage(&format!("Wrote {updated}"));
        }
        if let Some(cov) = coverage.as_mut() {
            let names: &[&str] = if cov.tracks() == 2 { &["coverage.plus.bw", "coverage.minus.bw"] } else { &["coverage.bw"] };
            for (track, name) in names.iter().enumerate() {
                let path = outfile(name);
                if cov.write(track, &path, &loaded.header)? {
                    log.stage(&format!("Wrote {path}"));
                } else {
                    log.stage(&format!("No assigned fragments for {name}; not written"));
                }
            }
        }
    }
    let rss = peak_rss_bytes().map(format_bytes).unwrap_or_else(|| "n/a".into());
    log.stage(&format!("Wrote {path} (peak memory {rss})"));
    Ok(())
}
