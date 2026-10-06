//! `assign`: a drop-in for `telescope assign` that reproduces the Python
//! reference's report exactly, in a fraction of the memory and time.
//!
//! The pieces are deliberately separable:
//!   * [`annotation`] - GTF loading and overlap rules (`--overlap_compat`);
//!   * [`loader`]     - BAM streaming into a compact score matrix;
//!   * [`model`]      - the reassignment model behind a trait;
//!   * [`report`]     - reassignment modes and the report writer;
//!   * [`numpy`]      - numeric primitives matched to numpy/scipy;
//!   * `pyset`/`pytree` - Python set and `intervaltree` ordering, which decide
//!     ties between equally-overlapped loci in Telescope-compatible mode.

pub mod annotation;
pub mod loader;
pub mod model;
pub mod numpy;
mod pyset;
mod pytree;
pub mod report;

use std::io;
use std::path::Path;

use clap::{Parser, ValueEnum};

use crate::logging::{Logger, Verbosity, format_bytes, peak_rss_bytes};
use annotation::{Annotation, OverlapCompat};
use loader::{LoadOptions, Stranded};
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
pub enum ModelArg {
    /// Telescope's EM, bit-for-bit
    Telescope,
}

/// Option names and defaults follow `telescope assign`.
#[derive(Parser, Debug)]
pub struct AssignArgs {
    /// Alignment file (BAM), collated so a fragment's alignments are adjacent
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

    /// Reassignment mode used for the final_count column
    #[arg(long = "reassign_mode", value_enum, default_value = "exclude")]
    pub reassign_mode: ReassignModeArg,

    /// Minimum probability for a high-confidence assignment
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
    let compat = match args.overlap_compat {
        OverlapCompatArg::Telescope => OverlapCompat::Telescope,
        OverlapCompatArg::Corrected => OverlapCompat::Corrected,
    };

    let annot = Annotation::from_gtf(&args.gtffile, &args.attribute, stranded != Stranded::None, compat)?;
    log.stage(&format!("Loaded {} features.", annot.loci.len()));

    let loaded = loader::load(
        &args.samfile,
        &annot,
        &LoadOptions {
            overlap_threshold: args.overlap_threshold,
            stranded,
            no_feature_key: args.no_feature_key.clone(),
        },
    )?;
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
            match compat {
                OverlapCompat::Telescope => "resolved in Telescope's order",
                OverlapCompat::Corrected => "resolved to the first locus in the GTF",
            }
        ));
    }
    if i.overlap_unique + i.overlap_ambig == 0 {
        log.stage("No alignments overlapping annotation");
        return Ok(());
    }

    let m = &loaded.matrix;
    // Telescope.get_random_seed()
    let seed = ((i.total_fragments % m.n_rows as u64) * m.n_cols as u64 % 4294967295) as u32;

    let model: Box<dyn ReassignmentModel> = match args.model {
        ModelArg::Telescope => Box::new(TelescopeEm {
            pi_prior: args.pi_prior,
            theta_prior: args.theta_prior,
            epsilon: args.em_epsilon,
            max_iter: args.max_iter,
            use_likelihood: args.use_likelihood,
        }),
    };
    let fit = model.fit(m, &mut |line| log.detail(line));
    log.stage(&format!(
        "EM {} after {} iterations. Final log-likelihood: {:.6}.",
        if fit.converged { "converged" } else { "terminated" },
        fit.iterations,
        fit.log_likelihood
    ));

    // Developer aid for chasing last-bit differences against the reference:
    // RUSTY_TELESCOPE_DUMP_PI=<path> writes feature, pi and pi_init as exact hex floats.
    if let Ok(dump) = std::env::var("RUSTY_TELESCOPE_DUMP_PI") {
        use std::io::Write;
        let mut f = io::BufWriter::new(std::fs::File::create(dump)?);
        for (j, name) in loaded.feat_names.iter().enumerate() {
            writeln!(f, "{name}\t{:016x}\t{:016x}", fit.pi[j].to_bits(), fit.pi_init[j].to_bits())?;
        }
    }

    let mode = match args.reassign_mode {
        ReassignModeArg::Exclude => ReassignMode::Exclude,
        ReassignModeArg::Choose => ReassignMode::Choose,
        ReassignModeArg::Average => ReassignMode::Average,
        ReassignModeArg::Conf => ReassignMode::Conf,
        ReassignModeArg::Unique => ReassignMode::Unique,
    };
    let path = Path::new(&args.outdir).join(format!("{}-telescope_report.tsv", args.exp_tag));
    report::write_report(
        &path.to_string_lossy(),
        &ReportInputs {
            matrix: m,
            fit: &fit,
            feat_names: &loaded.feat_names,
            feat_lengths: &loaded.feat_lengths,
            info: i,
            version: REFERENCE_VERSION,
            mode,
            conf_prob: args.conf_prob,
            seed,
        },
    )?;
    let rss = peak_rss_bytes().map(format_bytes).unwrap_or_else(|| "n/a".into());
    log.stage(&format!("Wrote {} (peak memory {rss})", path.display()));
    Ok(())
}
