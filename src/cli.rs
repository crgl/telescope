use clap::{Parser, Subcommand as ClapSubcommand, ValueEnum};

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum ModelType {
    Ribbonfish,
    Telescope,
}

/// Whether the EM E-step divides θ by annotation length (RSEM-style).
/// `Auto` resolves to ON for `ribbonfish` and OFF for `telescope`; explicit
/// `On` / `Off` override the model-default.
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum LengthCorrection {
    Auto,
    On,
    Off,
}

/// Whether to partition `__no_feature__` reads into cytogenetic bands.
/// `Auto` resolves to ON for `ribbonfish` and OFF for `telescope`; explicit
/// `On` / `Off` override the model-default. Mirrors `LengthCorrection`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum BandNoFeature {
    Auto,
    On,
    Off,
}

/// Library strandedness mode for `--stranded` on the `annotate` subcommand.
///
/// Controls whether read↔annotation overlaps are filtered by strand:
/// - `RF`: paired, read1 antisense to transcript (read1 strand opposite to annotation strand)
/// - `FR`: paired, read1 sense to transcript (read1 strand same as annotation strand)
/// - `F`: single-end, read sense to transcript (read strand same as annotation strand)
/// - `R`: single-end, read antisense to transcript (read strand opposite to annotation strand)
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum StrandedMode {
    Rf,
    Fr,
    F,
    R,
}

#[derive(Parser)]
#[command(name = "rusty_telescope")]
#[command(about = "Tag BAM reads with GTF annotation overlap")]
pub struct Cli {
    #[command(subcommand)]
    pub command: Subcommand,
}

#[derive(ClapSubcommand)]
pub enum Subcommand {
    /// Annotate BAM reads with GTF feature overlap, run EM, and write tagged
    /// BAM + summary TSVs.
    Annotate(AnnotateArgs),
    /// Sample reads from a BAM, summarize per-annotation strand concordance,
    /// and recommend a `--stranded` mode.
    DetectStrand(DetectStrandArgs),
    /// Legacy mode: reproduce `telescope assign` exactly (same options and
    /// report), with far lower memory use.
    Assign(crate::legacy::AssignArgs),
}

#[derive(Parser, Debug)]
pub struct AnnotateArgs {
    /// Path to input BAM file
    pub bam: String,

    /// Path to GTF annotation file
    #[arg(short, long, default_value = "reference/hg38.hml2_only.gtf")]
    pub gtf: String,

    /// GTF attribute field to use for tagging
    #[arg(short = 'f', long, default_value = "gene_id")]
    pub field: String,

    /// Include reads with no feature overlap in output
    #[arg(long)]
    pub include_no_feature: bool,

    /// Auto-normalize chromosome names (add/remove 'chr' prefix)
    #[arg(long)]
    pub normalize_chr: bool,

    /// Number of threads for annotation processing (0 = auto-detect)
    #[arg(short = 't', long, default_value_t = 0)]
    pub threads: usize,

    /// Treat each alignment independently (no QNAME grouping)
    #[arg(long)]
    pub single_end: bool,

    /// Output directory for results
    #[arg(short = 'o', long, default_value = "tele_out")]
    pub output_dir: String,

    /// Output all alignments, tagging representatives with ZR
    #[arg(long)]
    pub all_alignments: bool,

    /// Confidence threshold for annotation assignment (0.5–1.0).
    /// Only annotations with softmax probability >= this threshold get a ZF tag.
    /// A threshold of 1.0 corresponds to "unique" mode (only unambiguous assignments).
    #[arg(long, default_value_t = 0.9, value_parser = validate_confidence)]
    pub confidence: f64,

    /// Exclude discordant mate pairs (mates mapping to different chromosomes)
    #[arg(long)]
    pub exclude_discordant: bool,

    /// Minimum overlap (bp) between an alignment (or merged mate-pair) and an
    /// annotation for that overlap to be considered valid. Overlaps below this
    /// threshold are dropped before they enter the matrix or EM. Set to 0 to
    /// disable.
    #[arg(long, default_value_t = 30)]
    pub min_overlap: usize,

    /// Minimum overlap coefficient for a pair to be written to the Jaccard
    /// TSV. Applied before any sorting, so it also prunes in-memory tuples.
    #[arg(long, default_value_t = 0.05, value_parser = validate_similarity_threshold)]
    pub similarity_threshold: f64,

    /// Minimum number of read groups overlapping an annotation. Pairs in the
    /// Jaccard TSV are dropped when either annotation has fewer supporting
    /// read groups than this. Set to 0 to disable.
    #[arg(long, default_value_t = 10)]
    pub min_reads: usize,

    /// Skip the BAM streaming write step. All other outputs (summary, jaccard,
    /// run_summary) are still produced. Intended for profiling the compute
    /// phases without the I/O-dominated write cost.
    #[arg(long)]
    pub debug: bool,

    /// Suppress all stderr output (errors only)
    #[arg(long, conflicts_with = "verbose")]
    pub quiet: bool,

    /// Show additional detail at each pipeline stage
    #[arg(long, conflicts_with = "quiet")]
    pub verbose: bool,

    /// Prior-probability model. `ribbonfish` = per-fragment softmax with
    /// temperature 5 (default). `telescope` = linearly rescale AS scores
    /// so the global min maps to 0 and the global max maps to 100, then
    /// normalize via expm1 (global min contributes zero mass).
    #[arg(long, value_enum, default_value_t = ModelType::Ribbonfish)]
    pub model: ModelType,

    /// Maximum EM iterations for read-distribution estimation. EM stops at this
    /// count or once it converges (see `--min-mass-delta`), whichever comes first.
    #[arg(long, default_value_t = 100)]
    pub max_iter: usize,

    /// EM convergence threshold: stop once the largest per-annotation change in
    /// assigned fragment mass (Σγ) between consecutive iterations falls below this
    /// value. Default 1.0 (less than one fragment reassigned on any annotation).
    /// Set to 0 to disable early stopping and always run `--max-iter` iterations.
    #[arg(long, default_value_t = 1.0)]
    pub min_mass_delta: f64,

    /// Telescope-only. Number of ambiguously mapping reads added to each
    /// annotation's θ in the M-step as a regularization measure. Ignored for
    /// `--model ribbonfish`.
    #[arg(long, default_value_t = 200_000)]
    pub theta_prior: u64,

    /// Whether to divide θ by annotation length in the E-step (RSEM-style
    /// effective-length correction). `auto` (default) resolves to `on` for
    /// the `ribbonfish` model and `off` for `telescope`. With correction on,
    /// a fragment with otherwise equal evidence on two annotations favors
    /// the shorter one in proportion to its inverse length.
    #[arg(long, value_enum, default_value_t = LengthCorrection::Auto)]
    pub length_correction: LengthCorrection,

    /// Library strandedness. Drops overlaps where the read's strand is
    /// incompatible with the annotation's GTF strand. Use `detect-strand` to
    /// pick the right mode for an unknown library. When unset, all overlaps
    /// are kept (unstranded).
    #[arg(long, value_enum)]
    pub stranded: Option<StrandedMode>,

    /// Partition `__no_feature__` reads into cytogenetic bands from
    /// `--cytoband`. Each intergenic read is assigned to its best-overlapping
    /// band (e.g. `__no_feature_1p36.33__`), splitting the summary and Jaccard
    /// outputs by band; the run-summary stays aggregated. `auto` (default)
    /// resolves to `on` for `ribbonfish` and `off` for `telescope`.
    #[arg(long, value_enum, default_value_t = BandNoFeature::Auto)]
    pub band_no_feature: BandNoFeature,

    /// Path to a UCSC cytoBand file (tab-separated: chrom, chromStart,
    /// chromEnd, name, gieStain). Loaded only when `--band-no-feature` resolves
    /// to `on`.
    #[arg(long, default_value = "reference/cytoBand.txt")]
    pub cytoband: String,

    /// Skip Jaccard similarity computation. The {stem}_jaccard.tsv file is
    /// not produced. Saves the cost of building the SparseMatrix and the
    /// pairs buffer entirely.
    #[arg(long)]
    pub skip_jaccard: bool,

    /// Reduce peak RAM at the cost of some speed. Streams Jaccard pairs
    /// directly to the TSV (unsorted; no in-memory tuples buffer) and
    /// defaults --threads to 1 (overridable with --threads N). Incompatible
    /// with --all-alignments. Pass --skip-jaccard alongside to also skip
    /// computing Jaccard.
    #[arg(long)]
    pub low_memory: bool,
}

#[derive(Parser, Debug)]
pub struct DetectStrandArgs {
    /// Path to input BAM file
    pub bam: String,

    /// Path to GTF annotation file
    #[arg(short, long, default_value = "reference/hg38.hml2_only.gtf")]
    pub gtf: String,

    /// GTF attribute field to use for grouping rows in the output TSV
    #[arg(short = 'f', long, default_value = "gene_id")]
    pub field: String,

    /// Output directory for results
    #[arg(short = 'o', long, default_value = "tele_out")]
    pub output_dir: String,

    /// Number of annotation-overlapping records to sample before stopping
    #[arg(long, default_value_t = 10_000)]
    pub sample: usize,

    /// Auto-normalize chromosome names (add/remove 'chr' prefix)
    #[arg(long)]
    pub normalize_chr: bool,

    /// Suppress all stderr output (errors only)
    #[arg(long, conflicts_with = "verbose")]
    pub quiet: bool,

    /// Show additional detail
    #[arg(long, conflicts_with = "quiet")]
    pub verbose: bool,
}

fn validate_confidence(s: &str) -> Result<f64, String> {
    let val: f64 = s
        .parse()
        .map_err(|_| format!("'{s}' is not a valid number"))?;
    if (0.5..=1.0).contains(&val) {
        Ok(val)
    } else {
        Err("confidence must be between 0.5 and 1.0".to_string())
    }
}

fn validate_similarity_threshold(s: &str) -> Result<f64, String> {
    let val: f64 = s
        .parse()
        .map_err(|_| format!("'{s}' is not a valid number"))?;
    if (0.0..=1.0).contains(&val) {
        Ok(val)
    } else {
        Err("similarity-threshold must be between 0.0 and 1.0".to_string())
    }
}
