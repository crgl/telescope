# telescope_rs

A fast Rust toolkit for locus-level quantification of transposable elements from RNA-seq alignments. It has two faces:

- **`assign`** is a drop-in replacement for [Telescope](https://github.com/mlbendall/telescope)'s `telescope assign`. It reproduces the Python program's report byte for byte, at roughly a tenth of the time and memory, and with `--legacy` writes the same alignment files too. See [`assign`: Telescope-compatible mode](#assign-telescope-compatible-mode).
- **`annotate`** is a separate, multithreaded pipeline with its own overlap rules, prior models, confidence filtering and Jaccard similarity output. **`detect-strand`** is its helper for inferring library strandedness. Everything from [Options](#options) onward describes these two.

## Usage

```
telescope_rs <COMMAND> [OPTIONS]
```

| Command | Purpose |
|---------|---------|
| `assign` | Telescope-compatible reassignment: same options and report as `telescope assign`, optionally the same updated BAM |
| `gtf-ties` | List the regions of a GTF that two or more loci share, where `assign` picks one arbitrarily, and which locus wins |
| `annotate` | Tag a BAM with GTF overlap, run EM, write tagged BAM + summary TSVs |
| `detect-strand` | Sample reads, summarize strand concordance per annotation, recommend a `--stranded` mode |

## `assign`: Telescope-compatible mode

> **Reference version.** `assign` is built to match **Telescope v1.0.4.1 from the [hanalysis fork](https://github.com/hanalysis/telescope)** (commit `bb58c43`, the Python code in this repository's `telescope/` directory), run under Python 3.8 or newer. That fork differs from the original `mlbendall/telescope` in adding `--stranded_mode` and in how it handles unpaired reads within paired-end data, so results are not claimed to match other Telescope releases.

```bash
# Same call shape as `telescope assign`
telescope_rs assign alignments.bam annotation.gtf --outdir results --exp_tag sample1

# Also write a BAM of the alignments fragments were assigned to
telescope_rs assign alignments.bam annotation.gtf --updated_sam

# Stranded paired-end library (dUTP)
telescope_rs assign alignments.bam annotation.gtf --stranded_mode RF

# Match a Telescope install that runs on Python 3.7 or older
telescope_rs assign alignments.bam annotation.gtf --tie_hash python37

# Exactly the alignment files Telescope writes (all alignments, other.bam, tmp_tele.bam)
telescope_rs assign alignments.bam annotation.gtf --updated_sam --legacy

# Coverage of assigned fragments as bigWig, with no BAM kept
telescope_rs assign alignments.bam annotation.gtf --bigwig
```

Input is a SAM or BAM file (detected from its content) in which all alignments of a fragment are adjacent, as aligners write them. The run is single-threaded. Output is `<outdir>/<exp_tag>-telescope_report.tsv`, in Telescope's format, with Telescope's version string (`1.0.4.1`) in the header so files compare equal.

### Options shared with `telescope assign`

Names and defaults are Telescope's, underscores included.

| Flag | Default | Description |
|------|---------|-------------|
| `<SAMFILE>` | *(required)* | Alignment file, SAM or BAM |
| `<GTFFILE>` | *(required)* | Annotation file |
| `--attribute` | `locus` | GTF attribute that defines a locus |
| `--no_feature_key` | `__no_feature` | Name for alignments that overlap no feature |
| `--outdir` | `.` | Output directory |
| `--exp_tag` | `telescope` | Output file prefix |
| `--reassign_mode` | `exclude` | Mode behind the `final_count` column: `exclude`, `choose`, `average`, `conf`, `unique` |
| `--conf_prob` | `0.9` | Threshold for the `conf` mode and the `final_conf` column |
| `--overlap_threshold` | `0.2` | Fraction of a fragment that must lie within a feature |
| `--stranded_mode` | `None` | `None`, `RF`, `FR`, `R`, `F` |
| `--pi_prior` | `0` | Prior on pi |
| `--theta_prior` | `200000` | Prior on theta |
| `--em_epsilon` | `1e-7` | EM convergence cutoff |
| `--max_iter` | `100` | EM iteration cap |
| `--use_likelihood` | off | Converge on the change in log-likelihood |
| `--updated_sam` | off | Write the updated BAM (see below). Telescope's own file set needs `--legacy` |
| `--quiet`, `--debug` | off | Less or more log output |

Not available: `--ncpu`, `--tempdir`, `--logfile`, `--skip_em`, `--annotation_class`, `--overlap_mode`, the checkpoint file and `telescope resume`. Abbreviated option names are not accepted.

### Options beyond Telescope

| Flag | Default | Description |
|------|---------|-------------|
| `--legacy` | off | With `--updated_sam`, write exactly Telescope's files: all alignments, `other.bam`, `tmp_tele.bam` kept, compression level 6. Sets the defaults of the options marked *. The report is the same with or without it |
| `--model` | `telescope` | `telescope`: Telescope's EM. `no-theta`: the same EM without theta and its prior |
| `--overlap_compat` | `telescope` | Sets `--overlap_coords` and `--overlap_ties` together |
| `--overlap_coords` | follows `--overlap_compat` | `telescope`: overlaps measured one base to the right, and a GTF row bridging two intervals of its own locus aborts the run. `corrected`: true coordinates, no abort |
| `--overlap_ties` | follows `--overlap_compat` | Which locus wins when two overlap an alignment equally. `telescope`: the order Python's set iteration yields. `corrected`: the locus that appears first in the GTF |
| `--tie_hash` | `python38` | Which Python's tuple hash drives Telescope-style ties: `python38` (3.8 and newer) or `python37` (3.7 and older) |
| `--float_sums` | `numpy` | `numpy`: sums accumulated as numpy and scipy do, needed for bit-identical results. `sequential`: plain left to right |
| `--updated_sam_content` * | `assigned` (`all` with `--legacy`) | `assigned`: only the alignment each fragment was assigned to. `all`: what Telescope writes, including `other.bam` |
| `--compression_level` * | `1` (`6` with `--legacy`) | Compression of the output BAMs, 0 to 9. Lower is faster and larger |
| `--no-other` | off | Do not write `other.bam` even when all alignments are requested |
| `--bigwig` | off | Write coverage of assigned fragments as bigWig (see below) |

### Updated BAM (`--updated_sam`)

By default only `<exp_tag>-updated.bam` is written. It holds one alignment (both mates) per assigned fragment: fragments assigned to the no-feature key or to nothing are left out. Each record carries Telescope's tags: `ZF` is the alignment's feature, `ZT` is `PRI`, `ZB` lists the top-scoring feature(s), `XP` is the membership weight as a percentage, `YC` is a display colour, and MAPQ is the phred-scaled weight. These records are the same ones, with the same content, that Telescope writes for assigned fragments.

With `--legacy` (or `--updated_sam_content all`) three files are written, as Telescope does:

- `<exp_tag>-updated.bam`: every alignment of every fragment that overlaps the annotation. `ZT` is `PRI` for the best alignment per feature and `SEC` otherwise, and alignments the fragment was not assigned to are flagged secondary.
- `<exp_tag>-other.bam`: fragments that are unmapped or overlap no feature. This is most of the input and most of the run time; `--no-other` skips it.
- `<exp_tag>-tmp_tele.bam`: Telescope's intermediate file, which it leaves behind (`--legacy` only).

These match Telescope's files record for record as SAM text. The BAM bytes differ, since compression is not reproduced.

Telescope's code tries to add an `@PG` line to the updated BAM but the append has no effect, so the header equals the input's. This is reproduced.

On a 40-million-pair sample with a HERV annotation, `--updated_sam` takes about 29 s by default and about 137 s with `--legacy`; Python Telescope takes about 580 s.

### Coverage (`--bigwig`)

`--bigwig` writes per-base coverage of assigned fragments to `<exp_tag>-coverage.bw`, or to `<exp_tag>-coverage.plus.bw` and `<exp_tag>-coverage.minus.bw` when `--stranded_mode` is set. Each assigned fragment counts once over the reference bases its assigned alignment covers: aligned blocks only, so introns and deletions are excluded, and the overlap between two mates is not counted twice. Values are raw depth. It can be used with or without `--updated_sam`; on its own it leaves no BAM behind.

### How it relates to Python Telescope

What is reproduced, and why it matters:

- **Overlap arithmetic.** Every GTF row is used regardless of feature type, and overlaps are measured one base to the right of the true position.
- **Ties.** When two loci overlap an alignment equally, Telescope keeps whichever its interval tree returns first, which follows Python's hash-table layout. `assign` carries a port of CPython's `set` and of the `intervaltree` package to return the same order. Ties are common (hundreds of thousands of alignments per sample on a genes-plus-HERV annotation) and decide 1-3% of final counts.
- **Python version.** Python changed its tuple hash in 3.8, so Telescope under Python 3.7 and Telescope under Python 3.10 give different results on the same input (0.9-5.4% of counts in the runs compared). `--tie_hash` selects which one to match; both are reproduced byte for byte.
- **Arithmetic order.** Sums follow numpy's blocked pairwise scheme and scipy's sparse-matrix conventions, including dropping a candidate whose weight underflows to zero. Switching this off (`--float_sums sequential`) changed at most three counts in any run tested.
- **Failures.** Telescope crashes on an empty alignment file, on a mapped alignment with no `AS` tag, and on a GTF row missing the locus attribute or bridging two intervals of its locus. `assign` stops with an error in the same situations. It also stops if a rescaled score would not fit Telescope's 16-bit score matrix, rather than guess what Telescope would do.

Read names are not stored (a 128-bit hash stands in for Telescope's name-to-row dictionary), and only fragments that overlap the annotation are kept in memory.

### Validation

Compared against Python Telescope 1.0.4.1 (Python 3.10, numpy 1.26.4, scipy 1.15.2, macOS arm64). "Identical" means the report matches byte for byte and the fitted proportions match bit for bit; for `--updated_sam`, all three BAMs match as SAM text.

| Data | Comparisons | Result |
|------|-------------|--------|
| Telescope's bundled test data, all five reassign modes, BAM and SAM input | report + BAMs | identical |
| 10 CCLE RNA-seq runs (10% subsamples, HISAT2 `-k 100`), hg38 and T2T-CHM13, HERV-only and HERV+genes annotations, plus an HML2-only annotation | 50 reports | identical |
| The same 50 with `--updated_sam --legacy` | report + BAMs | identical |
| One run at full depth (40 million pairs), HERV annotation, `--updated_sam --legacy` | report + BAMs | identical |
| bowtie2 and STAR alignments of one run | report + BAMs | identical |
| Stranded paired-end library (ENCODE ENCSR000CON), all five `--stranded_mode` values | report + BAMs | identical |
| Single-end (read 1 only) unstranded and stranded, HISAT2 and bowtie2 | report + BAMs | identical |
| PacBio long reads (minimap2, 70,282 supplementary alignments) | report + BAMs | identical |
| Telescope under Python 3.7, with `--tie_hash python37` | 7 reports | identical |

Against the Python 3.7 environment (older numpy and scipy, run under Rosetta) a few dozen records per updated BAM differ in MAPQ only, 160 versus 255, where a membership weight lands exactly on 1 or one rounding step below it.

On the heaviest run (about 230 million alignment records, HERV+genes annotation) Python Telescope took 63 minutes and 25.3 GB; `assign` took 2.5 minutes and 3.4 GB. The other HERV+genes runs were 4-24 minutes and 4-10 GB against 17-44 seconds and 0.5-1 GB.

Not yet covered: coordinate-sorted input, input from Linux or other numpy/scipy versions, and full-depth (unsubsampled) data.

### Regression tests

`cargo test` checks `assign` against the reference without needing Python or any external data:

- `tests/legacy.rs` uses the test data Telescope ships.
- `tests/legacy_fixtures.rs` runs 23 cases over small real-data fixtures in `tests/data/legacy/` (about 12 MB), each with the report Python Telescope produced. They cover paired and single-end reads, stranded and unstranded libraries, HISAT2, bowtie2 and STAR, long reads with supplementary alignments, SAM input, a GTF with mangled quoting, every reassign mode, `--theta_prior 0`, and Telescope under both Python 3.7 and 3.10. One fixture is sized so that numpy's summation order and the zero-weight rule change the answer if they are not reproduced.

`tests/data/legacy/README.md` says what each fixture is for and how it was made.

## `gtf-ties`: where an annotation forces arbitrary choices

Wherever two or more loci cover the same bases, a read lying wholly inside the shared stretch overlaps them all equally, and Telescope (hence `assign`) gives it to whichever locus its interval tree yields first. That order follows Python's hash-table layout, not anything biological. `gtf-ties` lists every such stretch for a GTF, without needing any reads.

```bash
telescope_rs gtf-ties annotation.gtf                       # loci defined by gene_id, Python 3.8+ rules
telescope_rs gtf-ties annotation.gtf --attribute locus     # match an `assign` run that uses the default attribute
telescope_rs gtf-ties annotation.gtf --tie_hash python37   # Telescope under Python 3.7 or older
```

| Flag | Default | Description |
|------|---------|-------------|
| `<GTFFILE>` | *(required)* | Annotation file |
| `--attribute` | `gene_id` | GTF attribute that defines a locus. `assign` defaults to `locus`; use the same value for both |
| `--tie_hash` | `python38` | `python38`: Telescope under Python 3.8 or newer. `python37`: under 3.7 or older |
| `--same_strand` | off | Only report competition between loci on the same strand, as in an `assign` run with a `--stranded_mode` |
| `--outdir` | `.` | Output directory |
| `--exp_tag` | GTF file name | Output file prefix |

Two files are written:

- `<tag>-tie_regions.tsv`: one row per shared stretch, from one interval boundary to the next. Columns: `chrom`, `start`, `end` (GTF coordinates, inclusive), `length`, `winner`, `winner_strand`, `n_competing`, `losers` (each with its strand), and `winner_at_first_base`.
- `<tag>-tie_loci.tsv`: one row per locus involved, with its length, the bases it shares, the bases it wins and loses, and which loci it wins over or loses to (with base counts).

How to read the winner:

- `winner` is the locus that gets a read lying inside the stretch.
- `winner_at_first_base` is the locus that gets a read starting exactly on the stretch's first base. It usually equals `winner`, but Telescope's order also depends on the interval boundaries a read spans, so it can differ.
- A read that crosses into a neighbouring stretch can likewise resolve differently from both; this is not tabulated.
- A read extending beyond the shared bases overlaps the loci unequally and goes to the larger overlap, as usual.
- A locus that is shared along its whole length and never wins receives no read that lies wholly inside it. The run log counts these.

## `annotate` and `detect-strand`

The rest of this document covers the `annotate` pipeline and its `detect-strand` helper. None of it applies to `assign`.

### Examples

```bash
# Tag reads using default GTF and gene_id field
telescope_rs annotate input.bam

# Use a custom GTF and output directory
telescope_rs annotate input.bam --gtf annotations.gtf -o results/

# Keep unmatched reads in output
telescope_rs annotate input.bam --include-no-feature

# Tag with a different GTF attribute
telescope_rs annotate input.bam --field locus

# Single-end mode (no mate pairing, but still groups by QNAME)
telescope_rs annotate input.bam --single-end

# Only include uniquely mapped read pairs
telescope_rs annotate input.bam --confidence 1.0

# Loosen the confidence threshold
telescope_rs annotate input.bam --confidence 0.7

# Output all alignments with representative tags
telescope_rs annotate input.bam --all-alignments

# Drop mate pairs mapping to different chromosomes
telescope_rs annotate input.bam --exclude-discordant

# Use the telescope prior model (global linear rescale of AS scores)
telescope_rs annotate input.bam --model telescope

# Telescope with a lighter θ regularizer (default --theta-prior is 200000)
telescope_rs annotate input.bam --model telescope --theta-prior 1000

# Skip EM (use priors as posteriors)
telescope_rs annotate input.bam --max-iter 0

# Lower the overlap cutoff written to the TSV
telescope_rs annotate input.bam --similarity-threshold 0.01

# Limit to 4 threads
telescope_rs annotate input.bam --threads 4

# Run silently (errors only)
telescope_rs annotate input.bam --quiet

# Show detailed progress
telescope_rs annotate input.bam --verbose

# Restrict overlaps to strand-compatible features (paired, dUTP/TruSeq stranded)
telescope_rs annotate input.bam --stranded RF

# Partition no-feature reads into cytobands (default-on under ribbonfish)
telescope_rs annotate input.bam --cytoband reference/cytoBand.txt

# Disable banding to get the single __no_feature__ bucket back
telescope_rs annotate input.bam --band-no-feature off

# Infer the library strandedness from a sample of reads, then annotate accordingly
telescope_rs detect-strand input.bam
telescope_rs annotate input.bam --stranded "$(grep -E '^# recommended_mode' tele_out/input_strandedness.tsv | cut -d' ' -f3)"
```

## Options

### `annotate`

| Flag | Short | Default | Description |
|------|-------|---------|-------------|
| `<BAM>` | | *(required)* | Path to input BAM file |
| `--gtf` | `-g` | `reference/hg38.hml2_only.gtf` | Path to GTF annotation file |
| `--field` | `-f` | `gene_id` | GTF attribute field to use for tagging |
| `--output-dir` | `-o` | `tele_out` | Output directory for results |
| `--confidence` | | `0.9` | Confidence threshold for annotation assignment (0.5--1.0). A threshold of 1.0 corresponds to "unique" mode (only unambiguous assignments) |
| `--model` | | `ribbonfish` | Prior-probability model: `ribbonfish` (per-fragment softmax) or `telescope` (global linear rescale + expm1). See [Prior Models](#prior-models) |
| `--max-iter` | | `100` | Maximum EM iterations for read-distribution estimation. EM stops at this count or once it converges (`--min-mass-delta`), whichever comes first. `0` skips EM and uses the priors directly |
| `--min-mass-delta` | | `1.0` | EM convergence threshold: stop early once the largest per-annotation change in assigned fragment mass (Σγ) between consecutive iterations falls below this value (i.e. less than one fragment reassigned anywhere). `0` disables early stopping and always runs `--max-iter` iterations |
| `--length-correction` | | `auto` | Whether the EM E-step divides θ by annotation length (RSEM-style). Supported by both models. `auto` resolves to `on` for `ribbonfish` and `off` for `telescope`; explicit `on` / `off` override |
| `--theta-prior` | | `200000` | **Telescope only.** Number of ambiguously mapping reads added to each contended annotation's θ in the M-step as regularization. Ignored under `--model ribbonfish`. See [EM Reassignment](#em-reassignment) |
| `--min-overlap` | | `30` | Minimum combined overlap (bp) between an alignment unit and an annotation for that overlap to be considered valid. Set to 0 to disable |
| `--similarity-threshold` | | `0.05` | Minimum overlap coefficient retained in the Jaccard TSV (applied before any sort) |
| `--min-reads` | | `10` | Minimum number of read groups overlapping an annotation. Pairs in the Jaccard TSV are dropped when either annotation has fewer supporting read groups. Set to 0 to disable |
| `--single-end` | | off | Skip mate pairing (each BAM record is its own alignment unit, still grouped by QNAME) |
| `--all-alignments` | | off | Output all alignments, tagging representatives with ZR |
| `--include-no-feature` | | off | Include reads with no feature overlap in output |
| `--exclude-discordant` | | off | Exclude discordant mate pairs (mates mapping to different chromosomes) |
| `--normalize-chr` | | off | Auto-add/remove `chr` prefix when matching chromosome names between BAM and GTF |
| `--stranded` | | unset (unstranded) | Library strandedness: `RF`, `FR`, `F`, or `R`. See [Stranded Mode](#stranded-mode) |
| `--band-no-feature` | | `auto` | Partition `__no_feature__` reads into cytogenetic bands from `--cytoband`. `auto` resolves to `on` for `ribbonfish` and `off` for `telescope`; explicit `on` / `off` override. See [Band Mode](#band-mode-no-feature-partitioning) |
| `--cytoband` | | `reference/cytoBand.txt` | Path to a UCSC cytoBand file (chrom, chromStart, chromEnd, name, gieStain). Loaded only when `--band-no-feature` resolves to `on` |
| `--skip-jaccard` | | off | Skip Jaccard similarity computation. No `{stem}_jaccard.tsv` is produced. The biggest peak-RAM win on datasets with dense annotation co-occurrence |
| `--low-memory` | | off | Trade speed for peak RAM. Streams Jaccard pairs to the TSV unsorted (no in-memory tuples buffer, no sort) and defaults `--threads` to 1 (overridable via `--threads N`). Combine with `--skip-jaccard` to also skip computing Jaccard. Incompatible with `--all-alignments` |
| `--threads` | `-t` | `0` (auto) | Number of threads for parallel processing |
| `--debug` | | off | Skip the BAM streaming write step (summary, jaccard, run_summary still produced). Intended for profiling the compute phases |
| `--quiet` | | off | Suppress all stderr output (errors only) |
| `--verbose` | | off | Show additional detail at each pipeline stage |

`--quiet` and `--verbose` are mutually exclusive.

### `detect-strand`

| Flag | Short | Default | Description |
|------|-------|---------|-------------|
| `<BAM>` | | *(required)* | Path to input BAM file |
| `--gtf` | `-g` | `reference/hg38.hml2_only.gtf` | Path to GTF annotation file |
| `--field` | `-f` | `gene_id` | GTF attribute field to group rows by in the output TSV |
| `--output-dir` | `-o` | `tele_out` | Output directory for results |
| `--sample` | | `10000` | Number of annotation-overlapping records to sample before stopping |
| `--normalize-chr` | | off | Auto-add/remove `chr` prefix when matching chromosome names between BAM and GTF |
| `--quiet` | | off | Suppress all stderr output (errors only) |
| `--verbose` | | off | Show additional detail |

## Paired-End vs Single-End Mode

By default, telescope_rs operates in **paired-end mode**:

- Adjacent BAM records with the same QNAME and complementary segment flags (0x40 first-in-pair / 0x80 second-in-pair) are paired into a single alignment unit
- Mate pairs are validated per the SAM spec: flag correspondence (0x10/0x20) and RNEXT/PNEXT matching. Validation failures produce a warning and the mates are treated as unpaired
- The alignment spans of paired mates are merged before annotation overlap is computed; overlapping regions between the two mates are not double-counted (inclusion-exclusion)
- The alignment score (AS tag) of a mate pair is the sum of both mates' scores
- Both BAM records of a paired mate receive the same output tags (ZB, ZF, ZR)
- Discordant pairs (mates on different chromosomes) are processed normally but can be excluded with `--exclude-discordant`

In **single-end mode** (`--single-end`):

- No mate pairing is attempted; each BAM record is its own alignment unit
- Records are still grouped by QNAME (multimapped reads share a QNAME and compete during ambiguity resolution)

## Prior Models

The `--model` flag selects how per-fragment alignment scores are converted into prior probabilities over a fragment's candidate annotations. Priors are fixed after this step and fed into EM.

### `ribbonfish` (default)

Per-fragment softmax with a temperature of 5:

$$\pi_{r,a} = \frac{\exp((\text{AS}_{r,a} - \max_{a'} \text{AS}_{r,a'}) / 5)}{\sum_{a'} \exp((\text{AS}_{r,a'} - \max_{a'} \text{AS}_{r,a'}) / 5)}$$

Only each fragment's own AS scores matter; no information is shared across fragments at the prior stage.

### `telescope`

Every fragment's representative AS scores are **linearly rescaled** using the dataset-wide min/max so the global min maps to `0` and the global max maps to `100`. The rescaled scores are passed through `expm1` and normalized:

$$s'_{r,a} = \frac{\text{AS}_{r,a} - \text{AS}_{\min}}{\text{AS}_{\max} - \text{AS}_{\min}} \cdot 100 \qquad \pi_{r,a} \propto \text{expm1}(s'_{r,a})$$

An alignment tied with the global AS minimum contributes zero prior mass (`expm1(0) = 0`). If `AS_min == AS_max` the priors fall back to uniform.

## EM Reassignment

After priors are computed, telescope_rs runs expectation-maximization to share abundance information across fragments. The update depends on the selected `--model`.

EM runs until the first of two stopping conditions: it reaches `--max-iter` iterations, or it converges — the largest per-annotation change in assigned fragment mass (Σγ, the unweighted expected fragment count) between consecutive iterations drops below `--min-mass-delta` (default `1.0`, i.e. less than one fragment reassigned on any annotation). Setting `--min-mass-delta 0` disables early stopping and always runs the full `--max-iter`. Setting `--max-iter 0` returns γ = π, i.e. the confidence pick operates directly on the priors. The final per-fragment γ replaces the prior as the probability used by the confidence threshold. The actual number of iterations run, and why EM stopped (`converged` / `max-iter reached` / `skipped`), are reported in `*_run_summary.txt`.

### `ribbonfish` EM

Let $\pi_{r,a}$ be the fixed prior and $\theta_a$ be the current estimated read count for annotation $a$.

- **Init**: each non-empty fragment splits its `1.0` mass uniformly across its own candidate annotations.
- **E-step** (parallel over fragments): $\gamma_{r,a} = \theta_a \cdot \pi_{r,a} \,/\, \sum_{a' \in R_r} \theta_{a'} \cdot \pi_{r,a'}$
- **M-step** (parallel reduce): $\theta_a = \sum_r \gamma_{r,a}$

### `telescope` EM (π / θ split + per-fragment weights)

The telescope model splits the single read-distribution parameter into two:

- $\pi_a$ — the **total** fraction of reads from annotation $a$ (unique + ambiguous).
- $\theta_a$ — the fraction of **ambiguous** reads from annotation $a$.

Read groups are partitioned by their candidate count: a group with exactly one candidate annotation is **unique**, two or more is **ambiguous** (empty rows pass through unchanged). Unique groups are never iterated — they are not reassigned and their contribution to $\pi$ is fixed — which also makes each EM step cheaper.

Each fragment carries a **weight** $w_r = \mathrm{expm1}\!\big(\frac{\max_a \text{AS}_{r,a} - \text{AS}_{\min}}{\text{AS}_{\max}-\text{AS}_{\min}}\cdot 100\big)$ — the *best* alignment score transformed exactly like the telescope prior, but taken as the max and **before** normalization. The weight is therefore driven only by a fragment's best-scoring alignment and is not diluted by having several valid alignments. With no usable global AS bounds the weight falls back to `1.0` (unweighted); a fragment whose best score is the global minimum gets weight `0`.

Let $\pi^{\text{fix}}_a = \sum_{r \in \text{unique}(a)} w_r$ (constant across iterations).

- **E-step** (ambiguous fragments only): $\gamma_{r,a} \propto \pi_{r,a} \cdot \pi_a \cdot \theta_a$, normalized per fragment. The prior is multiplied by **both** π and θ.
- **M-step**: $\theta_a \propto \texttt{theta\_prior}\cdot\mathrm{expm1}(100) + \sum_{r \in \text{amb}} w_r\,\gamma_{r,a}$ and $\pi_a \propto \pi^{\text{fix}}_a + \sum_{r \in \text{amb}} w_r\,\gamma_{r,a}$, each normalized.

`--theta-prior` (default `200000`) adds that many pseudo-reads to **each contended annotation's** θ in every M-step, a regularizer that damps winner-take-all collapse. Note π is *not* regularized. Each regularizing read is treated as though it aligned **perfectly** — weighted by the maximum weight any fragment can take for the given AS bounds rather than as a raw unit count. Normally that maximum is $\mathrm{expm1}(100)$ (a fragment whose best score equals the global AS max), so even a small `--theta-prior` strongly smooths θ relative to the dataset's real weighted ambiguous mass. In the **degenerate case** where every read shares one AS value (global min == max), genuine fragment weights all fall back to `1.0`, and so does the regularizer's per-read weight — the `theta_prior` reads are then weighted *equally* to genuine ambiguous reads (not over-weighted by $\mathrm{expm1}(100)$).

### Length correction (`--length-correction`)

When length correction is **on**, the E-step substitutes $\theta_a / L_a$ for $\theta_a$ — RSEM-style effective-length normalization. Mechanically, given a fragment with otherwise equal evidence on two candidate annotations, the shorter annotation receives posterior mass in inverse proportion to its length: half the length → twice the γ.

The synthetic `__no_feature__` annotation is assigned a fixed length of $10^6$ bp so it competes meaningfully against typical real (much shorter) annotations rather than dominating them.

Default behavior depends on the prior model: `auto` resolves to **on** for `ribbonfish` and **off** for `telescope`. Pass `--length-correction on` or `off` to override. Length is the union of GTF intervals tagged with the attribute value (the `length` column in `summary.tsv`).

Both EM models support length correction. In the `telescope` model the substitution applies to the read-distribution term θ in the E-step (`prior · π_a · θ_a/L_a`); π and the per-fragment prior are left uncorrected, exactly as in `ribbonfish`. `telescope` defaults to **off** (`auto`); pass `--length-correction on` to enable it.

## Band Mode (no-feature partitioning)

By default a read overlapping no GTF feature collapses into a single synthetic
`__no_feature__` annotation. `--band-no-feature` instead **partitions** those
intergenic ("out of annotation") reads into cytogenetic bands read from a UCSC
`cytoBand.txt` file (`--cytoband`, default `reference/cytoBand.txt`), so the
summary and Jaccard outputs resolve at the band level.

- **Default depends on the model**: `auto` (the default) resolves to **on** for
  `--model ribbonfish` and **off** for `--model telescope`. Pass
  `--band-no-feature on` / `off` to override. With banding off, the legacy single
  `__no_feature__` bucket is restored.
- **Assignment**: each no-feature read (or mate pair) is assigned to the single
  band with the greatest overlap — one band per read. Reads that hit *any* GTF
  feature are unaffected; only the no-feature path is partitioned.
- **Banding never changes which reads are processed.** A read group is kept only
  when at least one of its alignments overlaps a *real* GTF feature; bands do not
  count toward that decision. So banding only **relabels** the `__no_feature__`
  reads that already survive (those sharing a group with a real-feature
  alignment) — read groups with no real-feature overlap are dropped exactly as
  with `--band-no-feature off`. The number of reads overlapping the annotation
  (and the dropped count) is therefore invariant to `--band-no-feature` and
  `--model`.
- **Naming**: bands are named `__no_feature_<chrom><band>__`, combining column 1
  (chromosome) and column 4 (band) of the cytoBand file — e.g. `chr1` + `p36.33`
  → `__no_feature_1p36.33__`. Bands are **always chr-normalized** (independent of
  `--normalize-chr`): `chr1` and `1` are the same chromosome, and the `chr`
  prefix is never included in the band name.
- **Length**: each band reports its real width (`chromEnd - chromStart`) in the
  `length` column, so length-corrected EM (`ribbonfish` default) uses the band's
  true length instead of the arbitrary `1000000`.
- **Outputs**: `summary.tsv` and `jaccard.tsv` split by band (band names appear
  as ordinary annotations); `run_summary.txt` keeps a **single aggregated**
  `No feature` count across all bands.
- **BAM output**: band reads are still "no feature" for the tagged BAM — excluded
  by default and written (with the band name in `ZB`/`ZF`) only under
  `--include-no-feature`, exactly like `__no_feature__`.
- **Fallback**: a surviving no-feature read that can't be placed in a band
  (unmapped, or on a contig absent from the cytoBand file) keeps the generic
  `__no_feature__` label.

Band mode adds band columns to the read-annotation matrix and the Jaccard matrix
for the surviving no-feature reads, a modest increase in memory and runtime. The
`--skip-jaccard`, `--low-memory`, and `--min-reads` mitigations apply to band
symbols unchanged.

## Ambiguity Resolution

When a read group overlaps multiple annotations, telescope_rs resolves ambiguity in three phases:

1. **Representative selection.** For each annotation a read group overlaps, the best representative alignment is picked (tiebreakers: highest AS → proper pair → greatest overlap → first in input order).
2. **Priors → EM posteriors.** The prior model (§ [Prior Models](#prior-models)) produces per-fragment priors, which are fed into the EM loop (§ [EM Reassignment](#em-reassignment)).
3. **Confidence pick.** The annotation with the highest posterior γ is selected per fragment. It receives a `ZF` tag only if γ ≥ the confidence threshold (default `0.9`).

Read groups with no confident annotation are excluded from output (or tagged `__no_feature__` with `--include-no-feature`). A threshold of `1.0` corresponds to "unique" mode — only read groups mapping to exactly one annotation pass.

## Stranded Mode

`--stranded <MODE>` (on `annotate`) drops read↔annotation overlaps where the read's strand is incompatible with the annotation's GTF strand (column 7). When unset, all overlaps are retained — equivalent to unstranded.

The four modes:

| Mode | Library | Compatibility |
|------|---------|---------------|
| `RF` | Paired, dUTP/TruSeq stranded — read1 antisense to transcript | Alignment-unit strand must be **opposite** the annotation strand |
| `FR` | Paired, sense-stranded — read1 sense to transcript | Alignment-unit strand must be **same** as the annotation strand |
| `F`  | Single-end, sense-stranded — read sense to transcript | Read strand must be **same** as the annotation strand |
| `R`  | Single-end, antisense — read antisense to transcript | Read strand must be **opposite** the annotation strand |

`F`/`R` and `FR`/`RF` share the same comparator internally; the four-value enum lets `detect-strand` pick a recommendation that matches RNA-seq convention based on whether the BAM is paired.

**Alignment-unit strand** is derived from BAM flags once per alignment unit (per record in single-end, per mate pair in paired-end):

- Paired with read1 present (`0x40`): use read1's mapped strand (`-` if `0x10` is set, else `+`).
- Paired with only read2 (`0x80`): use the *reverse* of read2's mapped strand. This recovers the implied read1 strand and lets the rules above apply uniformly.
- Single-end / neither bit set: use the read's own mapped strand.

**Annotations with unknown strand** (`.` in GTF column 7) are always retained. If a read has no determinable strand (e.g., unmapped), no filtering is applied to its overlaps.

Filtering happens before the read-annotation matrix is built, so dropped overlaps do not contribute mass to EM. A read group whose only overlaps are all strand-incompatible falls through to the no-feature path (excluded from output unless `--include-no-feature` is set).

## Detecting strandedness

`telescope_rs detect-strand <BAM>` samples the first `--sample` (default 10,000) BAM records that overlap a strand-known GTF feature, tallies how often each (annotation, read-role) hit lines up with the annotation's strand, and writes a TSV plus a recommendation:

```
annotation     strand  pct_r1_same  pct_r2_same
GENE_A         +       98.20        1.74
GENE_B         -       97.50        2.41
...
# recommended_mode: FR
```

Columns:
- `annotation` — value of the configured GTF attribute (default `gene_id`).
- `strand` — that annotation's GTF strand (`+`, `-`, or `.`).
- `pct_r1_same` — percent of sampled read1 (or single-end) hits whose mapped strand matches the annotation's strand.
- `pct_r2_same` — same for read2 hits. Empty when no read2 hit was seen for the annotation (e.g., single-end input).

Rows are sorted by total supporting reads, descending; ties broken alphabetically. The trailing `# recommended_mode:` line gives one of `RF`, `FR`, `F`, `R`, or `unstranded`. Recommendation rules use a fixed 80% decisiveness threshold:

- **No read2 seen** (single-end): `pct_r1_same ≥ 80%` → `F`; `pct_r1_same ≤ 20%` → `R`; otherwise `unstranded`.
- **Paired**: `pct_r1_same ≥ 80%` and `pct_r2_same ≤ 20%` → `FR`; the mirror case → `RF`; otherwise `unstranded`.

The summary is also printed to stderr (suppress with `--quiet`).

## Output

Results are written to the output directory (default `tele_out/`).

From `annotate`:

- **`{stem}_annotated.bam`** -- Tagged BAM file with annotation results
- **`{stem}_jaccard.tsv`** -- Jaccard similarity between annotations (sparse tuple format: `annotation_a`, `annotation_b`, `jaccard`, `overlap`; filtered by `--similarity-threshold` and `--min-reads`). Sorted alphabetically by default; written unsorted under `--low-memory`
- **`{stem}_summary.tsv`** -- Per-annotation summary with final, initial, and unique counts plus prior/posterior fractional-assignment mass, sorted by final count descending
- **`{stem}_run_summary.txt`** -- Run statistics including record counts, confident/ambiguous/no-feature breakdowns, elapsed time, and peak RAM

From `detect-strand`:

- **`{stem}_strandedness.tsv`** -- Per-annotation strand-concordance summary with a `# recommended_mode:` footer

The input BAM filename stem is used to name the output files. With `--debug`, the annotated BAM is not written, but all other outputs are still produced. `--stranded` does not change the output schema — only which overlaps populate the matrix.

### BAM tags

| Tag | Description |
|-----|-------------|
| `ZB` | Initial annotation name (GTF attribute value) -- set on every output alignment |
| `ZF` | Confident annotation assignment -- set only when the EM posterior γ meets the confidence threshold |
| `ZR` | Representative flag (`1` = representative, `0` = not). Only present with `--all-alignments` |

In paired-end mode, both mates of a pair receive identical tags.

### Summary TSV columns

| Column | Description |
|--------|-------------|
| `annotation` | Annotation name |
| `final_count` | Number of read groups confidently assigned to this annotation (γ ≥ threshold) |
| `initial_count` | Number of read groups with any alignment overlapping this annotation (multimapping allowed) |
| `unique_count` | Number of read groups whose row contained only this annotation (i.e. would pass a confidence threshold of 1.0) |
| `length` | Total length (bp) of this annotation, computed as the union of all GTF intervals tagged with the annotation's attribute value. The synthetic `__no_feature__` row reports a fixed `1000000`; per-band `__no_feature_<chrom><band>__` rows (band mode) report the band's real width |
| `prior_mass` | Fractional assignment by the **prior**: each fragment contributes its unit of mass split across its support per the prior distribution, summed per annotation (unweighted). Each column totals the number of annotated fragments |
| `posterior_mass` | Same fractional assignment by the **finalized EM posterior** (γ). For the `ribbonfish` model this equals the converged read distribution θ; for `telescope` it is the expected fragment count Σγ (not π, which additionally weights by per-fragment AS) |

### Run summary

The run summary (`{stem}_run_summary.txt`) is also printed to stderr unless `--quiet` is set:

```
telescope_rs run summary
==========================
Total BAM records:       1234567
Read groups:              617283
Confident assignments:    500000
Ambiguous (dropped):      100000
No feature:                17283
Pair validation fails:        12
Confidence threshold:       0.90
Threads used:                  8
Elapsed time:             12.34s
Peak RAM:                  1.2 GB
```

"Pair validation fails" only appears when failures occur. "Peak RAM" only appears when the runtime can read RSS for the process. "Threads used" reflects the effective thread count actually used by rayon (which is also what `--low-memory` and `--threads N` control).

## Logging

By default, telescope_rs prints progress checkpoints to stderr with elapsed time at each pipeline stage:

```
[   0.45s] Loaded GTF: 12345 features across 25 chromosomes
[   1.23s] Read 1234567 BAM records
[   1.50s] Paired into 617283 alignment units (600000 mate pairs, 17283 unpaired)
[   3.67s] Annotated 617283 alignment units in parallel
[   4.01s] Grouped into 500000 read groups
[   4.12s] Built read-annotation matrix
[   4.20s] Computed per-fragment priors
[   4.80s] Ran EM for 10 iterations
[   4.90s] Resolved ambiguity (400000 confident, 83000 ambiguous, 17000 no-feature)
[   4.95s] Computed 142 Jaccard pairs (threshold 0.05)
[   5.10s] Built 500000 output directives
[   6.78s] Wrote output to tele_out/
[   6.78s] Peak RAM: 1.2 GB
[   6.78s] Done in 6.78s
```

Use `--quiet` to suppress all stderr output, or `--verbose` for additional detail (e.g., number of annotated units, unique annotation count).

## Pipeline

1. Load GTF annotation index (with string interning for attribute values). When `--band-no-feature` is on, also load the cytoBand file, interning each band as a synthetic `__no_feature_<chrom><band>__` annotation
2. Stream BAM records sequentially, extracting lightweight fields including mate pairing info; track global min/max AS across the file (used by the telescope prior)
3. Pair mates: match adjacent BAM records with complementary first/last segment flags, validate per SAM spec, merge spans, and sum AS scores. In `--single-end` mode, skip pairing
4. Annotate each alignment unit with all overlapping GTF features (parallel). For mate pairs, overlapping regions between the two mates are not double-counted
5. Group alignment units by QNAME
6. Build a sparse read-annotation matrix: for each annotation a read group overlaps, select a representative alignment unit using tiebreaking: highest AS score → proper pair → greatest overlap → first in input order
7. Compute per-fragment priors under the selected `--model` (parallel)
8. Run EM until `--max-iter` or max-mass-delta convergence (`--min-mass-delta`), whichever comes first (sequential outer loop, parallel E-step and M-reduce)
9. Apply the confidence threshold to the final posterior γ and pick one annotation per fragment (parallel)
10. Build a sparse co-occurrence matrix and compute pairwise Jaccard similarity and overlap coefficient between annotations that share at least one read (filtered by `--similarity-threshold` and `--min-reads`)
11. Stream through the input BAM a second time, applying tags and writing the output BAM (both mates of a pair receive the same tags). Skipped when `--debug` is set
12. Write Jaccard TSV, summary TSV, and run summary to output directory

## Default Behavior

With no optional flags, `telescope_rs annotate`:

- Treats reads as **unstranded** — annotation column 7 is parsed but not used to filter overlaps. Pass `--stranded` to enable strand-aware filtering, or run `detect-strand` first to pick a mode
- Pairs adjacent BAM records into mate pairs using SAM flags (0x40/0x80), validating flag correspondence and RNEXT/PNEXT
- Merges mate pair spans before computing annotation overlap (no double-counting)
- Sums alignment scores across mates in a pair
- Groups alignment units by QNAME
- For each annotation, selects the best representative alignment unit per read group
- Drops alignment ↔ annotation overlaps below 30 bp (combined per mate pair); see `--min-overlap`
- Computes per-fragment priors using the `ribbonfish` softmax with temperature 5
- Runs EM for up to 100 iterations (stopping early on convergence) to share abundance information across fragments, with **length correction enabled** (annotation length divides θ in the E-step). The default flips to **off** under `--model telescope`
- Tags the highest-posterior annotation with `ZF` when γ ≥ 0.9 (on both mates)
- When a read overlaps multiple exons of the same gene, overlap is summed across all exons sharing the same attribute value
- Partitions no-feature reads into cytogenetic bands from `reference/cytoBand.txt` (banding defaults **on** under `ribbonfish`); pass `--band-no-feature off` for the single `__no_feature__` bucket. See [Band Mode](#band-mode-no-feature-partitioning)
- Reads with no confident annotation are excluded from output
- Computes Jaccard similarity only between annotation pairs that co-occur in at least one read group, keeping pairs with overlap coefficient ≥ 0.05
- Prints progress checkpoints and a run summary to stderr
- Writes output files to `tele_out/`
- Uses all available CPU cores for parallel processing

## Memory

`telescope_rs`'s peak RAM scales with the number of read groups, the average annotations per group, and (most steeply) annotation co-occurrence density — the Jaccard `SparseMatrix` is the largest single allocator on dense datasets, easily exceeding the rest of the pipeline combined. On a real 28M-group / 164M-overlap dataset the default mode peaked at ~60 GB; a structurally similar but sparser tester only peaks at ~8 GB.

Two flags address this:

- **`--skip-jaccard`** turns off Jaccard computation entirely. No `SparseMatrix`, no `{stem}_jaccard.tsv`. Use this when you don't need the similarity matrix.
- **`--low-memory`** still computes Jaccard but **streams** pairs straight to the TSV — no in-memory tuples buffer, no sort. The `_jaccard.tsv` is therefore unsorted; `sort` it post-hoc if you need it ordered. `--low-memory` also defaults `--threads` to 1 (overridable) and rejects `--all-alignments`, since the directive list scales by record count rather than group count there. Expect ~2× slowdown. Combine with `--skip-jaccard` to also skip computing Jaccard.
- **`--min-reads <N>`** (default 10) drops Jaccard pairs where either annotation has fewer than N supporting read groups. Independent of `--low-memory`; reduces noise from low-coverage annotations.

Several representation changes are always-on (no flag), so they help every run:
- The read-annotation matrix uses sorted `Vec<(Symbol, MatrixCell)>` rows instead of `HashMap`, removing per-row hash overhead at no cost for the typical row size.
- `OutputDirective` carries interned `Symbol`s instead of heap `String`s; annotation names are resolved to bytes only at BAM-write time.
- The EM step consumes the priors `Vec` and returns owned `FragmentPosteriors`, so the priors are freed before the resolution and output phases run.
