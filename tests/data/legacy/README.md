# Reference fixtures for `assign`

Small real-data inputs, with the outputs **Python Telescope v1.0.4.1 from the
[hanalysis fork](https://github.com/hanalysis/telescope)** (commit `bb58c43`, the
`telescope/` directory of this repository) produced for them. `tests/legacy_fixtures.rs`
runs every line of `cases.tsv` and requires a byte-for-byte match with `expected/`.

Reference environments: Python 3.10.21 with numpy 1.26.4, scipy 1.15.2, pysam 0.24.1 and
intervaltree 3.1.0 (`py310`), and Python 3.7.12 with numpy 1.21.6 and scipy 1.7.3 (`py37`).
Cases marked `rust` have no Python equivalent; their expected output was written by this
program and only guards against unintended change.

## What each fixture is for

| Fixture | Source | Fragments | What it exercises |
|---|---|---|---|
| `pe_hs1_herv.bam` + `hs1_herv_subset.gtf` | CCLE SRR8618300, HISAT2 `-k 100`, T2T-CHM13, paired, unstranded | 15,848 | Every fragment of that run's 10% subsample that overlaps the HERV annotation, plus a sample of the rest. 14,179 matrix rows (numpy sums restart every 8,192 values), rows of up to 67 candidates (pairwise summation), 840 candidates dropped at zero weight, 1,828 tied alignments, EM that hits the iteration cap. Reverting to plain left-to-right sums changes this report. Also carries the `--theta_prior 0`, `--use_likelihood`, reassign-mode, Python 3.7 and `no-theta` cases |
| `pe_hg38_hml2.bam` + `reference/hg38.hml2_only.gtf` | same run, hg38, chromosomes named `1`, `2`, ... | 1,901 | A GTF whose attribute quoting is mangled (`locus """"HML2.LTR1"""";`): locus names keep literal quotes, as they do in Telescope. Rows of up to 90 candidates |
| `pe_hisat2.bam` | same run, hg38 | 4,000 | Paired-end HISAT2: spliced alignments, non-positive scores, mixed pairs |
| `pe_bowtie2.bam` | same reads, bowtie2 `-k 100 --very-sensitive-local` | 2,500 | Local alignments with soft clips and positive scores; many secondary alignments |
| `pe_star.bam` | same reads, STAR 2.7.10b, up to 100 alignments | 2,500 | The pair's score on both mates; STAR's flags |
| `se_hisat2.bam`, `se_bowtie2.bam` | read 1 only of the same run | 2,500 each | Single-end input |
| `pe_stranded_hisat2.bam` | ENCODE ENCSR000CON (A549, reverse-stranded), HISAT2 `--rna-strandness RF` | 2,500 | `--stranded_mode RF` and `FR`, and the same file run unstranded |
| `se_stranded_hisat2.bam` | read 1 only of the same library, `--rna-strandness R` | 2,500 | `--stranded_mode R` and `F` |
| `long_minimap2.sam` | ENCODE ENCFF694INI (K562, PacBio Sequel), minimap2 `splice:hq` with pbmm2's Iso-Seq scoring, `--eqx -Y -N 100` | 1,168 | Long reads with 492 supplementary (chimeric) and 7,598 secondary records, `=`/`X` CIGAR operations, and SAM rather than BAM input |

The last eight share `hg38_window_herv_genes.gtf`: every row of the HERV+genes annotation
in hg38 chr19:12,000,000-13,000,000, a zinc-finger gene cluster where loci overlap. Each of
those fixtures has between 2 and 241 alignments that overlap two loci equally, so the tie
order matters in all of them. Switching to corrected ties or coordinates, or to the other
Python's hash, changes the report (checked on `pe_hisat2`, and for ties on `long_minimap2`).

Five cases also compare the updated BAM (`*.updated.tsv`: name, flag, reference, position,
MAPQ and the `ZF`, `ZT`, `ZB`, `XP`, `YC` tags of every record).

## How the files were made

- `build_fixtures.py` cut the alignment files and annotations out of full alignments. Whole
  fragments are kept in their original order. SEQ and QUAL are replaced by `*`, since
  Telescope never reads them. GTF rows keep only `gene_id`, `transcript_id` and `locus`.
- `make_expected.sh` ran the reference on each case.

Both need data that is not in the repository, so they are a record rather than part of the
test run. Add a case by appending a line to `cases.tsv` and rerunning `make_expected.sh`.

## What is not covered here

Coordinate-sorted input, other platforms' math libraries, and inputs large enough to stress
memory. The annotations are subsets, so tie outcomes here are those of the subset, not of
the full annotation.
