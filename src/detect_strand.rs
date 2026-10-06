//! Strandedness detection.
//!
//! Streams a BAM file and, for the first `sample` records that overlap a
//! GTF feature with a known strand, tallies how often each (annotation,
//! read-role) pair lines up with the annotation's strand. Writes a per-
//! annotation TSV plus a recommendation comment.
//!
//! Read-role:
//! - read1: BAM record has `0x40` (first segment).
//! - read2: BAM record has `0x80` (last segment) and not `0x40`.
//! - single: neither bit set; counted in the read1 column.

use std::collections::HashMap;
use std::fs;
use std::io::{self, Write};
use std::path::Path;

use noodles::{bam, sam::alignment::RecordBuf};

use crate::cli::DetectStrandArgs;
use crate::gtf::FeatureIndex;
use crate::intern::{Interner, Symbol};
use crate::logging::{Logger, Verbosity};
use crate::overlap::reference_span;

/// Threshold above which `pct_r1_same` (or the analogous read2 metric) is
/// considered decisive when picking a recommended mode.
const DECISIVE_PCT: f64 = 80.0;

#[derive(Default, Clone, Copy, Debug)]
struct Tally {
    same: u64,
    total: u64,
}

impl Tally {
    fn add(&mut self, same: bool) {
        if same {
            self.same += 1;
        }
        self.total += 1;
    }

    fn pct(&self) -> Option<f64> {
        if self.total == 0 {
            None
        } else {
            Some(100.0 * self.same as f64 / self.total as f64)
        }
    }
}

#[derive(Default, Clone, Copy, Debug)]
struct AnnotationTally {
    strand: Option<char>,
    r1: Tally,
    r2: Tally,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Recommendation {
    Fr,
    Rf,
    F,
    R,
    Unstranded,
}

impl Recommendation {
    fn as_str(self) -> &'static str {
        match self {
            Recommendation::Fr => "FR",
            Recommendation::Rf => "RF",
            Recommendation::F => "F",
            Recommendation::R => "R",
            Recommendation::Unstranded => "unstranded",
        }
    }
}

fn recommend(global_r1: Tally, global_r2: Tally) -> Recommendation {
    let r1 = global_r1.pct().unwrap_or(50.0);
    let r2_seen = global_r2.total > 0;
    let r2 = global_r2.pct().unwrap_or(50.0);

    if !r2_seen {
        if r1 >= DECISIVE_PCT {
            return Recommendation::F;
        }
        if r1 <= 100.0 - DECISIVE_PCT {
            return Recommendation::R;
        }
        return Recommendation::Unstranded;
    }

    if r1 >= DECISIVE_PCT && r2 <= 100.0 - DECISIVE_PCT {
        return Recommendation::Fr;
    }
    if r1 <= 100.0 - DECISIVE_PCT && r2 >= DECISIVE_PCT {
        return Recommendation::Rf;
    }
    Recommendation::Unstranded
}

pub fn run_detect_strand(args: DetectStrandArgs) -> io::Result<()> {
    let verbosity = if args.quiet {
        Verbosity::Quiet
    } else if args.verbose {
        Verbosity::Verbose
    } else {
        Verbosity::Normal
    };
    let log = Logger::new(verbosity);

    let index = FeatureIndex::from_file(&args.gtf, &args.field)?;
    log.stage(&format!(
        "Loaded GTF: {} features across {} chromosomes",
        index.num_features(),
        index.num_chromosomes(),
    ));

    let mut reader = bam::io::reader::Builder.build_from_path(&args.bam)?;
    let header = reader.read_header()?;

    let mut tallies: HashMap<Symbol, AnnotationTally> = HashMap::new();
    let mut global_r1 = Tally::default();
    let mut global_r2 = Tally::default();
    let mut sampled: usize = 0;
    let mut total_seen: usize = 0;
    let mut record_buf = RecordBuf::default();

    'outer: while sampled < args.sample {
        match reader.read_record_buf(&header, &mut record_buf) {
            Ok(0) => break 'outer,
            Ok(_) => {
                total_seen += 1;
                let Some((chrom, start, stop)) = reference_span(&record_buf, &header) else {
                    continue;
                };
                let flags = record_buf.flags();
                let is_first = flags.is_first_segment();
                let is_last = flags.is_last_segment();
                let is_reverse = flags.is_reverse_complemented();
                let read_strand = if is_reverse { '-' } else { '+' };
                // Read1 if first-segment is set, OR if neither first nor last
                // is set (true single-end, no role bit).
                let read_role = if is_last && !is_first {
                    ReadRole::ReadTwo
                } else {
                    ReadRole::ReadOne
                };

                let hits = index.find(chrom, start, stop, args.normalize_chr);
                let mut contributed = false;
                for iv in &hits {
                    let Some(feat_strand) = iv.val.strand else {
                        continue;
                    };
                    let Some(&sym) = iv.val.attributes.get(&args.field) else {
                        continue;
                    };
                    let same = read_strand == feat_strand;
                    let entry = tallies.entry(sym).or_insert(AnnotationTally {
                        strand: Some(feat_strand),
                        ..AnnotationTally::default()
                    });
                    match read_role {
                        ReadRole::ReadOne => {
                            entry.r1.add(same);
                            global_r1.add(same);
                        }
                        ReadRole::ReadTwo => {
                            entry.r2.add(same);
                            global_r2.add(same);
                        }
                    }
                    contributed = true;
                }
                if contributed {
                    sampled += 1;
                }
            }
            Err(e) => return Err(e),
        }
    }

    drop(reader);
    drop(record_buf);

    log.stage(&format!(
        "Sampled {} annotation-overlapping records (scanned {})",
        sampled, total_seen,
    ));

    let recommendation = recommend(global_r1, global_r2);

    fs::create_dir_all(&args.output_dir)?;
    let stem = Path::new(&args.bam)
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("output");
    let out_path = Path::new(&args.output_dir).join(format!("{stem}_strandedness.tsv"));
    write_tsv(&out_path, &tallies, &index.interner, recommendation)?;
    log.stage(&format!("Wrote {}", out_path.display()));

    if verbosity != Verbosity::Quiet {
        eprintln!(
            "rusty_telescope detect-strand summary\n\
             =====================================\n\
             Sampled records:        {}\n\
             read1 same-strand pct:  {:.2}\n\
             read2 same-strand pct:  {}\n\
             Recommended mode:       {}",
            sampled,
            global_r1.pct().unwrap_or(0.0),
            global_r2
                .pct()
                .map(|v| format!("{v:.2}"))
                .unwrap_or_else(|| "n/a (single-end?)".to_string()),
            recommendation.as_str(),
        );
    }

    Ok(())
}

#[derive(Clone, Copy)]
enum ReadRole {
    ReadOne,
    ReadTwo,
}

fn write_tsv(
    path: &Path,
    tallies: &HashMap<Symbol, AnnotationTally>,
    interner: &Interner,
    rec: Recommendation,
) -> io::Result<()> {
    let mut rows: Vec<(&str, AnnotationTally)> = tallies
        .iter()
        .map(|(sym, t)| (interner.resolve(*sym), *t))
        .collect();
    // Sort by total reads supporting the annotation, descending; tie-break on
    // annotation name for deterministic output.
    rows.sort_by(|a, b| {
        let totals_a = a.1.r1.total + a.1.r2.total;
        let totals_b = b.1.r1.total + b.1.r2.total;
        totals_b.cmp(&totals_a).then_with(|| a.0.cmp(b.0))
    });

    let mut f = fs::File::create(path)?;
    writeln!(f, "annotation\tstrand\tpct_r1_same\tpct_r2_same")?;
    for (name, t) in rows {
        let strand_col = match t.strand {
            Some(c) => c.to_string(),
            None => ".".to_string(),
        };
        let r1_col = t.r1.pct().map(|v| format!("{v:.2}")).unwrap_or_default();
        let r2_col = t.r2.pct().map(|v| format!("{v:.2}")).unwrap_or_default();
        writeln!(f, "{name}\t{strand_col}\t{r1_col}\t{r2_col}")?;
    }
    writeln!(f, "# recommended_mode: {}", rec.as_str())?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn t(same: u64, total: u64) -> Tally {
        Tally { same, total }
    }

    #[test]
    fn test_recommend_se_f() {
        // SE: r2 not seen, r1 same dominates.
        assert_eq!(recommend(t(95, 100), Tally::default()), Recommendation::F);
    }

    #[test]
    fn test_recommend_se_r() {
        assert_eq!(recommend(t(5, 100), Tally::default()), Recommendation::R);
    }

    #[test]
    fn test_recommend_se_unstranded() {
        assert_eq!(
            recommend(t(50, 100), Tally::default()),
            Recommendation::Unstranded
        );
    }

    #[test]
    fn test_recommend_fr() {
        assert_eq!(recommend(t(95, 100), t(5, 100)), Recommendation::Fr);
    }

    #[test]
    fn test_recommend_rf() {
        assert_eq!(recommend(t(5, 100), t(95, 100)), Recommendation::Rf);
    }

    #[test]
    fn test_recommend_pe_unstranded_when_both_concordant() {
        // Both reads same-strand to the annotation → not RNA strand-specific.
        assert_eq!(
            recommend(t(95, 100), t(95, 100)),
            Recommendation::Unstranded
        );
    }

    #[test]
    fn test_recommend_threshold_boundary() {
        // 80/20 boundary: should still be FR.
        assert_eq!(recommend(t(80, 100), t(20, 100)), Recommendation::Fr);
    }
}
