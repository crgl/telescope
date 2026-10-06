//! GTF loading and fragment/feature overlap, as Telescope's
//! `_AnnotationIntervalTree` does it.
//!
//! `OverlapCompat::Telescope` reproduces the reference quirks:
//!   * every GTF row is used, whatever its feature type;
//!   * a row spanning 1-based `[start, end]` is stored as `[start, end + 1)`
//!     and a 0-based alignment block `[s, e)` is queried as `[s, e + 1)`, so
//!     overlaps are measured one base to the right and can exceed the true
//!     overlap by one;
//!   * a row that bridges two existing intervals of its own locus aborts the
//!     run (Telescope's `assert len(mergeable) == 1`);
//!   * when two loci overlap an alignment equally, the winner is whichever
//!     Python's set iteration yields first. That needs the ported interval
//!     tree in [`super::pytree`], so this mode builds and queries that tree.
//!
//! `OverlapCompat::Corrected` keeps the same model but with true coordinates,
//! no abort, and ties going to the locus that appears first in the GTF. Both are selected at run time, so the quirks stay optional.

use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::{self, BufRead, BufReader};

use regex::Regex;
use rust_lapper::{Interval, Lapper};

use super::pyset::PySet;
use super::pytree::{Iv, IvTable, PyIntervalTree};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OverlapCompat {
    Telescope,
    Corrected,
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub struct Feature {
    pub locus: u32,
    pub strand: u8,
}

pub struct Annotation {
    /// Locus names in first-seen GTF order; the index is the locus id.
    pub loci: Vec<String>,
    /// Summed interval length per locus (Telescope's `feature_length`).
    pub lengths: Vec<u64>,
    index: Index,
    stranded: bool,
}

enum Index {
    Telescope { trees: HashMap<String, PyIntervalTree>, ivs: IvTable },
    Corrected(HashMap<String, Lapper<u64, Feature>>),
}

/// Reusable buffers for [`Annotation::best_feature`].
#[derive(Default)]
pub struct Scratch {
    totals: Vec<(u32, u64)>,
    result: PySet,
    points: PySet,
}

/// Best feature for one alignment, plus whether the choice hinged on a tie.
pub struct Hit {
    pub locus: u32,
    pub overlap: u64,
    pub tied: bool,
}

/// One locus's intervals on one chromosome, kept disjoint and keyed by
/// begin: (end, strand, id in the ported tree).
type LocusIntervals = BTreeMap<u64, (u64, u8, u32)>;

fn invalid(msg: String) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg)
}

impl Annotation {
    pub fn from_gtf(
        path: &str,
        attribute: &str,
        stranded: bool,
        compat: OverlapCompat,
    ) -> io::Result<Self> {
        // Same pattern Telescope hands to re.findall.
        let attr_re = Regex::new(r#"(\w+)\s+"(.+?)";"#).expect("static regex");
        let mut locus_ids: HashMap<String, u32> = HashMap::new();
        let mut loci: Vec<String> = Vec::new();
        let mut chrom_ids: HashMap<String, u32> = HashMap::new();
        let mut chrom_names: Vec<String> = Vec::new();
        let mut merged: HashMap<(u32, u32), LocusIntervals> = HashMap::new();
        // Telescope mode mirrors every add/remove onto the ported tree.
        let mut trees: Vec<PyIntervalTree> = Vec::new();
        let mut iv_table = IvTable::default();

        let reader = BufReader::new(File::open(path)?);
        for (lineno, line) in reader.lines().enumerate() {
            let line = line?;
            if line.starts_with('#') {
                continue;
            }
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() != 9 {
                return Err(invalid(format!(
                    "{path}:{}: expected 9 tab-separated fields, found {}",
                    lineno + 1,
                    f.len()
                )));
            }
            let mut key_val: Option<&str> = None;
            for cap in attr_re.captures_iter(f[8]) {
                if &cap[1] == attribute {
                    key_val = Some(cap.get(2).unwrap().as_str());
                }
            }
            let Some(key_val) = key_val else {
                return Err(invalid(format!(
                    "{path}:{}: missing attribute \"{attribute}\"",
                    lineno + 1
                )));
            };
            let parse = |s: &str| -> io::Result<u64> {
                s.trim()
                    .parse::<u64>()
                    .map_err(|_| invalid(format!("{path}:{}: bad coordinate '{s}'", lineno + 1)))
            };
            let (start, end) = (parse(f[3])?, parse(f[4])?);
            let (mut begin, mut stop) = match compat {
                OverlapCompat::Telescope => (start, end + 1),
                OverlapCompat::Corrected => (start.saturating_sub(1), end),
            };
            let strand = f[6].as_bytes().first().copied().unwrap_or(b'.');

            let locus = *locus_ids.entry(key_val.to_string()).or_insert_with(|| {
                loci.push(key_val.to_string());
                (loci.len() - 1) as u32
            });
            let chrom = *chrom_ids.entry(f[0].to_string()).or_insert_with(|| {
                chrom_names.push(f[0].to_string());
                trees.push(PyIntervalTree::default());
                (chrom_names.len() - 1) as u32
            });

            let ivs = merged.entry((chrom, locus)).or_default();
            // Intervals within a locus are disjoint, so those overlapping
            // [begin, stop) are the last few that start before `stop`.
            let hits: Vec<u64> = ivs
                .range(..stop)
                .rev()
                .take_while(|(_, (e, _, _))| *e > begin)
                .map(|(b, _)| *b)
                .collect();
            if hits.len() > 1 && compat == OverlapCompat::Telescope {
                return Err(invalid(format!(
                    "{path}:{}: row overlaps {} existing intervals of locus {key_val}; \
                     Telescope aborts here (assert len(mergeable) == 1)",
                    lineno + 1,
                    hits.len()
                )));
            }
            for b in hits {
                let (e, _, old) = ivs.remove(&b).unwrap();
                begin = begin.min(b);
                stop = stop.max(e);
                if compat == OverlapCompat::Telescope {
                    trees[chrom as usize].remove_interval(old, &iv_table);
                    iv_table.release(old);
                }
            }
            let mut id = 0;
            if compat == OverlapCompat::Telescope {
                id = iv_table.push(Iv::new(begin, stop, locus, strand));
                trees[chrom as usize].add(id, &iv_table);
            }
            ivs.insert(begin, (stop, strand, id));
        }

        let mut lengths = vec![0u64; loci.len()];
        let mut per_chrom: Vec<Vec<Interval<u64, Feature>>> = vec![Vec::new(); chrom_names.len()];
        for ((chrom, locus), ivs) in merged {
            for (begin, (stop, strand, _)) in ivs {
                lengths[locus as usize] += stop - begin;
                per_chrom[chrom as usize].push(Interval {
                    start: begin,
                    stop,
                    val: Feature { locus, strand },
                });
            }
        }
        iv_table.shrink();
        trees.iter_mut().for_each(PyIntervalTree::shrink);
        let index = match compat {
            OverlapCompat::Telescope => Index::Telescope {
                trees: chrom_names.into_iter().zip(trees).collect(),
                ivs: iv_table,
            },
            OverlapCompat::Corrected => Index::Corrected(
                chrom_names
                    .into_iter()
                    .zip(per_chrom)
                    .map(|(name, ivs)| (name, Lapper::new(ivs)))
                    .collect(),
            ),
        };
        Ok(Annotation { loci, lengths, index, stranded })
    }

    /// Telescope's `intersect_blocks` + `most_common()[0]`: total overlap per
    /// locus over the fragment's merged blocks, then the largest. Ties go to
    /// the locus encountered first.
    pub fn best_feature(
        &self,
        chrom: &str,
        blocks: &[(i64, i64)],
        frag_strand: u8,
        scratch: &mut Scratch,
    ) -> Option<Hit> {
        let Scratch { totals, result, points } = scratch;
        totals.clear();
        let stranded = self.stranded;
        let tally = |totals: &mut Vec<(u32, u64)>, locus: u32, strand: u8, ov: u64| {
            if stranded && strand != frag_strand {
                return;
            }
            match totals.iter_mut().find(|(l, _)| *l == locus) {
                Some((_, tot)) => *tot += ov,
                None => totals.push((locus, ov)),
            }
        };
        match &self.index {
            Index::Telescope { trees, ivs } => {
                let tree = trees.get(chrom)?;
                for &(bs, be) in blocks {
                    let (qs, qe) = (bs.max(0) as u64, (be + 1).max(0) as u64);
                    tree.overlap(qs, qe, result, points, ivs);
                    for id in result.iter() {
                        let iv = ivs.get(id);
                        tally(totals, iv.locus, iv.strand, iv.end.min(qe) - iv.begin.max(qs));
                    }
                }
            }
            Index::Corrected(chroms) => {
                let tree = chroms.get(chrom)?;
                for &(bs, be) in blocks {
                    let (qs, qe) = (bs.max(0) as u64, be.max(0) as u64);
                    if qe <= qs {
                        continue;
                    }
                    let first = totals.len();
                    for iv in tree.find(qs, qe) {
                        tally(totals, iv.val.locus, iv.val.strand, iv.stop.min(qe) - iv.start.max(qs));
                    }
                    // deterministic: loci new to this block enter in GTF order
                    totals[first..].sort_by_key(|t| t.0);
                }
            }
        }
        let best = totals.iter().map(|&(_, ov)| ov).max()?;
        let mut tied = totals.iter().filter(|&&(_, ov)| ov == best);
        let locus = tied.next().unwrap().0;
        Some(Hit { locus, overlap: best, tied: tied.next().is_some() })
    }
}
