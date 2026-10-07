//! Second pass of `--updated_sam` (Telescope's `Telescope.update_sam`).
//!
//! Reads back the tagged BAM the loader wrote and annotates each alignment
//! with the model's verdict:
//!   * `ZT:SEC` alignments (not the best for their feature): flagged
//!     secondary, MAPQ 0, grey `YC` colour;
//!   * `ZT:PRI` alignments: MAPQ = phred of the membership weight, `XP` =
//!     weight as a percentage, primary and vermilion if the fragment was
//!     assigned to that feature, otherwise secondary (yellow at >= 20 %,
//!     pale green below).
//!
//! `Content::All` writes every record of every overlapping fragment, as
//! Telescope does. `Content::Assigned` writes only the alignment a fragment
//! was actually assigned to (and nothing for unassigned fragments).

use std::collections::HashMap;
use std::io;

use super::loader::{Aln, Loaded, MateKey, Rec, TAG_FEATURE, TAG_RANK, classify, name_hash};
use super::rawbam::{CompressionLevel, Raw, RawReader, RawWriter, set_flag, set_mapq, set_tag_str, set_tag_u8};
use super::model::ModelFit;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Content {
    All,
    Assigned,
}

const TAG_PERCENT: [u8; 2] = *b"XP";
const TAG_COLOUR: [u8; 2] = *b"YC";
const SECONDARY: u16 = 0x100;
// colors.py: (248,248,248), D2PAL['vermilion'], D2PAL['yellow'], GPAL[2]
const GREY: &str = "248,248,248";
const VERMILION: &str = "217,95,2";
const YELLOW: &str = "230,171,2";
const PALE_GREEN: &str = "209,236,228";

/// `helpers.phred`
fn phred(p: f64) -> u8 {
    if p < 1.0 { (-10.0 * (1.0 - p).log10()).round_ties_even() as u8 } else { 255 }
}

fn invalid(msg: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg.to_string())
}

/// `assigned` is the reassignment matrix for the chosen mode
/// ([`super::report::reassign_entries`]).
pub fn write(
    tagged_path: &str,
    out_path: &str,
    loaded: &Loaded,
    fit: &ModelFit,
    assigned: &[f64],
    content: Content,
) -> io::Result<()> {
    let m = &loaded.matrix;
    let rows = loaded.rows.as_ref().expect("row index kept for updated BAM");
    let cols: HashMap<&[u8], u32> =
        loaded.feat_names.iter().enumerate().map(|(j, n)| (n.as_bytes(), j as u32)).collect();

    let mut reader = RawReader::open(tagged_path)?;
    let mut out = RawWriter::create(out_path, &loaded.header, CompressionLevel::default())?;

    let mut recs: Vec<Rec> = Vec::new();
    let mut raw: Vec<Vec<u8>> = Vec::new();
    let mut n = 0usize;
    let mut cur_name: Vec<u8> = Vec::new();
    let mut alns: Vec<Aln> = Vec::new();
    let mut cache: Vec<(MateKey, usize, bool)> = Vec::new();
    let mut record: Vec<u8> = Vec::new();

    let mut flush = |name: &[u8], recs: &[Rec], raw: &mut [Vec<u8>]| -> io::Result<()> {
        // Telescope re-reads its temporary BAM with the same fragment logic.
        classify(recs, &mut alns, &mut cache);
        let row = *rows.get(&name_hash(name)).ok_or_else(|| invalid("tagged read missing from matrix"))? as usize;
        for &aln in &alns {
            let members = || std::iter::once(aln.r1).chain(aln.r2);
            if recs[aln.r1].is_unmapped() {
                if content == Content::All {
                    for r in members() {
                        out.write(&raw[r])?;
                    }
                }
                continue;
            }
            let rank = Raw(&raw[aln.r1]).aux_str(TAG_RANK).ok_or_else(|| invalid("Missing ZT tag"))?;
            let (mapq, percent, secondary, colour) = if rank == b"SEC" {
                if content == Content::Assigned {
                    continue;
                }
                (0, None, true, GREY)
            } else {
                let feat = Raw(&raw[aln.r1]).aux_str(TAG_FEATURE).ok_or_else(|| invalid("Missing ZF tag"))?;
                let col = *cols.get(feat).ok_or_else(|| invalid("unknown feature in ZF tag"))?;
                let entry = m.row(row).find(|&k| m.indices[k] == col);
                let prob = entry.map_or(0.0, |k| fit.z[k]);
                let is_assigned = entry.is_some_and(|k| assigned[k] > 0.0);
                if content == Content::Assigned && !(is_assigned && col != 0) {
                    continue;
                }
                let percent = (prob * 100.0).round_ties_even() as u8;
                let colour = if is_assigned {
                    VERMILION
                } else if prob >= 0.2 {
                    YELLOW
                } else {
                    PALE_GREEN
                };
                (phred(prob), Some(percent), !is_assigned, colour)
            };
            for r in members() {
                let rec = &mut raw[r];
                let bits = Raw(rec).flag();
                set_flag(rec, if secondary { bits | SECONDARY } else { bits & !SECONDARY });
                set_mapq(rec, mapq);
                if let Some(p) = percent {
                    set_tag_u8(rec, TAG_PERCENT, p);
                }
                set_tag_str(rec, TAG_COLOUR, colour.as_bytes());
                out.write(rec)?;
            }
        }
        Ok(())
    };

    while reader.read(&mut record)? {
        let name = Raw(&record).name();
        if n > 0 && name != cur_name.as_slice() {
            flush(&cur_name, &recs[..n], &mut raw[..n])?;
            n = 0;
        }
        if n == 0 {
            cur_name.clear();
            cur_name.extend_from_slice(name);
        }
        if n == recs.len() {
            recs.push(Rec::default());
            raw.push(Vec::new());
        }
        recs[n].fill(Raw(&record));
        raw[n].clone_from(&record);
        n += 1;
    }
    if n > 0 {
        flush(&cur_name, &recs[..n], &mut raw[..n])?;
    }
    out.finish()
}

#[cfg(test)]
mod tests {
    use super::phred;

    #[test]
    fn phred_matches_helper() {
        assert_eq!(phred(0.9), 10);
        assert_eq!(phred(0.999999), 60);
        assert_eq!(phred(0.0), 0);
        assert_eq!(phred(1.0), 255);
    }
}
