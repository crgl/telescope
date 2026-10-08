//! Per-base coverage of assigned fragments, written as bigWig (`--bigwig`).
//!
//! Each assigned fragment contributes the reference bases its assigned
//! alignment covers: aligned blocks only (introns and deletions excluded),
//! with the two mates of a pair merged so their overlap counts once. Depth
//! is the number of such fragments over each base, unnormalised.

use std::collections::HashMap;
use std::io;

use bigtools::beddata::BedParserStreamingIterator;
use bigtools::{BigWigWrite, Value};
use noodles::sam::Header;

/// Collects covered blocks per strand track and reference sequence.
pub struct Coverage {
    /// [track][reference id] -> blocks (0-based, half-open)
    blocks: Vec<Vec<Vec<(u32, u32)>>>,
}

impl Coverage {
    /// `tracks` is 1 for a single file, 2 to split plus (0) and minus (1).
    pub fn new(tracks: usize, n_refs: usize) -> Self {
        Coverage { blocks: vec![vec![Vec::new(); n_refs]; tracks] }
    }

    pub fn tracks(&self) -> usize {
        self.blocks.len()
    }

    /// Adds one fragment. `blocks` are sorted by start; overlapping or
    /// touching ones are merged here so no base is counted twice.
    pub fn add(&mut self, track: usize, ref_id: usize, blocks: &[(i64, i64)]) {
        let Some(dst) = self.blocks[track].get_mut(ref_id) else { return };
        let mut cur: Option<(i64, i64)> = None;
        for &(s, e) in blocks {
            match cur {
                Some((cs, ce)) if s <= ce => cur = Some((cs, ce.max(e))),
                Some((cs, ce)) => {
                    dst.push((cs.max(0) as u32, ce.max(0) as u32));
                    cur = Some((s, e));
                }
                None => cur = Some((s, e)),
            }
        }
        if let Some((cs, ce)) = cur {
            dst.push((cs.max(0) as u32, ce.max(0) as u32));
        }
    }

    /// Writes one track; returns false (and writes nothing) if it is empty.
    pub fn write(&mut self, track: usize, path: &str, header: &Header) -> io::Result<bool> {
        let refs: Vec<(String, u32)> = header
            .reference_sequences()
            .iter()
            .map(|(name, seq)| (String::from_utf8_lossy(name.as_ref()).into_owned(), usize::from(seq.length()) as u32))
            .collect();
        let per_ref = std::mem::take(&mut self.blocks[track]);
        if per_ref.iter().all(Vec::is_empty) {
            return Ok(false);
        }
        let chrom_sizes: HashMap<String, u32> = refs.iter().cloned().collect();
        let intervals = per_ref
            .into_iter()
            .enumerate()
            .flat_map(|(i, blocks)| depth_intervals(blocks).into_iter().map(move |v| (i, v)))
            .map(|(i, (start, end, depth))| (refs[i].0.clone(), Value { start, end, value: depth as f32 }));
        let runtime = tokio::runtime::Builder::new_current_thread().build()?;
        BigWigWrite::create_file(path, chrom_sizes)?
            .write(BedParserStreamingIterator::wrap_infallible_iter(intervals, true), runtime)
            .map_err(|e| io::Error::other(e.to_string()))?;
        Ok(true)
    }
}

/// Runs of constant, non-zero depth from a set of blocks: `(start, end, depth)`.
fn depth_intervals(blocks: Vec<(u32, u32)>) -> Vec<(u32, u32, u32)> {
    let mut starts: Vec<u32> = blocks.iter().map(|b| b.0).collect();
    let mut ends: Vec<u32> = blocks.iter().map(|b| b.1).collect();
    drop(blocks);
    starts.sort_unstable();
    ends.sort_unstable();
    let mut out: Vec<(u32, u32, u32)> = Vec::new();
    let (mut i, mut j, mut depth, mut prev) = (0, 0, 0u32, 0u32);
    while j < ends.len() {
        // next position where the depth changes
        let pos = if i < starts.len() && starts[i] < ends[j] { starts[i] } else { ends[j] };
        if depth > 0 && pos > prev {
            match out.last_mut() {
                Some(last) if last.1 == prev && last.2 == depth => last.1 = pos,
                _ => out.push((prev, pos, depth)),
            }
        }
        while i < starts.len() && starts[i] == pos {
            depth += 1;
            i += 1;
        }
        while j < ends.len() && ends[j] == pos {
            depth -= 1;
            j += 1;
        }
        prev = pos;
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn depth_from_overlapping_blocks() {
        // 10-20, 15-30, 30-40 and a second copy of 15-30
        let got = depth_intervals(vec![(15, 30), (10, 20), (30, 40), (15, 30)]);
        assert_eq!(got, vec![(10, 15, 1), (15, 20, 3), (20, 30, 2), (30, 40, 1)]);
    }

    #[test]
    fn mates_are_merged_before_counting() {
        let mut c = Coverage::new(1, 1);
        // mate blocks 100-150 and 130-180 (overlap), then a spliced block 300-320
        c.add(0, 0, &[(100, 150), (130, 180), (300, 320)]);
        assert_eq!(c.blocks[0][0], vec![(100, 180), (300, 320)]);
        assert_eq!(depth_intervals(c.blocks[0][0].clone()), vec![(100, 180, 1), (300, 320, 1)]);
    }
}
