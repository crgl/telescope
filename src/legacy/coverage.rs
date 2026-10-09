//! Per-base coverage of assigned fragments, written as bigWig (`--bigwig`).
//!
//! Each assigned alignment contributes the reference bases it covers:
//! aligned blocks only (introns and deletions excluded), with the two mates
//! of a pair merged so their overlap counts once. An alignment weighs what
//! the reassignment gave it: 1 when the whole fragment was assigned to it,
//! a fraction when the mode splits a fragment across candidates. Depth is
//! the sum of those weights over each base, unnormalised.

use std::collections::HashMap;
use std::io;

use bigtools::beddata::BedParserStreamingIterator;
use bigtools::{BigWigWrite, Value};
use noodles::sam::Header;

/// Collects covered blocks per strand track and reference sequence.
pub struct Coverage {
    /// [track][reference id] -> blocks (0-based, half-open) with their weight
    blocks: Vec<Vec<Vec<(u32, u32, f32)>>>,
}

impl Coverage {
    /// `tracks` is 1 for a single file, 2 to split plus (0) and minus (1).
    pub fn new(tracks: usize, n_refs: usize) -> Self {
        Coverage { blocks: vec![vec![Vec::new(); n_refs]; tracks] }
    }

    pub fn tracks(&self) -> usize {
        self.blocks.len()
    }

    /// Adds one alignment with the share of its fragment assigned to it (1
    /// unless the reassignment mode splits fragments). `blocks` are sorted by
    /// start; overlapping or touching ones are merged here so no base is
    /// counted twice.
    pub fn add(&mut self, track: usize, ref_id: usize, blocks: &[(i64, i64)], weight: f32) {
        let Some(dst) = self.blocks[track].get_mut(ref_id) else { return };
        let mut cur: Option<(i64, i64)> = None;
        for &(s, e) in blocks {
            match cur {
                Some((cs, ce)) if s <= ce => cur = Some((cs, ce.max(e))),
                Some((cs, ce)) => {
                    dst.push((cs.max(0) as u32, ce.max(0) as u32, weight));
                    cur = Some((s, e));
                }
                None => cur = Some((s, e)),
            }
        }
        if let Some((cs, ce)) = cur {
            dst.push((cs.max(0) as u32, ce.max(0) as u32, weight));
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
            .map(|(i, (start, end, depth))| (refs[i].0.clone(), Value { start, end, value: depth }));
        let runtime = tokio::runtime::Builder::new_current_thread().build()?;
        BigWigWrite::create_file(path, chrom_sizes)?
            .write(BedParserStreamingIterator::wrap_infallible_iter(intervals, true), runtime)
            .map_err(|e| io::Error::other(e.to_string()))?;
        Ok(true)
    }
}

/// Runs of constant, non-zero depth from weighted blocks: `(start, end, depth)`.
fn depth_intervals(blocks: Vec<(u32, u32, f32)>) -> Vec<(u32, u32, f32)> {
    // (position, change in depth, change in number of open blocks)
    let mut events: Vec<(u32, f64, i32)> = Vec::with_capacity(blocks.len() * 2);
    for &(s, e, w) in &blocks {
        events.push((s, w as f64, 1));
        events.push((e, -(w as f64), -1));
    }
    drop(blocks);
    events.sort_by(|a, b| a.0.cmp(&b.0).then(a.2.cmp(&b.2)));
    let mut out: Vec<(u32, u32, f32)> = Vec::new();
    let (mut depth, mut open, mut prev) = (0.0f64, 0i32, 0u32);
    let mut i = 0;
    while i < events.len() {
        let pos = events[i].0;
        if open > 0 && pos > prev {
            let d = depth as f32;
            match out.last_mut() {
                _ if d <= 0.0 => {} // weights too small to register
                Some(last) if last.1 == prev && last.2 == d => last.1 = pos,
                _ => out.push((prev, pos, d)),
            }
        }
        while i < events.len() && events[i].0 == pos {
            depth += events[i].1;
            open += events[i].2;
            i += 1;
        }
        if open == 0 {
            depth = 0.0; // clear rounding left over from fractional weights
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
        let got = depth_intervals(vec![(15, 30, 1.0), (10, 20, 1.0), (30, 40, 1.0), (15, 30, 1.0)]);
        assert_eq!(got, vec![(10, 15, 1.0), (15, 20, 3.0), (20, 30, 2.0), (30, 40, 1.0)]);
    }

    #[test]
    fn fractional_weights_add_up() {
        let got = depth_intervals(vec![(0, 10, 0.25), (5, 15, 0.5), (20, 30, 0.25)]);
        assert_eq!(got, vec![(0, 5, 0.25), (5, 10, 0.75), (10, 15, 0.5), (20, 30, 0.25)]);
    }

    #[test]
    fn mates_are_merged_before_counting() {
        let mut c = Coverage::new(1, 1);
        // mate blocks 100-150 and 130-180 (overlap), then a spliced block 300-320
        c.add(0, 0, &[(100, 150), (130, 180), (300, 320)], 1.0);
        assert_eq!(c.blocks[0][0], vec![(100, 180, 1.0), (300, 320, 1.0)]);
        assert_eq!(depth_intervals(c.blocks[0][0].clone()), vec![(100, 180, 1.0), (300, 320, 1.0)]);
    }
}
