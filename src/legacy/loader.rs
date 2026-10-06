//! Streaming BAM loader reproducing Telescope's `_load_sequential` and
//! `_mapping_to_matrix`.
//!
//! Memory stays proportional to the fragments that overlap the annotation:
//! each contributes one `(row, column, score)` triple per feature, and read
//! names are never kept (a 128-bit hash stands in for Telescope's
//! name -> row dictionary).

use std::collections::HashMap;
use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};
use std::io;

use noodles::{
    bam,
    sam::alignment::{
        RecordBuf,
        record::{cigar::op::Kind, data::field::Tag},
    },
};

use super::annotation::{Annotation, Scratch};
use super::model::ScoreMatrix;

/// `telescope.utils.BIG_INT`
const BIG_INT: i64 = (1 << 32) - 1;
/// Marks an alignment that overlaps no feature (Telescope's `no_feature_key`).
const NO_FEATURE: u32 = u32::MAX;

/// Library orientation for strand-aware assignment (`--stranded_mode`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Stranded {
    None,
    Rf,
    R,
    Fr,
    F,
}

impl Stranded {
    fn as_str(self) -> &'static str {
        match self {
            Stranded::None => "None",
            Stranded::Rf => "RF",
            Stranded::R => "R",
            Stranded::Fr => "FR",
            Stranded::F => "F",
        }
    }

    /// Fragment strand exactly as `Assigner._assign_pair_threshold` derives it
    /// (it indexes the mode *string*, including for `"None"`).
    fn frag_strand(self, r1_reversed: bool, is_paired: bool) -> u8 {
        let s = self.as_str().as_bytes();
        let is_f = if is_paired { s[s.len() - 1] == b'F' } else { s[0] == b'F' };
        match (r1_reversed, is_f) {
            (true, true) | (false, false) => b'+',
            _ => b'-',
        }
    }
}

/// Counters Telescope reports in the `## RunInfo` header line.
#[derive(Default, Debug, Clone)]
pub struct RunInfo {
    pub annotated_features: u64,
    pub total_fragments: u64,
    pub pair_mapped: u64,
    pub pair_mixed: u64,
    pub single_mapped: u64,
    pub unmapped: u64,
    pub unique: u64,
    pub ambig: u64,
    pub overlap_unique: u64,
    pub overlap_ambig: u64,
}

pub struct Loaded {
    pub matrix: ScoreMatrix,
    /// Feature name per matrix column; column 0 is the no-feature key.
    pub feat_names: Vec<String>,
    pub feat_lengths: Vec<u64>,
    pub info: RunInfo,
    /// Alignments whose best feature was one of several equally-overlapped loci.
    pub tied_alignments: u64,
}

#[derive(Default)]
struct Rec {
    flag: u16,
    ref_id: i32,
    start: i32,
    mate_ref: i32,
    mate_start: i32,
    tlen_abs: u32,
    score: Option<i64>,
    blocks: Vec<(i64, i64)>,
}

impl Rec {
    fn is_paired(&self) -> bool {
        self.flag & 0x1 != 0
    }
    fn is_proper(&self) -> bool {
        self.flag & 0x2 != 0
    }
    fn is_unmapped(&self) -> bool {
        self.flag & 0x4 != 0
    }
    fn is_reverse(&self) -> bool {
        self.flag & 0x10 != 0
    }
    fn is_read1(&self) -> bool {
        self.flag & 0x40 != 0
    }
    /// `alignment.readkey` minus the query name (constant within a bundle).
    fn readkey(&self) -> MateKey {
        (self.is_read1(), self.ref_id, self.start, self.mate_ref, self.mate_start, self.tlen_abs)
    }
    /// `alignment.matekey`: the key this record's mate would have.
    fn matekey(&self) -> MateKey {
        (!self.is_read1(), self.mate_ref, self.mate_start, self.ref_id, self.start, self.tlen_abs)
    }

    fn fill(&mut self, rec: &RecordBuf, as_tag: &Tag) {
        self.flag = u16::from(rec.flags());
        self.ref_id = rec.reference_sequence_id().map_or(-1, |i| i as i32);
        self.start = rec.alignment_start().map_or(-1, |p| usize::from(p) as i32 - 1);
        self.mate_ref = rec.mate_reference_sequence_id().map_or(-1, |i| i as i32);
        self.mate_start = rec.mate_alignment_start().map_or(-1, |p| usize::from(p) as i32 - 1);
        self.tlen_abs = rec.template_length().unsigned_abs();
        self.score = rec.data().get(as_tag).and_then(|v| v.as_int());
        // pysam's get_blocks(): gapless aligned runs, 0-based half-open.
        self.blocks.clear();
        let mut pos = self.start as i64;
        for op in rec.cigar().as_ref() {
            let len = op.len() as i64;
            match op.kind() {
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                    self.blocks.push((pos, pos + len));
                    pos += len;
                }
                Kind::Deletion | Kind::Skip => pos += len,
                _ => {}
            }
        }
    }
}

type MateKey = (bool, i32, i32, i32, i32, u32);

/// One of Telescope's `AlignedPair`s: indices into the current bundle.
#[derive(Clone, Copy)]
struct Aln {
    r1: usize,
    r2: Option<usize>,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Code {
    SingleUnmapped,
    SingleMapped,
    PairUnmapped,
    PairMapped,
    PairMixed,
}

/// `fetch_fragments_seq` for one bundle of same-named records. The whole
/// bundle is classified from its first record.
fn classify(recs: &[Rec], out: &mut Vec<Aln>, cache: &mut Vec<(MateKey, usize, bool)>) -> Code {
    out.clear();
    let first = &recs[0];
    let singles = |out: &mut Vec<Aln>| out.extend((0..recs.len()).map(|r1| Aln { r1, r2: None }));
    if !first.is_paired() {
        singles(out);
        return if first.is_unmapped() { Code::SingleUnmapped } else { Code::SingleMapped };
    }
    if !first.is_proper() {
        if recs.len() == 2 && recs.iter().all(Rec::is_unmapped) {
            out.push(Aln { r1: 0, r2: Some(1) });
            return Code::PairUnmapped;
        }
        singles(out);
        return Code::PairMixed;
    }
    // pair_bundle: match each record to a cached mate by key; leftovers are
    // emitted as singles in dict (insertion) order.
    cache.clear();
    for (i, r) in recs.iter().enumerate() {
        if !r.is_paired() {
            out.push(Aln { r1: i, r2: None });
            continue;
        }
        let want = r.matekey();
        if let Some(slot) = cache.iter_mut().find(|(k, _, alive)| *alive && *k == want) {
            slot.2 = false;
            let mate = slot.1;
            out.push(if r.is_read1() { Aln { r1: i, r2: Some(mate) } } else { Aln { r1: mate, r2: Some(i) } });
        } else {
            let key = r.readkey();
            match cache.iter_mut().find(|(k, _, alive)| *alive && *k == key) {
                Some(slot) => slot.1 = i,
                None => cache.push((key, i, true)),
            }
        }
    }
    out.extend(cache.iter().filter(|c| c.2).map(|c| Aln { r1: c.1, r2: None }));
    Code::PairMapped
}

/// `helpers.merge_blocks(blocks, 1)` in place.
fn merge_blocks(blocks: &mut Vec<(i64, i64)>) {
    if blocks.len() <= 1 {
        return;
    }
    blocks.sort_by_key(|b| b.0);
    let mut w = 0;
    for i in 1..blocks.len() {
        if blocks[i].0 - blocks[w].1 > 1 {
            w += 1;
            blocks[w] = blocks[i];
        } else {
            blocks[w].1 = blocks[w].1.max(blocks[i].1);
        }
    }
    blocks.truncate(w + 1);
}

fn name_hash(name: &[u8]) -> u128 {
    let mut a = DefaultHasher::new();
    name.hash(&mut a);
    let mut b = DefaultHasher::new();
    (0x9e37_79b9_7f4a_7c15u64, name).hash(&mut b);
    ((a.finish() as u128) << 64) | b.finish() as u128
}

fn invalid(msg: String) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg)
}

pub struct LoadOptions {
    pub overlap_threshold: f64,
    pub stranded: Stranded,
    pub no_feature_key: String,
}

struct Accum<'a> {
    annot: &'a Annotation,
    opts: &'a LoadOptions,
    ref_names: Vec<String>,
    info: RunInfo,
    min_as: i64,
    max_as: i64,
    nofeat: [u64; 2],
    feat: [u64; 2],
    tied: u64,
    rows: HashMap<u128, u32>,
    /// locus id -> matrix column (0 = unassigned; column 0 is no-feature)
    col_of_locus: Vec<u32>,
    col_locus: Vec<u32>,
    triples: Vec<(u32, u32, i32)>,
    // scratch
    alns: Vec<Aln>,
    cache: Vec<(MateKey, usize, bool)>,
    blocks: Vec<(i64, i64)>,
    hits: Scratch,
    scored: Vec<(u32, i64, i64)>,
    by_feat: Vec<(u32, i64, i64)>,
}

impl Accum<'_> {
    fn bundle(&mut self, name: &[u8], recs: &[Rec]) -> io::Result<()> {
        self.info.total_fragments += 1;
        let code = classify(recs, &mut self.alns, &mut self.cache);
        match code {
            Code::SingleUnmapped | Code::PairUnmapped => {
                self.info.unmapped += 1;
                return Ok(());
            }
            Code::SingleMapped => self.info.single_mapped += 1,
            Code::PairMapped => self.info.pair_mapped += 1,
            Code::PairMixed => self.info.pair_mixed += 1,
        }

        // (feature, alnscore, alnscore + alnlen) per mapped alignment
        self.scored.clear();
        let qname = || String::from_utf8_lossy(name).into_owned();
        for k in 0..self.alns.len() {
            let aln = self.alns[k];
            let r1 = &recs[aln.r1];
            if r1.is_unmapped() {
                continue;
            }
            let need = |r: &Rec| {
                r.score.ok_or_else(|| invalid(format!("read {}: alignment has no AS tag", qname())))
            };
            let mut score = need(r1)?;
            self.blocks.clear();
            self.blocks.extend_from_slice(&r1.blocks);
            if let Some(i2) = aln.r2 {
                score += need(&recs[i2])?;
                self.blocks.extend_from_slice(&recs[i2].blocks);
            }
            merge_blocks(&mut self.blocks);
            let alnlen: i64 = self.blocks.iter().map(|b| b.1 - b.0).sum();
            self.min_as = self.min_as.min(score);
            self.max_as = self.max_as.max(score);

            let strand = self.opts.stranded.frag_strand(r1.is_reverse(), aln.r2.is_some());
            let chrom = usize::try_from(r1.ref_id).ok().and_then(|i| self.ref_names.get(i));
            let hit = chrom.and_then(|c| self.annot.best_feature(c, &self.blocks, strand, &mut self.hits));
            let feat = match hit {
                Some(h) if h.overlap as f64 > alnlen as f64 * self.opts.overlap_threshold => {
                    self.tied += h.tied as u64;
                    h.locus
                }
                _ => NO_FEATURE,
            };
            self.scored.push((feat, score, score + alnlen));
        }
        if self.scored.is_empty() {
            return Err(invalid(format!(
                "read {}: no mapped alignment in a fragment classed as mapped (Telescope fails here too)",
                qname()
            )));
        }
        let ambig = (self.scored.len() > 1) as usize;
        if self.scored.iter().all(|s| s.0 == NO_FEATURE) {
            self.nofeat[ambig] += 1;
            return Ok(());
        }
        self.feat[ambig] += 1;

        // process_overlap_frag: best alignment per feature (first on ties),
        // then features ordered by that alignment's score, descending.
        self.by_feat.clear();
        for &(feat, score, total) in &self.scored {
            match self.by_feat.iter_mut().find(|e| e.0 == feat) {
                Some(e) if total > e.2 => *e = (feat, score, total),
                Some(_) => {}
                None => self.by_feat.push((feat, score, total)),
            }
        }
        self.by_feat.sort_by_key(|e| std::cmp::Reverse(e.1));

        let next_row = self.rows.len() as u32;
        let row = *self.rows.entry(name_hash(name)).or_insert(next_row);
        for &(feat, _, total) in &self.by_feat {
            let col = if feat == NO_FEATURE {
                0
            } else {
                if self.col_of_locus[feat as usize] == 0 {
                    self.col_locus.push(feat);
                    self.col_of_locus[feat as usize] = self.col_locus.len() as u32;
                }
                self.col_of_locus[feat as usize]
            };
            let total = i32::try_from(total).map_err(|_| invalid("alignment score out of range".into()))?;
            self.triples.push((row, col, total));
        }
        Ok(())
    }
}

pub fn load(bam_path: &str, annot: &Annotation, opts: &LoadOptions) -> io::Result<Loaded> {
    let mut reader = bam::io::reader::Builder.build_from_path(bam_path)?;
    let header = reader.read_header()?;
    let ref_names: Vec<String> = header
        .reference_sequences()
        .keys()
        .map(|k| String::from_utf8_lossy(k.as_ref()).into_owned())
        .collect();
    let as_tag = Tag::ALIGNMENT_SCORE;

    let mut acc = Accum {
        annot,
        opts,
        ref_names,
        info: RunInfo { annotated_features: annot.loci.len() as u64, ..Default::default() },
        min_as: BIG_INT,
        max_as: -BIG_INT,
        nofeat: [0; 2],
        feat: [0; 2],
        tied: 0,
        rows: HashMap::new(),
        col_of_locus: vec![0; annot.loci.len()],
        col_locus: Vec::new(),
        triples: Vec::new(),
        alns: Vec::new(),
        cache: Vec::new(),
        blocks: Vec::new(),
        hits: Scratch::default(),
        scored: Vec::new(),
        by_feat: Vec::new(),
    };

    // Records of the current bundle; slots are reused across bundles.
    let mut pool: Vec<Rec> = Vec::new();
    let mut n = 0usize;
    let mut cur_name: Vec<u8> = Vec::new();
    let mut record = RecordBuf::default();
    loop {
        if reader.read_record_buf(&header, &mut record)? == 0 {
            break;
        }
        let name: &[u8] = record.name().map_or(&[][..], |nm| nm.as_ref());
        if n > 0 && name != cur_name.as_slice() {
            acc.bundle(&cur_name, &pool[..n])?;
            n = 0;
        }
        if n == 0 {
            cur_name.clear();
            cur_name.extend_from_slice(name);
        }
        if n == pool.len() {
            pool.push(Rec::default());
        }
        pool[n].fill(&record, &as_tag);
        n += 1;
    }
    if n > 0 {
        acc.bundle(&cur_name, &pool[..n])?;
    }

    // _mapping_to_matrix: value = (AS - minAS + 1) + alnlen, max per cell.
    let Accum { mut info, min_as, nofeat, feat, tied, rows, col_locus, mut triples, .. } = acc;
    info.unique = nofeat[0] + feat[0];
    info.ambig = nofeat[1] + feat[1];
    let n_rows = rows.len();
    drop(rows);
    let n_cols = col_locus.len() + 1;

    triples.sort_unstable();
    let mut indptr = Vec::with_capacity(n_rows + 1);
    let mut indices: Vec<u32> = Vec::with_capacity(triples.len());
    let mut data: Vec<u16> = Vec::with_capacity(triples.len());
    indptr.push(0usize);
    let mut cur_row = 0u32;
    let mut last: Option<(u32, u32)> = None;
    for &(row, col, total) in &triples {
        let v = total as i64 - min_as + 1;
        let v = u16::try_from(v).map_err(|_| {
            invalid(format!("rescaled score {v} does not fit Telescope's uint16 score matrix"))
        })?;
        if last == Some((row, col)) {
            // sorted ascending, so a later duplicate is the larger score
            *data.last_mut().unwrap() = v;
            continue;
        }
        while cur_row < row {
            indptr.push(indices.len());
            cur_row += 1;
        }
        indices.push(col);
        data.push(v);
        last = Some((row, col));
    }
    while indptr.len() < n_rows + 1 {
        indptr.push(indices.len());
    }
    drop(triples);

    info.overlap_unique = indptr.windows(2).filter(|w| w[1] - w[0] == 1).count() as u64;
    info.overlap_ambig = n_rows as u64 - info.overlap_unique;

    let mut feat_names = Vec::with_capacity(n_cols);
    let mut feat_lengths = Vec::with_capacity(n_cols);
    feat_names.push(opts.no_feature_key.clone());
    feat_lengths.push(0);
    for &l in &col_locus {
        feat_names.push(annot.loci[l as usize].clone());
        feat_lengths.push(annot.lengths[l as usize]);
    }

    Ok(Loaded {
        matrix: ScoreMatrix { n_rows, n_cols, indptr, indices, data },
        feat_names,
        feat_lengths,
        info,
        tied_alignments: tied,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn merge_blocks_matches_helper() {
        let mut b = vec![(4, 9), (10, 14), (1, 3)];
        merge_blocks(&mut b);
        assert_eq!(b, vec![(1, 14)]);
        let mut b = vec![(10, 20), (22, 30)];
        merge_blocks(&mut b);
        assert_eq!(b, vec![(10, 20), (22, 30)]);
    }

    #[test]
    fn frag_strand_follows_reference_logic() {
        assert_eq!(Stranded::Rf.frag_strand(true, true), b'+');
        assert_eq!(Stranded::Rf.frag_strand(false, true), b'-');
        assert_eq!(Stranded::F.frag_strand(true, false), b'+');
        assert_eq!(Stranded::None.frag_strand(false, true), b'+');
    }
}
