use crate::RawRecord;
use crate::overlap::AnnotationOverlap;

/// A single alignment unit (either one unpaired read or a merged mate pair)
/// with lightweight metadata needed for analysis.
pub struct AlignmentEntry {
    pub qname: Option<String>,
    pub input_order: usize,
    /// Second BAM line index for paired mates (both get tagged in output).
    pub mate_input_order: Option<usize>,
    pub annotations: Vec<AnnotationOverlap>,
    pub alignment_score: Option<i64>,
    pub is_proper_pair: bool,
}

/// A group of alignment units sharing the same QNAME.
#[allow(dead_code)]
pub struct ReadGroup {
    pub qname: String,
    pub alignments: Vec<AlignmentEntry>,
}

/// Intermediate representation after mate pairing, before annotation.
/// Carries spans from one or two BAM records for merged annotation.
/// Spans are stored inline (max 2) to avoid per-unit heap allocation.
pub struct AlignmentUnit<'a> {
    pub qname: Option<String>,
    pub input_order: usize,
    pub mate_input_order: Option<usize>,
    pub alignment_score: Option<i64>,
    pub is_proper_pair: bool,
    /// Same convention as `AlignmentEntry::strand`.
    pub strand: Option<char>,
    spans_buf: [(&'a str, usize, usize); 2],
    span_count: u8,
}

impl<'a> AlignmentUnit<'a> {
    pub fn spans(&self) -> &[(&'a str, usize, usize)] {
        &self.spans_buf[..self.span_count as usize]
    }
}

/// Categorized pairing statistics.
pub struct PairingStats {
    pub valid_pairs: usize,
    pub unpaired: usize,
    pub flag_failures: usize,
    pub ref_failures: usize,
    pub pos_failures: usize,
}

impl PairingStats {
    fn new() -> Self {
        PairingStats {
            valid_pairs: 0,
            unpaired: 0,
            flag_failures: 0,
            ref_failures: 0,
            pos_failures: 0,
        }
    }

    pub fn total_failures(&self) -> usize {
        self.flag_failures + self.ref_failures + self.pos_failures
    }
}

/// Result of the mate pairing phase.
pub struct PairingResult<'a> {
    pub units: Vec<AlignmentUnit<'a>>,
    pub stats: PairingStats,
}

/// Sum two optional alignment scores.
fn sum_scores(a: Option<i64>, b: Option<i64>) -> Option<i64> {
    match (a, b) {
        (Some(x), Some(y)) => Some(x + y),
        (Some(x), None) | (None, Some(x)) => Some(x),
        (None, None) => None,
    }
}

/// Check whether two adjacent records can form a mate pair.
/// Requires: same QNAME, both paired (0x1), one first_segment + one last_segment.
fn can_pair(a: &RawRecord, b: &RawRecord) -> bool {
    if !a.is_paired || !b.is_paired {
        return false;
    }
    let a_name = a.qname.as_deref();
    let b_name = b.qname.as_deref();
    if a_name.is_none() || a_name != b_name {
        return false;
    }
    (a.is_first_segment && b.is_last_segment) || (a.is_last_segment && b.is_first_segment)
}

enum ValidationResult {
    Valid,
    FlagMismatch,
    RefMismatch,
    PosMismatch,
}

/// Validate a mate pair per spec.md:
/// - Flag correspondence: if read1 has mate_reverse (0x20), read2 should have reverse (0x10)
///   (skipped when either mate is unmapped, since strand bits are undefined)
/// - RNEXT/PNEXT of read1 should match RNAME/POS of read2, and vice versa
fn validate_pair(a: &RawRecord, b: &RawRecord) -> ValidationResult {
    // Flag correspondence: 0x20 on one ↔ 0x10 on the other
    // Skip for unmapped mates — strand bits are undefined per SAM spec
    let a_unmapped = a.ref_id.is_none();
    let b_unmapped = b.ref_id.is_none();
    if !a_unmapped
        && !b_unmapped
        && (!a.is_mate_reverse && b.is_reverse || !b.is_mate_reverse && a.is_reverse)
    {
        return ValidationResult::FlagMismatch;
    }

    // RNEXT/PNEXT matching
    if a.mate_ref_id != b.ref_id || b.mate_ref_id != a.ref_id {
        return ValidationResult::RefMismatch;
    }

    // Compare mate_alignment_start with the other's alignment start (both 0-based).
    // Uses alignment_start_0based (from POS) rather than span start (from CIGAR),
    // because records with CIGAR=* have no span but still have a valid POS.
    if a.mate_alignment_start != b.alignment_start_0based
        || b.mate_alignment_start != a.alignment_start_0based
    {
        return ValidationResult::PosMismatch;
    }

    ValidationResult::Valid
}

const EMPTY_SPAN: (&str, usize, usize) = ("", 0, 0);

/// The strand a mapped record aligns to: `'-'` if `0x10` is set, else `'+'`.
/// Returns `None` for records with no reference (unmapped).
fn record_strand(r: &RawRecord) -> Option<char> {
    r.ref_id?;
    Some(if r.is_reverse { '-' } else { '+' })
}

/// Strand for a single (unpaired) record's alignment unit.
///
/// Conventions:
/// - read1 (first segment): use the read's own mapped strand.
/// - read2-only (last segment): use the *reverse* of the read's mapped strand,
///   so the unit's strand reflects the implied read1 strand.
/// - neither bit set (true single-end): use the read's own mapped strand.
fn single_unit_strand(r: &RawRecord) -> Option<char> {
    let s = record_strand(r)?;
    if r.is_last_segment && !r.is_first_segment {
        Some(if s == '+' { '-' } else { '+' })
    } else {
        Some(s)
    }
}

/// Strand for a mated unit: prefer read1, else fall back to single_unit_strand
/// (which reverses read2 if needed).
fn mated_unit_strand(a: &RawRecord, b: &RawRecord) -> Option<char> {
    let r1 = if a.is_first_segment {
        Some(a)
    } else if b.is_first_segment {
        Some(b)
    } else {
        None
    };
    if let Some(r1) = r1 {
        return record_strand(r1);
    }
    // No first-segment mate present — derive from whichever last-segment we have.
    if a.is_last_segment {
        single_unit_strand(a)
    } else if b.is_last_segment {
        single_unit_strand(b)
    } else {
        record_strand(a).or_else(|| record_strand(b))
    }
}

/// Build an AlignmentUnit from a single RawRecord (unpaired).
fn single_unit<'a>(raw: RawRecord<'a>) -> AlignmentUnit<'a> {
    let (spans_buf, span_count) = match raw.span {
        Some(s) => ([s, EMPTY_SPAN], 1),
        None => ([EMPTY_SPAN, EMPTY_SPAN], 0),
    };
    let strand = single_unit_strand(&raw);
    AlignmentUnit {
        qname: raw.qname,
        input_order: raw.input_order,
        mate_input_order: None,
        alignment_score: raw.alignment_score,
        is_proper_pair: raw.is_proper_pair,
        strand,
        spans_buf,
        span_count,
    }
}

/// Build an AlignmentUnit from two mated RawRecords.
fn mated_unit<'a>(a: RawRecord<'a>, b: RawRecord<'a>) -> AlignmentUnit<'a> {
    let input_order = a.input_order.min(b.input_order);
    let mate_input_order = a.input_order.max(b.input_order);
    let alignment_score = sum_scores(a.alignment_score, b.alignment_score);
    let is_proper_pair = a.is_proper_pair; // both mates should agree
    let strand = mated_unit_strand(&a, &b);
    let qname = a.qname.or(b.qname);

    let mut spans_buf = [EMPTY_SPAN, EMPTY_SPAN];
    let mut span_count: u8 = 0;
    if let Some(s) = a.span {
        spans_buf[span_count as usize] = s;
        span_count += 1;
    }
    if let Some(s) = b.span {
        spans_buf[span_count as usize] = s;
        span_count += 1;
    }

    AlignmentUnit {
        qname,
        input_order,
        mate_input_order: Some(mate_input_order),
        alignment_score,
        is_proper_pair,
        strand,
        spans_buf,
        span_count,
    }
}

/// Pair adjacent mate records (PE mode) or pass through as singles (SE mode).
///
/// In PE mode, scans for adjacent records where one is first_segment and the
/// other is last_segment with the same QNAME. Validates flag correspondence
/// and RNEXT/PNEXT matching per the SAM spec. Failed validations produce a
/// warning and two unpaired units.
///
/// Note: this function was investigated for parallelization via chunk-based
/// splitting at safe pair boundaries. However, the per-record work (with
/// inline span storage in AlignmentUnit) is too cheap (~150ns/record) and
/// memory-bound. Chunking + merging overhead consistently exceeded the
/// parallelism benefit by ~50% on the long_JEG3_tester.bam benchmark.
///
/// In SE mode, every record becomes its own AlignmentUnit (no pairing attempted).
pub fn pair_mates(records: Vec<RawRecord>, single_end: bool) -> PairingResult {
    let mut units = Vec::with_capacity(records.len());
    let mut stats = PairingStats::new();
    let mut warned = false;

    let mut iter = records.into_iter().peekable();
    while let Some(record) = iter.next() {
        if !single_end
            && record.is_paired
            && !record.is_mate_unmapped
            && let Some(next) = iter.peek()
            && can_pair(&record, next)
        {
            let next = iter.next().unwrap();
            match validate_pair(&record, &next) {
                ValidationResult::Valid => {
                    stats.valid_pairs += 1;
                    units.push(mated_unit(record, next));
                }
                failure => {
                    match failure {
                        ValidationResult::FlagMismatch => stats.flag_failures += 1,
                        ValidationResult::RefMismatch => stats.ref_failures += 1,
                        ValidationResult::PosMismatch => stats.pos_failures += 1,
                        ValidationResult::Valid => unreachable!(),
                    }
                    if !warned {
                        eprintln!(
                            "Warning: mate pair validation failed for '{}' — \
                             treating as unpaired (flag/position mismatch)",
                            record.qname.as_deref().unwrap_or("?")
                        );
                        warned = true;
                    }
                    units.push(single_unit(record));
                    units.push(single_unit(next));
                }
            }
            continue;
        }
        stats.unpaired += 1;
        units.push(single_unit(record));
    }

    PairingResult { units, stats }
}

/// Group alignment entries by QNAME using a single-pass run-length scan.
///
/// Relies on the invariant that entries with the same QNAME are adjacent in
/// `entries`. This holds because the tool only runs on aligner-native or
/// name-sorted BAMs (the same invariant `pair_mates` depends on), and
/// `pair_mates` preserves input order. Start a new `ReadGroup` whenever the
/// QNAME differs from the previously opened group.
///
/// Entries with no QNAME are skipped with a warning to stderr.
pub fn group_by_qname(entries: Vec<AlignmentEntry>) -> Vec<ReadGroup> {
    let mut groups: Vec<ReadGroup> = Vec::new();
    let mut warned_no_name = false;

    for mut entry in entries {
        let qname = match entry.qname.take() {
            Some(name) => name,
            None => {
                if !warned_no_name {
                    eprintln!("Warning: skipping record(s) with no QNAME");
                    warned_no_name = true;
                }
                continue;
            }
        };

        match groups.last_mut() {
            Some(g) if g.qname == qname => g.alignments.push(entry),
            _ => groups.push(ReadGroup {
                qname,
                alignments: vec![entry],
            }),
        }
    }

    groups
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_entry(name: Option<&str>, input_order: usize) -> AlignmentEntry {
        AlignmentEntry {
            qname: name.map(|n| n.to_string()),
            input_order,
            mate_input_order: None,
            annotations: vec![],
            alignment_score: None,
            is_proper_pair: false,
        }
    }

    #[test]
    fn test_group_by_qname() {
        let entries = vec![
            make_entry(Some("read1"), 0),
            make_entry(Some("read1"), 1),
            make_entry(Some("read2"), 2),
        ];
        let groups = group_by_qname(entries);
        assert_eq!(groups.len(), 2);
        assert_eq!(groups[0].qname, "read1");
        assert_eq!(groups[0].alignments.len(), 2);
        assert_eq!(groups[1].qname, "read2");
        assert_eq!(groups[1].alignments.len(), 1);
    }

    #[test]
    fn test_group_skips_no_name() {
        let entries = vec![make_entry(None, 0), make_entry(Some("read1"), 1)];
        let groups = group_by_qname(entries);
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].qname, "read1");
    }

    fn make_raw(
        input_order: usize,
        is_first: bool,
        is_last: bool,
        is_reverse: bool,
        mapped: bool,
    ) -> RawRecord<'static> {
        RawRecord {
            qname: Some("r".to_string()),
            input_order,
            alignment_score: None,
            is_proper_pair: true,
            span: None,
            is_paired: true,
            is_first_segment: is_first,
            is_last_segment: is_last,
            is_reverse,
            is_mate_reverse: false,
            is_mate_unmapped: false,
            mate_ref_id: None,
            mate_alignment_start: None,
            ref_id: if mapped { Some(0) } else { None },
            alignment_start_0based: None,
        }
    }

    #[test]
    fn test_single_unit_strand_se_forward() {
        // SE record (no first/last bits): strand == read's mapped strand.
        let r = make_raw(0, false, false, false, true);
        assert_eq!(single_unit_strand(&r), Some('+'));
    }

    #[test]
    fn test_single_unit_strand_se_reverse() {
        let r = make_raw(0, false, false, true, true);
        assert_eq!(single_unit_strand(&r), Some('-'));
    }

    #[test]
    fn test_single_unit_strand_read1_uses_own() {
        // First-segment record: use its own strand.
        let r = make_raw(0, true, false, false, true);
        assert_eq!(single_unit_strand(&r), Some('+'));
        let r = make_raw(0, true, false, true, true);
        assert_eq!(single_unit_strand(&r), Some('-'));
    }

    #[test]
    fn test_single_unit_strand_read2_only_inverts() {
        // Last-segment-only record (read1 missing): unit strand is the reverse.
        let r = make_raw(0, false, true, true, true);
        assert_eq!(single_unit_strand(&r), Some('+'));
        let r = make_raw(0, false, true, false, true);
        assert_eq!(single_unit_strand(&r), Some('-'));
    }

    #[test]
    fn test_single_unit_strand_unmapped() {
        let r = make_raw(0, true, false, false, false);
        assert_eq!(single_unit_strand(&r), None);
    }

    #[test]
    fn test_mated_unit_strand_uses_read1() {
        // read1 forward, read2 reverse → unit strand = '+'
        let r1 = make_raw(0, true, false, false, true);
        let r2 = make_raw(1, false, true, true, true);
        assert_eq!(mated_unit_strand(&r1, &r2), Some('+'));
        // Argument order should not matter.
        assert_eq!(mated_unit_strand(&r2, &r1), Some('+'));
    }

    #[test]
    fn test_sum_scores() {
        assert_eq!(sum_scores(Some(100), Some(50)), Some(150));
        assert_eq!(sum_scores(Some(100), None), Some(100));
        assert_eq!(sum_scores(None, Some(50)), Some(50));
        assert_eq!(sum_scores(None, None), None);
    }
}
