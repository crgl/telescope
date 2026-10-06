use std::cmp;
use std::collections::HashMap;

use noodles::sam::{Header, alignment::RecordBuf};

use crate::gtf::FeatureIndex;
use crate::intern::Symbol;

/// Result of annotating a BAM record against the GTF feature index.
pub const NO_FEATURE: &str = "__no_feature__";

/// Per-annotation overlap info for a single alignment.
#[derive(Debug, Clone)]
pub struct AnnotationOverlap {
    pub annotation: Symbol,
    pub overlap_bp: usize,
}

/// Strand-compatibility filter applied during overlap.
///
/// Built from the `--stranded` library mode and the alignment-unit's strand:
/// FR/F libraries use `Same(s_r)` (annotation must be on the same strand as the
/// read), RF/R use `Opposite(s_r)`. Annotations whose own strand is unknown
/// (`None`) always pass.
#[derive(Debug, Clone, Copy)]
pub enum StrandFilter {
    Same(char),
    Opposite(char),
}

impl StrandFilter {
    /// Whether a feature with the given strand is compatible.
    /// Features with unknown strand are always retained.
    fn accepts(&self, feature_strand: Option<char>) -> bool {
        let Some(s) = feature_strand else { return true };
        match *self {
            StrandFilter::Same(r) => s == r,
            StrandFilter::Opposite(r) => s != r,
        }
    }
}

/// Extract reference span from a BAM record.
/// Returns (chrom_name, start, stop) in 0-based half-open coordinates [start, stop)
/// matching rust-lapper's coordinate system.
/// Returns None for unmapped reads or records missing position/CIGAR.
/// Borrows chrom name from the header to avoid per-record allocation.
pub fn reference_span<'a>(
    record: &RecordBuf,
    header: &'a Header,
) -> Option<(&'a str, usize, usize)> {
    let ref_id = record.reference_sequence_id()?;
    let alignment_start = record.alignment_start()?;
    let span = record.cigar().alignment_span();

    if span == 0 {
        return None;
    }

    // Get chromosome name from header
    let (name, _) = header.reference_sequences().get_index(ref_id)?;
    let chrom = std::str::from_utf8(name.as_ref()).ok()?;

    // alignment_start is 1-based Position. Convert to 0-based for rust-lapper.
    let start_0based = usize::from(alignment_start) - 1;
    let stop = start_0based + span;

    Some((chrom, start_0based, stop))
}

/// Annotate a genomic span against the GTF feature index.
/// Returns all overlapping annotations with their summed overlap amounts.
/// Takes pre-extracted coordinates instead of a RecordBuf, enabling parallel use.
///
/// If `strand_filter` is `Some`, features whose GTF strand is incompatible
/// are excluded. Features with unknown strand always pass.
///
/// `min_overlap` filters the per-annotation overlap *after* exons sharing a
/// symbol have been summed, so a read crossing two short exons of the same
/// gene can still pass even if no individual exon overlap exceeds the
/// threshold. Set to `0` to disable.
#[allow(clippy::too_many_arguments)]
pub fn annotate_from_span(
    chrom: &str,
    start: usize,
    stop: usize,
    index: &FeatureIndex,
    field: &str,
    normalize_chr: bool,
    strand_filter: Option<StrandFilter>,
    min_overlap: usize,
) -> Vec<AnnotationOverlap> {
    let hits = index.find(chrom, start, stop, normalize_chr);

    if hits.is_empty() {
        return Vec::new();
    }

    let mut overlap_by_symbol: HashMap<Symbol, usize> = HashMap::new();

    for iv in &hits {
        let overlap = cmp::min(stop, iv.stop).saturating_sub(cmp::max(start, iv.start));
        if overlap == 0 {
            continue;
        }
        if let Some(filter) = strand_filter
            && !filter.accepts(iv.val.strand)
        {
            continue;
        }

        if let Some(&sym) = iv.val.attributes.get(field) {
            *overlap_by_symbol.entry(sym).or_insert(0) += overlap;
        }
    }

    overlap_by_symbol
        .into_iter()
        .filter(|(_, bp)| *bp >= min_overlap)
        .map(|(annotation, overlap_bp)| AnnotationOverlap {
            annotation,
            overlap_bp,
        })
        .collect()
}

/// Annotate from multiple spans (e.g., a mate pair), avoiding double-counting
/// overlapping regions. For same-chromosome spans, uses inclusion-exclusion:
///   overlap(feature, span1 ∪ span2) = overlap(feature, span1) + overlap(feature, span2)
///                                     - overlap(feature, span1 ∩ span2)
/// For different-chromosome spans, overlaps are summed independently (no intersection).
///
/// `min_overlap` is applied to the final combined per-symbol overlap. The
/// inner per-span calls use `min_overlap = 0` so a pair where each mate
/// individually falls below the threshold but together exceeds it still
/// passes (e.g. 20bp + 15bp at the threshold of 30bp passes).
#[allow(clippy::too_many_arguments)]
pub fn annotate_from_spans(
    spans: &[(&str, usize, usize)],
    index: &FeatureIndex,
    field: &str,
    normalize_chr: bool,
    strand_filter: Option<StrandFilter>,
    min_overlap: usize,
) -> Vec<AnnotationOverlap> {
    match spans.len() {
        0 => Vec::new(),
        1 => annotate_from_span(
            spans[0].0,
            spans[0].1,
            spans[0].2,
            index,
            field,
            normalize_chr,
            strand_filter,
            min_overlap,
        ),
        _ => {
            let (chrom1, start1, stop1) = spans[0];
            let (chrom2, start2, stop2) = spans[1];

            // Get overlaps for each span independently. Use min_overlap = 0
            // here; we apply the user's threshold to the combined value below.
            let overlaps1 = annotate_from_span(
                chrom1,
                start1,
                stop1,
                index,
                field,
                normalize_chr,
                strand_filter,
                0,
            );
            let overlaps2 = annotate_from_span(
                chrom2,
                start2,
                stop2,
                index,
                field,
                normalize_chr,
                strand_filter,
                0,
            );

            // Merge into a single map
            let mut overlap_by_symbol: HashMap<Symbol, usize> = HashMap::new();
            for ao in &overlaps1 {
                *overlap_by_symbol.entry(ao.annotation).or_insert(0) += ao.overlap_bp;
            }
            for ao in &overlaps2 {
                *overlap_by_symbol.entry(ao.annotation).or_insert(0) += ao.overlap_bp;
            }

            // Subtract double-counted region for same-chromosome overlapping spans
            if chrom1 == chrom2 {
                let inter_start = cmp::max(start1, start2);
                let inter_stop = cmp::min(stop1, stop2);
                if inter_start < inter_stop {
                    // The spans overlap — subtract the intersection's contribution
                    let intersection_overlaps = annotate_from_span(
                        chrom1,
                        inter_start,
                        inter_stop,
                        index,
                        field,
                        normalize_chr,
                        strand_filter,
                        0,
                    );
                    for ao in &intersection_overlaps {
                        if let Some(val) = overlap_by_symbol.get_mut(&ao.annotation) {
                            *val = val.saturating_sub(ao.overlap_bp);
                        }
                    }
                }
            }

            // Final threshold filter on the combined per-annotation overlap.
            overlap_by_symbol
                .into_iter()
                .filter(|(_, bp)| *bp >= min_overlap.max(1))
                .map(|(annotation, overlap_bp)| AnnotationOverlap {
                    annotation,
                    overlap_bp,
                })
                .collect()
        }
    }
}

/// Returns ALL annotations a record overlaps, with their summed overlap amounts.
/// Overlaps from multiple exons sharing the same attribute value are summed.
pub fn annotate_record_all(
    record: &RecordBuf,
    header: &Header,
    index: &FeatureIndex,
    field: &str,
    normalize_chr: bool,
) -> Vec<AnnotationOverlap> {
    let Some((chrom, start, stop)) = reference_span(record, header) else {
        return Vec::new();
    };

    annotate_from_span(chrom, start, stop, index, field, normalize_chr, None, 0)
}

/// Annotate a BAM record with the best-overlapping GTF feature.
///
/// Returns the resolved annotation string from the feature with the
/// greatest total overlap, or `__no_feature__` if no features overlap.
#[allow(dead_code)]
pub fn annotate_record(
    record: &RecordBuf,
    header: &Header,
    index: &FeatureIndex,
    field: &str,
    normalize_chr: bool,
) -> String {
    let overlaps = annotate_record_all(record, header, index, field, normalize_chr);

    overlaps
        .into_iter()
        .max_by_key(|o| o.overlap_bp)
        .map(|o| index.interner.resolve(o.annotation).to_string())
        .unwrap_or_else(|| NO_FEATURE.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gtf::FeatureIndex;
    use noodles::core::Position;
    use noodles::sam::{
        Header,
        alignment::{
            RecordBuf,
            record::cigar::{Op, op::Kind},
            record_buf::Cigar,
        },
        header::record::value::{Map, map::ReferenceSequence},
    };

    fn build_header_with_ref(name: &str, len: usize) -> Header {
        let mut header = Header::default();
        let map = Map::<ReferenceSequence>::new(std::num::NonZeroUsize::new(len).unwrap());
        header.reference_sequences_mut().insert(name.into(), map);
        header
    }

    fn build_record(ref_id: usize, start_1based: usize, cigar_ops: Vec<Op>) -> RecordBuf {
        let cigar: Cigar = cigar_ops.into_iter().collect();
        RecordBuf::builder()
            .set_reference_sequence_id(ref_id)
            .set_alignment_start(Position::new(start_1based).unwrap())
            .set_cigar(cigar)
            .build()
    }

    #[test]
    fn test_annotate_overlap() {
        let gtf = "chr1\tsrc\texon\t100\t200\t.\t+\t.\tgene_id \"GENE_A\"";
        let index = FeatureIndex::from_str(gtf, "gene_id");
        let header = build_header_with_ref("chr1", 1000);

        let record = build_record(0, 150, vec![Op::new(Kind::Match, 100)]);
        let result = annotate_record(&record, &header, &index, "gene_id", false);
        assert_eq!(result, "GENE_A");
    }

    #[test]
    fn test_strand_filter_keeps_compatible() {
        let gtf = "\
chr1\tsrc\texon\t100\t200\t.\t+\t.\tgene_id \"PLUS\"
chr1\tsrc\texon\t100\t200\t.\t-\t.\tgene_id \"MINUS\"";
        let index = FeatureIndex::from_str(gtf, "gene_id");

        // Same('+'): keep PLUS, drop MINUS.
        let hits = annotate_from_span(
            "chr1",
            120,
            180,
            &index,
            "gene_id",
            false,
            Some(StrandFilter::Same('+')),
            0,
        );
        assert_eq!(hits.len(), 1);
        assert_eq!(index.interner.resolve(hits[0].annotation), "PLUS");

        // Opposite('+'): keep MINUS, drop PLUS.
        let hits = annotate_from_span(
            "chr1",
            120,
            180,
            &index,
            "gene_id",
            false,
            Some(StrandFilter::Opposite('+')),
            0,
        );
        assert_eq!(hits.len(), 1);
        assert_eq!(index.interner.resolve(hits[0].annotation), "MINUS");

        // No filter: both.
        let hits = annotate_from_span("chr1", 120, 180, &index, "gene_id", false, None, 0);
        assert_eq!(hits.len(), 2);
    }

    #[test]
    fn test_min_overlap_drops_short() {
        let gtf = "chr1\tsrc\texon\t100\t200\t.\t+\t.\tgene_id \"GENE_A\"";
        let index = FeatureIndex::from_str(gtf, "gene_id");

        // Read [0, 120) overlaps GENE_A by 21 bp (positions 99..120).
        // min_overlap = 30 → dropped. min_overlap = 20 → kept.
        let hits = annotate_from_span("chr1", 0, 120, &index, "gene_id", false, None, 30);
        assert_eq!(hits.len(), 0);
        let hits = annotate_from_span("chr1", 0, 120, &index, "gene_id", false, None, 20);
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].overlap_bp, 21);
    }

    #[test]
    fn test_min_overlap_summed_across_exons_passes() {
        // Two short exons of the same gene; an alignment spanning both should
        // pass min_overlap=30 even though no individual exon overlap reaches it.
        let gtf = "\
chr1\tsrc\texon\t100\t120\t.\t+\t.\tgene_id \"GENE_A\"
chr1\tsrc\texon\t150\t170\t.\t+\t.\tgene_id \"GENE_A\"";
        let index = FeatureIndex::from_str(gtf, "gene_id");
        // Read [99, 170): hits both exons. Exon 1 = 21 bp, exon 2 = 21 bp.
        // Combined = 42, each individually = 21.
        let hits = annotate_from_span("chr1", 99, 170, &index, "gene_id", false, None, 30);
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].overlap_bp, 42);
    }

    #[test]
    fn test_strand_filter_unknown_annotation_passes() {
        let gtf = "chr1\tsrc\texon\t100\t200\t.\t.\t.\tgene_id \"UNK\"";
        let index = FeatureIndex::from_str(gtf, "gene_id");

        for filter in [StrandFilter::Same('+'), StrandFilter::Opposite('+')] {
            let hits =
                annotate_from_span("chr1", 120, 180, &index, "gene_id", false, Some(filter), 0);
            assert_eq!(hits.len(), 1, "unknown strand should always pass");
        }
    }

    #[test]
    fn test_annotate_no_overlap() {
        let gtf = "chr1\tsrc\texon\t100\t200\t.\t+\t.\tgene_id \"GENE_A\"";
        let index = FeatureIndex::from_str(gtf, "gene_id");
        let header = build_header_with_ref("chr1", 1000);

        let record = build_record(0, 300, vec![Op::new(Kind::Match, 100)]);
        let result = annotate_record(&record, &header, &index, "gene_id", false);
        assert_eq!(result, NO_FEATURE);
    }

    #[test]
    fn test_annotate_unmapped() {
        let gtf = "chr1\tsrc\texon\t100\t200\t.\t+\t.\tgene_id \"GENE_A\"";
        let index = FeatureIndex::from_str(gtf, "gene_id");
        let header = build_header_with_ref("chr1", 1000);

        let record = RecordBuf::default();
        let result = annotate_record(&record, &header, &index, "gene_id", false);
        assert_eq!(result, NO_FEATURE);
    }

    #[test]
    fn test_annotate_summed_overlap() {
        let gtf = "\
chr1\tsrc\texon\t100\t200\t.\t+\t.\tgene_id \"GENE_A\"
chr1\tsrc\texon\t300\t400\t.\t+\t.\tgene_id \"GENE_A\"
chr1\tsrc\texon\t150\t350\t.\t+\t.\tgene_id \"GENE_B\"";

        let index = FeatureIndex::from_str(gtf, "gene_id");
        let header = build_header_with_ref("chr1", 1000);

        let record = build_record(0, 100, vec![Op::new(Kind::Match, 301)]);
        let result = annotate_record(&record, &header, &index, "gene_id", false);
        assert_eq!(result, "GENE_A");
    }

    #[test]
    fn test_annotate_missing_attribute() {
        let gtf = "chr1\tsrc\texon\t100\t200\t.\t+\t.\tgene_id \"GENE_A\"";
        let index = FeatureIndex::from_str(gtf, "gene_id");
        let header = build_header_with_ref("chr1", 1000);

        let record = build_record(0, 150, vec![Op::new(Kind::Match, 20)]);
        let result = annotate_record(&record, &header, &index, "nonexistent_field", false);
        assert_eq!(result, NO_FEATURE);
    }

    #[test]
    fn test_annotate_record_all_multi_gene() {
        let gtf = "\
chr1\tsrc\texon\t100\t200\t.\t+\t.\tgene_id \"GENE_A\"
chr1\tsrc\texon\t150\t250\t.\t+\t.\tgene_id \"GENE_B\"";

        let index = FeatureIndex::from_str(gtf, "gene_id");
        let header = build_header_with_ref("chr1", 1000);

        // Read spans [150, 219] (1-based, 70bp) = [149, 219) 0-based
        // GENE_A [99, 200): overlap = min(219,200) - max(149,99) = 200-149 = 51
        // GENE_B [149, 250): overlap = min(219,250) - max(149,149) = 219-149 = 70
        let record = build_record(0, 150, vec![Op::new(Kind::Match, 70)]);
        let overlaps = annotate_record_all(&record, &header, &index, "gene_id", false);

        assert_eq!(overlaps.len(), 2);
        let a = overlaps
            .iter()
            .find(|o| index.interner.resolve(o.annotation) == "GENE_A")
            .unwrap();
        let b = overlaps
            .iter()
            .find(|o| index.interner.resolve(o.annotation) == "GENE_B")
            .unwrap();
        assert_eq!(a.overlap_bp, 51);
        assert_eq!(b.overlap_bp, 70);
    }

    #[test]
    fn test_annotate_record_all_no_overlap() {
        let gtf = "chr1\tsrc\texon\t100\t200\t.\t+\t.\tgene_id \"GENE_A\"";
        let index = FeatureIndex::from_str(gtf, "gene_id");
        let header = build_header_with_ref("chr1", 1000);

        let record = build_record(0, 300, vec![Op::new(Kind::Match, 50)]);
        let overlaps = annotate_record_all(&record, &header, &index, "gene_id", false);
        assert!(overlaps.is_empty());
    }
}
