use std::fs;
use std::path::PathBuf;
use std::process::Command;

use noodles::{
    bam,
    core::Position,
    sam::{
        Header,
        alignment::{
            RecordBuf,
            io::Write as AlignmentWrite,
            record::Flags,
            record::cigar::{Op, op::Kind},
            record::data::field::Tag,
            record_buf::{Cigar, data::field::Value},
        },
        header::record::value::{Map, map::ReferenceSequence},
    },
};

/// Build a header with a single reference sequence.
fn build_header(ref_name: &str, ref_len: usize) -> Header {
    let mut header = Header::default();
    let map = Map::<ReferenceSequence>::new(std::num::NonZeroUsize::new(ref_len).unwrap());
    header
        .reference_sequences_mut()
        .insert(ref_name.into(), map);
    header
}

/// Build a header with two reference sequences.
fn build_header_2refs(name1: &str, len1: usize, name2: &str, len2: usize) -> Header {
    let mut header = Header::default();
    let map1 = Map::<ReferenceSequence>::new(std::num::NonZeroUsize::new(len1).unwrap());
    header.reference_sequences_mut().insert(name1.into(), map1);
    let map2 = Map::<ReferenceSequence>::new(std::num::NonZeroUsize::new(len2).unwrap());
    header.reference_sequences_mut().insert(name2.into(), map2);
    header
}

/// Build a proper mate pair (read1 + read2) with correct flags and RNEXT/PNEXT.
/// Both mates map to the same reference. Returns (read1_record, read2_record).
fn mate_pair(
    name: &str,
    ref_id: usize,
    start1_1based: usize,
    len1: usize,
    start2_1based: usize,
    len2: usize,
) -> (RecordBuf, RecordBuf) {
    // read1: forward strand, mate on reverse strand
    let flags1 = Flags::SEGMENTED
        | Flags::PROPERLY_SEGMENTED
        | Flags::MATE_REVERSE_COMPLEMENTED
        | Flags::FIRST_SEGMENT;
    // read2: reverse strand, mate on forward strand
    let flags2 = Flags::SEGMENTED
        | Flags::PROPERLY_SEGMENTED
        | Flags::REVERSE_COMPLEMENTED
        | Flags::LAST_SEGMENT;

    let cigar1: Cigar = vec![Op::new(Kind::Match, len1)].into_iter().collect();
    let cigar2: Cigar = vec![Op::new(Kind::Match, len2)].into_iter().collect();

    let r1 = RecordBuf::builder()
        .set_name(name)
        .set_flags(flags1)
        .set_reference_sequence_id(ref_id)
        .set_alignment_start(Position::new(start1_1based).unwrap())
        .set_cigar(cigar1)
        .set_mate_reference_sequence_id(ref_id)
        .set_mate_alignment_start(Position::new(start2_1based).unwrap())
        .build();

    let r2 = RecordBuf::builder()
        .set_name(name)
        .set_flags(flags2)
        .set_reference_sequence_id(ref_id)
        .set_alignment_start(Position::new(start2_1based).unwrap())
        .set_cigar(cigar2)
        .set_mate_reference_sequence_id(ref_id)
        .set_mate_alignment_start(Position::new(start1_1based).unwrap())
        .build();

    (r1, r2)
}

/// Build a proper mate pair with AS tags on both mates.
#[allow(clippy::too_many_arguments)]
fn mate_pair_with_as(
    name: &str,
    ref_id: usize,
    start1_1based: usize,
    len1: usize,
    start2_1based: usize,
    len2: usize,
    as1: i32,
    as2: i32,
) -> (RecordBuf, RecordBuf) {
    let (mut r1, mut r2) = mate_pair(name, ref_id, start1_1based, len1, start2_1based, len2);
    r1.data_mut()
        .insert(Tag::new(b'A', b'S'), Value::Int32(as1));
    r2.data_mut()
        .insert(Tag::new(b'A', b'S'), Value::Int32(as2));
    (r1, r2)
}

/// Build a discordant mate pair where mates map to different references.
fn discordant_mate_pair(
    name: &str,
    ref_id1: usize,
    start1_1based: usize,
    len1: usize,
    ref_id2: usize,
    start2_1based: usize,
    len2: usize,
) -> (RecordBuf, RecordBuf) {
    let flags1 = Flags::SEGMENTED | Flags::MATE_REVERSE_COMPLEMENTED | Flags::FIRST_SEGMENT;
    let flags2 = Flags::SEGMENTED | Flags::REVERSE_COMPLEMENTED | Flags::LAST_SEGMENT;

    let cigar1: Cigar = vec![Op::new(Kind::Match, len1)].into_iter().collect();
    let cigar2: Cigar = vec![Op::new(Kind::Match, len2)].into_iter().collect();

    let r1 = RecordBuf::builder()
        .set_name(name)
        .set_flags(flags1)
        .set_reference_sequence_id(ref_id1)
        .set_alignment_start(Position::new(start1_1based).unwrap())
        .set_cigar(cigar1)
        .set_mate_reference_sequence_id(ref_id2)
        .set_mate_alignment_start(Position::new(start2_1based).unwrap())
        .build();

    let r2 = RecordBuf::builder()
        .set_name(name)
        .set_flags(flags2)
        .set_reference_sequence_id(ref_id2)
        .set_alignment_start(Position::new(start2_1based).unwrap())
        .set_cigar(cigar2)
        .set_mate_reference_sequence_id(ref_id1)
        .set_mate_alignment_start(Position::new(start1_1based).unwrap())
        .build();

    (r1, r2)
}

/// Build an aligned record with a name, at a given 1-based position.
fn named_record(name: &str, ref_id: usize, start_1based: usize, match_len: usize) -> RecordBuf {
    let cigar: Cigar = vec![Op::new(Kind::Match, match_len)].into_iter().collect();
    RecordBuf::builder()
        .set_name(name)
        .set_flags(Flags::SEGMENTED)
        .set_reference_sequence_id(ref_id)
        .set_alignment_start(Position::new(start_1based).unwrap())
        .set_cigar(cigar)
        .build()
}

/// Build an unmapped record (RNAME = '*', no reference sequence id).
fn unmapped_record(name: &str) -> RecordBuf {
    RecordBuf::builder()
        .set_name(name)
        .set_flags(Flags::UNMAPPED)
        .build()
}

/// Build a "placed unmapped" record: the 0x4 unmapped flag is set and there is
/// no alignment (CIGAR '*'), but RNAME/POS are populated (typically the mapped
/// mate's coordinates). Such a record has a reference id, so it is dropped only
/// by the unmapped flag, not by an RNAME='*' check.
fn placed_unmapped_record(name: &str, ref_id: usize, start_1based: usize) -> RecordBuf {
    RecordBuf::builder()
        .set_name(name)
        .set_flags(Flags::SEGMENTED | Flags::UNMAPPED | Flags::LAST_SEGMENT)
        .set_reference_sequence_id(ref_id)
        .set_alignment_start(Position::new(start_1based).unwrap())
        .build()
}

/// Build a SE record on a specific strand.
fn named_record_strand(
    name: &str,
    ref_id: usize,
    start_1based: usize,
    match_len: usize,
    reverse: bool,
) -> RecordBuf {
    let mut flags = Flags::SEGMENTED;
    if reverse {
        flags |= Flags::REVERSE_COMPLEMENTED;
    }
    let cigar: Cigar = vec![Op::new(Kind::Match, match_len)].into_iter().collect();
    RecordBuf::builder()
        .set_name(name)
        .set_flags(flags)
        .set_reference_sequence_id(ref_id)
        .set_alignment_start(Position::new(start_1based).unwrap())
        .set_cigar(cigar)
        .build()
}

/// Build a mate pair where the read1/read2 strands can be flipped independently.
/// Keeps the SAM-spec invariant that `0x20` on one mate matches `0x10` on the other.
#[allow(clippy::too_many_arguments)]
fn mate_pair_strand(
    name: &str,
    ref_id: usize,
    start1_1based: usize,
    len1: usize,
    start2_1based: usize,
    len2: usize,
    r1_reverse: bool,
    r2_reverse: bool,
) -> (RecordBuf, RecordBuf) {
    let mut flags1 = Flags::SEGMENTED | Flags::PROPERLY_SEGMENTED | Flags::FIRST_SEGMENT;
    if r1_reverse {
        flags1 |= Flags::REVERSE_COMPLEMENTED;
    }
    if r2_reverse {
        flags1 |= Flags::MATE_REVERSE_COMPLEMENTED;
    }
    let mut flags2 = Flags::SEGMENTED | Flags::PROPERLY_SEGMENTED | Flags::LAST_SEGMENT;
    if r2_reverse {
        flags2 |= Flags::REVERSE_COMPLEMENTED;
    }
    if r1_reverse {
        flags2 |= Flags::MATE_REVERSE_COMPLEMENTED;
    }
    let cigar1: Cigar = vec![Op::new(Kind::Match, len1)].into_iter().collect();
    let cigar2: Cigar = vec![Op::new(Kind::Match, len2)].into_iter().collect();
    let r1 = RecordBuf::builder()
        .set_name(name)
        .set_flags(flags1)
        .set_reference_sequence_id(ref_id)
        .set_alignment_start(Position::new(start1_1based).unwrap())
        .set_cigar(cigar1)
        .set_mate_reference_sequence_id(ref_id)
        .set_mate_alignment_start(Position::new(start2_1based).unwrap())
        .build();
    let r2 = RecordBuf::builder()
        .set_name(name)
        .set_flags(flags2)
        .set_reference_sequence_id(ref_id)
        .set_alignment_start(Position::new(start2_1based).unwrap())
        .set_cigar(cigar2)
        .set_mate_reference_sequence_id(ref_id)
        .set_mate_alignment_start(Position::new(start1_1based).unwrap())
        .build();
    (r1, r2)
}

/// Build an aligned record with a name, position, flags, and an AS tag.
fn named_record_with_as(
    name: &str,
    ref_id: usize,
    start_1based: usize,
    match_len: usize,
    flags: Flags,
    as_score: i32,
) -> RecordBuf {
    let cigar: Cigar = vec![Op::new(Kind::Match, match_len)].into_iter().collect();
    let mut record = RecordBuf::builder()
        .set_name(name)
        .set_flags(flags)
        .set_reference_sequence_id(ref_id)
        .set_alignment_start(Position::new(start_1based).unwrap())
        .set_cigar(cigar)
        .build();
    record
        .data_mut()
        .insert(Tag::new(b'A', b'S'), Value::Int32(as_score));
    record
}

/// Write BAM bytes from a header and a list of records.
fn create_bam(header: &Header, records: &[RecordBuf]) -> Vec<u8> {
    let mut buf = Vec::new();
    {
        let mut writer = bam::io::Writer::new(&mut buf);
        writer.write_header(header).unwrap();
        for record in records {
            writer.write_alignment_record(header, record).unwrap();
        }
    }
    buf
}

/// Write content to a temp file with a unique name and return the path.
fn temp_file(label: &str, suffix: &str, content: &[u8]) -> PathBuf {
    let path = std::env::temp_dir().join(format!(
        "rusty_telescope_{}_{}{}",
        label,
        std::process::id(),
        suffix
    ));
    fs::write(&path, content).unwrap();
    path
}

/// Create a unique temp directory for output.
fn temp_dir(label: &str) -> PathBuf {
    let path =
        std::env::temp_dir().join(format!("rusty_telescope_{}_{}", label, std::process::id()));
    let _ = fs::remove_dir_all(&path);
    path
}

/// Read ZB (initial annotation) tag values from a BAM file on disk.
fn read_zb_tags_from_file(path: &std::path::Path) -> Vec<String> {
    let data = fs::read(path).unwrap();
    let mut reader = bam::io::Reader::new(&data[..]);
    let header = reader.read_header().unwrap();
    let zb_tag = Tag::new(b'Z', b'B');

    reader
        .record_bufs(&header)
        .map(|r| {
            let record = r.unwrap();
            match record.data().get(&zb_tag) {
                Some(Value::String(s)) => s.to_string(),
                other => panic!("Expected ZB string tag, got {:?}", other),
            }
        })
        .collect()
}

/// Read ZF (confident annotation) tag values from a BAM file on disk.
/// Returns None for records without a ZF tag.
fn read_zf_tags_from_file(path: &std::path::Path) -> Vec<Option<String>> {
    let data = fs::read(path).unwrap();
    let mut reader = bam::io::Reader::new(&data[..]);
    let header = reader.read_header().unwrap();
    let zf_tag = Tag::new(b'Z', b'F');

    reader
        .record_bufs(&header)
        .map(|r| {
            let record = r.unwrap();
            match record.data().get(&zf_tag) {
                Some(Value::String(s)) => Some(s.to_string()),
                None => None,
                other => panic!("Expected ZF string tag or absent, got {:?}", other),
            }
        })
        .collect()
}

/// Read ZR (representative) tag values from a BAM file on disk.
fn read_zr_tags_from_file(path: &std::path::Path) -> Vec<String> {
    let data = fs::read(path).unwrap();
    let mut reader = bam::io::Reader::new(&data[..]);
    let header = reader.read_header().unwrap();
    let zr_tag = Tag::new(b'Z', b'R');

    reader
        .record_bufs(&header)
        .map(|r| {
            let record = r.unwrap();
            match record.data().get(&zr_tag) {
                Some(Value::String(s)) => s.to_string(),
                other => panic!("Expected ZR string tag, got {:?}", other),
            }
        })
        .collect()
}

/// Get the output BAM path derived from the input BAM path and output dir.
fn output_bam_path(output_dir: &std::path::Path, input_bam: &std::path::Path) -> PathBuf {
    let stem = input_bam.file_stem().unwrap().to_str().unwrap();
    output_dir.join(format!("{}_annotated.bam", stem))
}

/// Get the output Jaccard TSV path derived from the input BAM path.
fn output_jaccard_path(output_dir: &std::path::Path, input_bam: &std::path::Path) -> PathBuf {
    let stem = input_bam.file_stem().unwrap().to_str().unwrap();
    output_dir.join(format!("{}_jaccard.tsv", stem))
}

/// Get the output summary TSV path derived from the input BAM path.
fn output_summary_path(output_dir: &std::path::Path, input_bam: &std::path::Path) -> PathBuf {
    let stem = input_bam.file_stem().unwrap().to_str().unwrap();
    output_dir.join(format!("{}_summary.tsv", stem))
}

// ============================================================
// Tests using --single-end (preserves per-alignment behavior)
// ============================================================

#[test]
fn test_annotation_tagging_single_end() {
    let gtf_content = "\
chr1\tsrc\texon\t100\t200\t.\t+\t.\tgene_id \"GENE_A\"
chr1\tsrc\texon\t500\t600\t.\t+\t.\tgene_id \"GENE_B\"";

    let header = build_header("chr1", 10000);
    let records = vec![
        named_record("r1", 0, 150, 20),
        named_record("r2", 0, 550, 20),
        named_record("r3", 0, 300, 20),
    ];

    let bam_path = temp_file("tagging_se", ".bam", &create_bam(&header, &records));
    let gtf_path = temp_file("tagging_se", ".gtf", gtf_content.as_bytes());
    let out_dir = temp_dir("tagging_se");

    let output = Command::new(env!("CARGO_BIN_EXE_rusty_telescope"))
        .args([
            "annotate",
            bam_path.to_str().unwrap(),
            // Default --min-overlap is 30; existing fixtures use 20bp matches.
            // Disable the filter here so legacy semantics are preserved.
            "--min-overlap",
            "0",
            "--gtf",
            gtf_path.to_str().unwrap(),
            "--output-dir",
            out_dir.to_str().unwrap(),
            "--single-end",
        ])
        .output()
        .expect("failed to run binary");

    assert!(
        output.status.success(),
        "binary failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let bam_out = output_bam_path(&out_dir, &bam_path);
    // ZB = initial annotation
    let zb_tags = read_zb_tags_from_file(&bam_out);
    assert_eq!(zb_tags, vec!["GENE_A", "GENE_B"]);
    // ZF = confident annotation (single annotation per group → always confident)
    let zf_tags = read_zf_tags_from_file(&bam_out);
    assert_eq!(
        zf_tags,
        vec![Some("GENE_A".to_string()), Some("GENE_B".to_string())]
    );
}

#[test]
fn test_include_no_feature_single_end() {
    let gtf_content = "chr1\tsrc\texon\t100\t200\t.\t+\t.\tgene_id \"GENE_A\"";
    let header = build_header("chr1", 10000);
    let records = vec![
        named_record("r1", 0, 150, 20),
        named_record("r2", 0, 500, 20),
    ];

    let bam_path = temp_file("incl_nf_se", ".bam", &create_bam(&header, &records));
    let gtf_path = temp_file("incl_nf_se", ".gtf", gtf_content.as_bytes());
    let out_dir = temp_dir("incl_nf_se");

    let output = Command::new(env!("CARGO_BIN_EXE_rusty_telescope"))
        .args([
            "annotate",
            bam_path.to_str().unwrap(),
            // Default --min-overlap is 30; existing fixtures use 20bp matches.
            // Disable the filter here so legacy semantics are preserved.
            "--min-overlap",
            "0",
            "--gtf",
            gtf_path.to_str().unwrap(),
            "--output-dir",
            out_dir.to_str().unwrap(),
            "--single-end",
            "--include-no-feature",
        ])
        .output()
        .expect("failed to run binary");

    assert!(output.status.success());
    let bam_out = output_bam_path(&out_dir, &bam_path);
    // r2's group overlaps no real GTF feature and is dropped (banding is auto-on
    // here but only relabels no-feature reads in groups that DO overlap a real
    // feature, so r2 is not resurrected). --include-no-feature only affects
    // __no_feature__ records from mixed groups.
    let zb_tags = read_zb_tags_from_file(&bam_out);
    assert_eq!(zb_tags, vec!["GENE_A"]);
    let zf_tags = read_zf_tags_from_file(&bam_out);
    assert_eq!(zf_tags, vec![Some("GENE_A".to_string())]);
}

#[test]
fn test_normalize_chr_single_end() {
    let gtf_content = "1\tsrc\texon\t100\t200\t.\t+\t.\tgene_id \"GENE_A\"";
    let header = build_header("chr1", 10000);
    let records = vec![named_record("r1", 0, 150, 20)];

    let bam_path = temp_file("norm_chr_se", ".bam", &create_bam(&header, &records));
    let gtf_path = temp_file("norm_chr_se", ".gtf", gtf_content.as_bytes());
    let out_dir = temp_dir("norm_chr_se");

    // Without --normalize-chr
    let output = Command::new(env!("CARGO_BIN_EXE_rusty_telescope"))
        .args([
            "annotate",
            bam_path.to_str().unwrap(),
            // Default --min-overlap is 30; existing fixtures use 20bp matches.
            // Disable the filter here so legacy semantics are preserved.
            "--min-overlap",
            "0",
            "--gtf",
            gtf_path.to_str().unwrap(),
            "--output-dir",
            out_dir.to_str().unwrap(),
            "--single-end",
            "--include-no-feature",
        ])
        .output()
        .expect("failed to run binary");

    assert!(output.status.success());
    let bam_out = output_bam_path(&out_dir, &bam_path);
    // The lone group overlaps no real GTF feature (chr-prefix mismatch) and is
    // dropped, so the output BAM is empty even with --include-no-feature. Banding
    // is auto-on and the read does fall in a band, but a band never makes a
    // bandless-for-real-features group eligible. The chr-prefix warning still
    // fires during overlap lookup, which runs before grouping.
    let zb_tags = read_zb_tags_from_file(&bam_out);
    assert!(zb_tags.is_empty());

    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("normalize-chr"));

    // With --normalize-chr
    let out_dir2 = temp_dir("norm_chr_se2");
    let output = Command::new(env!("CARGO_BIN_EXE_rusty_telescope"))
        .args([
            "annotate",
            bam_path.to_str().unwrap(),
            // Default --min-overlap is 30; existing fixtures use 20bp matches.
            // Disable the filter here so legacy semantics are preserved.
            "--min-overlap",
            "0",
            "--gtf",
            gtf_path.to_str().unwrap(),
            "--output-dir",
            out_dir2.to_str().unwrap(),
            "--single-end",
            "--normalize-chr",
        ])
        .output()
        .expect("failed to run binary");

    assert!(output.status.success());
    let bam_out2 = output_bam_path(&out_dir2, &bam_path);
    let zb_tags = read_zb_tags_from_file(&bam_out2);
    assert_eq!(zb_tags, vec!["GENE_A"]);
}

#[test]
fn test_overlap_summing_single_end() {
    let gtf_content = "\
chr1\tsrc\texon\t100\t200\t.\t+\t.\tgene_id \"GENE_A\"
chr1\tsrc\texon\t300\t400\t.\t+\t.\tgene_id \"GENE_A\"
chr1\tsrc\texon\t150\t350\t.\t+\t.\tgene_id \"GENE_B\"";

    let header = build_header("chr1", 10000);
    let records = vec![named_record("r1", 0, 100, 301)];

    let bam_path = temp_file("summing_se", ".bam", &create_bam(&header, &records));
    let gtf_path = temp_file("summing_se", ".gtf", gtf_content.as_bytes());
    let out_dir = temp_dir("summing_se");

    let output = Command::new(env!("CARGO_BIN_EXE_rusty_telescope"))
        .args([
            "annotate",
            bam_path.to_str().unwrap(),
            // Default --min-overlap is 30; existing fixtures use 20bp matches.
            // Disable the filter here so legacy semantics are preserved.
            "--min-overlap",
            "0",
            "--gtf",
            gtf_path.to_str().unwrap(),
            "--output-dir",
            out_dir.to_str().unwrap(),
            "--single-end",
            "--confidence",
            "0.5",
        ])
        .output()
        .expect("failed to run binary");

    assert!(output.status.success());
    let bam_out = output_bam_path(&out_dir, &bam_path);
    let zb_tags = read_zb_tags_from_file(&bam_out);
    // Single-end, single alignment overlapping both GENE_A and GENE_B.
    // This is one read group with two annotations. With default confidence=0.5,
    // softmax over equal scores (no AS tag → 0) gives 0.5 each, both pass.
    // But ZB shows the annotation of the confident representative.
    // Actually with single-end and no AS tag, this alignment is representative
    // for both GENE_A and GENE_B. Softmax of a single score (same alignment
    // for both) uses score 0 twice → 0.5 each. Both pass at threshold 0.5.
    // In default mode, the confident representative for each annotation is
    // the same alignment (idx 0). We output it once for the highest-overlap
    // confident annotation.
    // GENE_A has more overlap than GENE_B, so ZB=GENE_A, ZF=GENE_A.
    assert_eq!(zb_tags.len(), 1);
    assert_eq!(zb_tags[0], "GENE_A");
}

// ============================================================
// Paired-end tests
// ============================================================

#[test]
fn test_paired_end_grouping() {
    let gtf_content = "chr1\tsrc\texon\t100\t200\t.\t+\t.\tgene_id \"GENE_A\"";
    let header = build_header("chr1", 10000);

    // Two alignments of the same read pair, both overlap GENE_A
    let records = vec![
        named_record("read1", 0, 150, 20),
        named_record("read1", 0, 160, 20),
    ];

    let bam_path = temp_file("paired_grp", ".bam", &create_bam(&header, &records));
    let gtf_path = temp_file("paired_grp", ".gtf", gtf_content.as_bytes());
    let out_dir = temp_dir("paired_grp");

    let output = Command::new(env!("CARGO_BIN_EXE_rusty_telescope"))
        .args([
            "annotate",
            bam_path.to_str().unwrap(),
            // Default --min-overlap is 30; existing fixtures use 20bp matches.
            // Disable the filter here so legacy semantics are preserved.
            "--min-overlap",
            "0",
            "--gtf",
            gtf_path.to_str().unwrap(),
            "--output-dir",
            out_dir.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run binary");

    assert!(
        output.status.success(),
        "binary failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let bam_out = output_bam_path(&out_dir, &bam_path);
    // Only one representative should be output, with both ZB and ZF
    let zb_tags = read_zb_tags_from_file(&bam_out);
    assert_eq!(zb_tags, vec!["GENE_A"]);
    let zf_tags = read_zf_tags_from_file(&bam_out);
    assert_eq!(zf_tags, vec![Some("GENE_A".to_string())]);
}

#[test]
fn test_representative_tiebreak_as_score() {
    let gtf_content = "chr1\tsrc\texon\t100\t200\t.\t+\t.\tgene_id \"GENE_A\"";
    let header = build_header("chr1", 10000);

    // Two alignments same QNAME, same overlap, different AS scores
    let records = vec![
        named_record_with_as("read1", 0, 150, 20, Flags::SEGMENTED, 50),
        named_record_with_as("read1", 0, 150, 20, Flags::SEGMENTED, 100),
    ];

    let bam_path = temp_file("tiebreak_as", ".bam", &create_bam(&header, &records));
    let gtf_path = temp_file("tiebreak_as", ".gtf", gtf_content.as_bytes());
    let out_dir = temp_dir("tiebreak_as");

    let output = Command::new(env!("CARGO_BIN_EXE_rusty_telescope"))
        .args([
            "annotate",
            bam_path.to_str().unwrap(),
            // Default --min-overlap is 30; existing fixtures use 20bp matches.
            // Disable the filter here so legacy semantics are preserved.
            "--min-overlap",
            "0",
            "--gtf",
            gtf_path.to_str().unwrap(),
            "--output-dir",
            out_dir.to_str().unwrap(),
            "--all-alignments",
        ])
        .output()
        .expect("failed to run binary");

    assert!(output.status.success());

    let bam_out = output_bam_path(&out_dir, &bam_path);
    // With --all-alignments, both records should be present with ZB tags
    let zb_tags = read_zb_tags_from_file(&bam_out);
    assert_eq!(zb_tags.len(), 2);

    // Check ZR tags: second alignment (AS=100) should be representative
    let zr_tags = read_zr_tags_from_file(&bam_out);
    assert_eq!(zr_tags, vec!["0", "1"]);

    // ZF should only be on the confident representative (single annotation → prob 1.0)
    let zf_tags = read_zf_tags_from_file(&bam_out);
    assert_eq!(zf_tags, vec![None, Some("GENE_A".to_string())]);
}

#[test]
fn test_jaccard_output() {
    let gtf_content = "\
chr1\tsrc\texon\t100\t200\t.\t+\t.\tgene_id \"GENE_A\"
chr1\tsrc\texon\t500\t600\t.\t+\t.\tgene_id \"GENE_B\"";

    let header = build_header("chr1", 10000);
    // read1 overlaps both GENE_A and GENE_B
    // read2 overlaps only GENE_A
    let records = vec![
        named_record("read1", 0, 150, 20),
        named_record("read1", 0, 550, 20),
        named_record("read2", 0, 150, 20),
    ];

    let bam_path = temp_file("jaccard_out", ".bam", &create_bam(&header, &records));
    let gtf_path = temp_file("jaccard_out", ".gtf", gtf_content.as_bytes());
    let out_dir = temp_dir("jaccard_out");

    let output = Command::new(env!("CARGO_BIN_EXE_rusty_telescope"))
        .args([
            "annotate",
            bam_path.to_str().unwrap(),
            // Default --min-overlap is 30; existing fixtures use 20bp matches.
            // Disable the filter here so legacy semantics are preserved.
            "--min-overlap",
            "0",
            "--gtf",
            gtf_path.to_str().unwrap(),
            "--output-dir",
            out_dir.to_str().unwrap(),
            // Default --min-reads is 10; this two-read fixture would otherwise be filtered.
            "--min-reads",
            "1",
        ])
        .output()
        .expect("failed to run binary");

    assert!(output.status.success());

    let jaccard_path = output_jaccard_path(&out_dir, &bam_path);
    assert!(jaccard_path.exists(), "Jaccard TSV file should exist");

    let content = fs::read_to_string(&jaccard_path).unwrap();
    let lines: Vec<&str> = content.lines().collect();
    assert_eq!(lines[0], "annotation_a\tannotation_b\tjaccard\toverlap");
    assert!(lines.len() >= 2, "Should have at least one tuple");

    // GENE_A has {read1, read2}, GENE_B has {read1}
    // Jaccard = 1/2 = 0.5
    let parts: Vec<&str> = lines[1].split('\t').collect();
    assert_eq!(parts[0], "GENE_A");
    assert_eq!(parts[1], "GENE_B");
    let jaccard: f64 = parts[2].parse().unwrap();
    let overlap: f64 = parts[3].parse().unwrap();
    assert!((jaccard - 0.5).abs() < 1e-10);
    assert!((overlap - 1.0).abs() < 1e-10);
}

#[test]
fn test_output_dir_creation() {
    let gtf_content = "chr1\tsrc\texon\t100\t200\t.\t+\t.\tgene_id \"GENE_A\"";
    let header = build_header("chr1", 10000);
    let records = vec![named_record("r1", 0, 150, 20)];

    let bam_path = temp_file("outdir", ".bam", &create_bam(&header, &records));
    let gtf_path = temp_file("outdir", ".gtf", gtf_content.as_bytes());
    let out_dir = temp_dir("outdir_new");

    // Ensure the directory doesn't exist
    let _ = fs::remove_dir_all(&out_dir);
    assert!(!out_dir.exists());

    let output = Command::new(env!("CARGO_BIN_EXE_rusty_telescope"))
        .args([
            "annotate",
            bam_path.to_str().unwrap(),
            // Default --min-overlap is 30; existing fixtures use 20bp matches.
            // Disable the filter here so legacy semantics are preserved.
            "--min-overlap",
            "0",
            "--gtf",
            gtf_path.to_str().unwrap(),
            "--output-dir",
            out_dir.to_str().unwrap(),
            "--single-end",
        ])
        .output()
        .expect("failed to run binary");

    assert!(output.status.success());
    assert!(out_dir.exists(), "Output directory should be created");
    assert!(output_bam_path(&out_dir, &bam_path).exists());
    assert!(output_jaccard_path(&out_dir, &bam_path).exists());
    assert!(output_summary_path(&out_dir, &bam_path).exists());
}

#[test]
fn test_all_alignments_flag() {
    let gtf_content = "chr1\tsrc\texon\t100\t200\t.\t+\t.\tgene_id \"GENE_A\"";
    let header = build_header("chr1", 10000);
    let records = vec![
        named_record("read1", 0, 150, 20),
        named_record("read1", 0, 160, 20),
    ];

    let bam_path = temp_file("all_alns", ".bam", &create_bam(&header, &records));
    let gtf_path = temp_file("all_alns", ".gtf", gtf_content.as_bytes());
    let out_dir = temp_dir("all_alns");

    let output = Command::new(env!("CARGO_BIN_EXE_rusty_telescope"))
        .args([
            "annotate",
            bam_path.to_str().unwrap(),
            // Default --min-overlap is 30; existing fixtures use 20bp matches.
            // Disable the filter here so legacy semantics are preserved.
            "--min-overlap",
            "0",
            "--gtf",
            gtf_path.to_str().unwrap(),
            "--output-dir",
            out_dir.to_str().unwrap(),
            "--all-alignments",
        ])
        .output()
        .expect("failed to run binary");

    assert!(output.status.success());

    let bam_out = output_bam_path(&out_dir, &bam_path);
    let zb_tags = read_zb_tags_from_file(&bam_out);
    assert_eq!(zb_tags.len(), 2); // both alignments present

    let zr_tags = read_zr_tags_from_file(&bam_out);
    // One should be "1" (representative), one "0"
    assert!(zr_tags.contains(&"1".to_string()));
    assert!(zr_tags.contains(&"0".to_string()));

    // ZF should only be on the representative (single annotation → confident)
    let zf_tags = read_zf_tags_from_file(&bam_out);
    let confident_count = zf_tags.iter().filter(|t| t.is_some()).count();
    assert_eq!(confident_count, 1);
}

// ============================================================
// Confidence / ambiguity resolution tests
// ============================================================

#[test]
fn test_confidence_filters_ambiguous() {
    let gtf_content = "\
chr1\tsrc\texon\t100\t200\t.\t+\t.\tgene_id \"GENE_A\"
chr1\tsrc\texon\t500\t600\t.\t+\t.\tgene_id \"GENE_B\"";

    let header = build_header("chr1", 10000);
    // read1 has two alignments with similar AS scores overlapping different annotations
    let records = vec![
        named_record_with_as("read1", 0, 150, 20, Flags::SEGMENTED, 100),
        named_record_with_as("read1", 0, 550, 20, Flags::SEGMENTED, 100),
    ];

    let bam_path = temp_file("conf_ambig", ".bam", &create_bam(&header, &records));
    let gtf_path = temp_file("conf_ambig", ".gtf", gtf_content.as_bytes());
    let out_dir = temp_dir("conf_ambig");

    // With high confidence threshold, equal scores → both probabilities = 0.5, neither passes 0.9
    let output = Command::new(env!("CARGO_BIN_EXE_rusty_telescope"))
        .args([
            "annotate",
            bam_path.to_str().unwrap(),
            // Default --min-overlap is 30; existing fixtures use 20bp matches.
            // Disable the filter here so legacy semantics are preserved.
            "--min-overlap",
            "0",
            "--gtf",
            gtf_path.to_str().unwrap(),
            "--output-dir",
            out_dir.to_str().unwrap(),
            "--confidence",
            "0.9",
            "--include-no-feature",
        ])
        .output()
        .expect("failed to run binary");

    assert!(
        output.status.success(),
        "binary failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let bam_out = output_bam_path(&out_dir, &bam_path);
    let zb_tags = read_zb_tags_from_file(&bam_out);
    // Neither annotation is confident → __no_feature__
    assert_eq!(zb_tags, vec!["__no_feature__"]);
    let zf_tags = read_zf_tags_from_file(&bam_out);
    assert_eq!(zf_tags, vec![Some("__no_feature__".to_string())]);
}

#[test]
fn test_unique_mode() {
    let gtf_content = "\
chr1\tsrc\texon\t100\t200\t.\t+\t.\tgene_id \"GENE_A\"
chr1\tsrc\texon\t500\t600\t.\t+\t.\tgene_id \"GENE_B\"";

    let header = build_header("chr1", 10000);
    // read1: maps to two annotations (ambiguous, close AS scores)
    // read2: maps to only GENE_A (unique)
    let records = vec![
        named_record_with_as("read1", 0, 150, 20, Flags::SEGMENTED, 10),
        named_record_with_as("read1", 0, 550, 20, Flags::SEGMENTED, 5),
        named_record_with_as("read2", 0, 150, 20, Flags::SEGMENTED, 100),
    ];

    let bam_path = temp_file("unique_mode", ".bam", &create_bam(&header, &records));
    let gtf_path = temp_file("unique_mode", ".gtf", gtf_content.as_bytes());
    let out_dir = temp_dir("unique_mode");

    let output = Command::new(env!("CARGO_BIN_EXE_rusty_telescope"))
        .args([
            "annotate",
            bam_path.to_str().unwrap(),
            // Default --min-overlap is 30; existing fixtures use 20bp matches.
            // Disable the filter here so legacy semantics are preserved.
            "--min-overlap",
            "0",
            "--gtf",
            gtf_path.to_str().unwrap(),
            "--output-dir",
            out_dir.to_str().unwrap(),
            "--confidence",
            "1.0",
            "--max-iter",
            "1",
        ])
        .output()
        .expect("failed to run binary");

    assert!(
        output.status.success(),
        "binary failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let bam_out = output_bam_path(&out_dir, &bam_path);
    // Only read2 (unique mapping to GENE_A) should appear
    let zb_tags = read_zb_tags_from_file(&bam_out);
    assert_eq!(zb_tags, vec!["GENE_A"]);
    let zf_tags = read_zf_tags_from_file(&bam_out);
    assert_eq!(zf_tags, vec![Some("GENE_A".to_string())]);
}

#[test]
fn test_summary_tsv() {
    let gtf_content = "\
chr1\tsrc\texon\t100\t200\t.\t+\t.\tgene_id \"GENE_A\"
chr1\tsrc\texon\t500\t600\t.\t+\t.\tgene_id \"GENE_B\"";

    let header = build_header("chr1", 10000);
    // read1: maps to both GENE_A and GENE_B (ambiguous, equal AS)
    // read2: maps to GENE_A only
    // read3: maps to GENE_A only
    let records = vec![
        named_record_with_as("read1", 0, 150, 20, Flags::SEGMENTED, 100),
        named_record_with_as("read1", 0, 550, 20, Flags::SEGMENTED, 100),
        named_record("read2", 0, 150, 20),
        named_record("read3", 0, 150, 20),
    ];

    let bam_path = temp_file("summary_out", ".bam", &create_bam(&header, &records));
    let gtf_path = temp_file("summary_out", ".gtf", gtf_content.as_bytes());
    let out_dir = temp_dir("summary_out");

    let output = Command::new(env!("CARGO_BIN_EXE_rusty_telescope"))
        .args([
            "annotate",
            bam_path.to_str().unwrap(),
            // Default --min-overlap is 30; existing fixtures use 20bp matches.
            // Disable the filter here so legacy semantics are preserved.
            "--min-overlap",
            "0",
            "--gtf",
            gtf_path.to_str().unwrap(),
            "--output-dir",
            out_dir.to_str().unwrap(),
            "--confidence",
            "0.5",
        ])
        .output()
        .expect("failed to run binary");

    assert!(
        output.status.success(),
        "binary failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let summary_path = output_summary_path(&out_dir, &bam_path);
    assert!(summary_path.exists(), "Summary TSV file should exist");

    let content = fs::read_to_string(&summary_path).unwrap();
    let lines: Vec<&str> = content.lines().collect();
    assert_eq!(
        lines[0],
        "annotation\tfinal_count\tinitial_count\tunique_count\tlength\tprior_mass\tposterior_mass"
    );

    // GENE_A: initial=3 (read1+read2+read3)
    // GENE_B: initial=1 (read1)
    // read1 maps to both GENE_A and GENE_B with equal AS=100 → softmax 0.5 each
    // Single-winner resolution: tiebreak picks one (GENE_A or GENE_B) for read1
    // read2 and read3 map only to GENE_A → prob 1.0 → confident
    // So GENE_A: final=3 (or 2 if read1 went to GENE_B), GENE_B: final=0 (or 1)
    // Sorted by (final DESC, initial DESC): GENE_A first

    assert!(lines.len() >= 3, "Should have header + 2 annotation rows");

    let gene_a_parts: Vec<&str> = lines[1].split('\t').collect();
    assert_eq!(gene_a_parts[0], "GENE_A");
    assert_eq!(gene_a_parts[2], "3"); // initial_count always 3
    // GENE_A unique: read2 and read3 each have a row of size 1 → 2.
    // read1's row has size 2 (both GENE_A and GENE_B) → contributes 0 to unique.
    assert_eq!(gene_a_parts[3], "2");

    let gene_b_parts: Vec<&str> = lines[2].split('\t').collect();
    assert_eq!(gene_b_parts[0], "GENE_B");
    assert_eq!(gene_b_parts[2], "1"); // initial_count always 1
    assert_eq!(gene_b_parts[3], "0"); // never the sole row entry

    // Final counts should sum to 3 (each read pair assigned to exactly one annotation)
    let gene_a_final: usize = gene_a_parts[1].parse().unwrap();
    let gene_b_final: usize = gene_b_parts[1].parse().unwrap();
    assert_eq!(gene_a_final + gene_b_final, 3);
    // GENE_B can have at most 1 final (only read1 could go there)
    assert!(gene_b_final <= 1);

    // prior_mass / posterior_mass (cols 5, 6): unweighted Σ prob per annotation.
    // read1 splits 0.5/0.5 across {A,B}; read2 & read3 are unique on A.
    // Prior: A = 0.5 + 1 + 1 = 2.5, B = 0.5. Each column totals 3.0 (one unit
    // per non-empty fragment).
    let a_prior: f64 = gene_a_parts[5].parse().unwrap();
    let a_post: f64 = gene_a_parts[6].parse().unwrap();
    let b_prior: f64 = gene_b_parts[5].parse().unwrap();
    let b_post: f64 = gene_b_parts[6].parse().unwrap();
    assert!((a_prior - 2.5).abs() < 1e-6, "GENE_A prior_mass: {a_prior}");
    assert!((b_prior - 0.5).abs() < 1e-6, "GENE_B prior_mass: {b_prior}");
    assert!(
        (a_prior + b_prior - 3.0).abs() < 1e-6,
        "prior_mass total should equal 3 fragments",
    );
    assert!(
        (a_post + b_post - 3.0).abs() < 1e-6,
        "posterior_mass total should equal 3 fragments",
    );
}

#[test]
fn test_summary_mass_columns_max_iter_zero() {
    // With --max-iter 0 the EM is skipped, so the finalized posteriors equal the
    // priors and prior_mass must equal posterior_mass for every annotation. A
    // uniquely-mapped fragment contributes exactly 1.0 to its annotation.
    let gtf_content = "\
chr1\tsrc\texon\t100\t200\t.\t+\t.\tgene_id \"GENE_A\"
chr1\tsrc\texon\t500\t600\t.\t+\t.\tgene_id \"GENE_B\"";
    let header = build_header("chr1", 10000);
    let records = vec![
        named_record_with_as("read1", 0, 150, 20, Flags::SEGMENTED, 100),
        named_record_with_as("read1", 0, 550, 20, Flags::SEGMENTED, 100),
        named_record("read2", 0, 150, 20),
    ];
    let bam_path = temp_file("mass_mi0", ".bam", &create_bam(&header, &records));
    let gtf_path = temp_file("mass_mi0", ".gtf", gtf_content.as_bytes());
    let out_dir = temp_dir("mass_mi0");

    let output = Command::new(env!("CARGO_BIN_EXE_rusty_telescope"))
        .args([
            "annotate",
            bam_path.to_str().unwrap(),
            "--min-overlap",
            "0",
            "--gtf",
            gtf_path.to_str().unwrap(),
            "--output-dir",
            out_dir.to_str().unwrap(),
            "--max-iter",
            "0",
        ])
        .output()
        .expect("failed to run binary");
    assert!(
        output.status.success(),
        "binary failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let summary = fs::read_to_string(output_summary_path(&out_dir, &bam_path)).unwrap();
    let lines: Vec<&str> = summary.lines().collect();
    let mut mass: std::collections::HashMap<&str, (f64, f64)> = Default::default();
    for line in &lines[1..] {
        let p: Vec<&str> = line.split('\t').collect();
        mass.insert(p[0], (p[5].parse().unwrap(), p[6].parse().unwrap()));
    }
    for (ann, (prior, posterior)) in &mass {
        assert!(
            (prior - posterior).abs() < 1e-9,
            "max-iter 0: prior_mass should equal posterior_mass for {ann}: {prior} vs {posterior}",
        );
    }
    // read2 is unique on GENE_A → contributes 1.0; read1 splits 0.5 → 1.5 total.
    assert!((mass["GENE_A"].1 - 1.5).abs() < 1e-6, "{:?}", mass["GENE_A"]);
    assert!((mass["GENE_B"].1 - 0.5).abs() < 1e-6, "{:?}", mass["GENE_B"]);
}

// ============================================================
// Proper paired-end (mate pairing) tests
// ============================================================

#[test]
fn test_pe_span_merging_no_double_count() {
    // Two mates overlap the same feature; their spans overlap in the middle.
    // Feature: GENE_A at [99, 200) (GTF 100-200).
    // Read1: [99, 170) = 71bp. Read2: [130, 200) = 70bp.
    // Overlap of read1 with feature: min(170,200)-max(99,99) = 71
    // Overlap of read2 with feature: min(200,200)-max(130,99) = 70
    // Intersection of reads: [130, 170) = 40bp
    // Overlap of intersection with feature: min(170,200)-max(130,99) = 40
    // Merged overlap = 71 + 70 - 40 = 101bp (the union covers the whole feature)
    let gtf_content = "chr1\tsrc\texon\t100\t200\t.\t+\t.\tgene_id \"GENE_A\"";
    let header = build_header("chr1", 10000);

    let (r1, r2) = mate_pair("read1", 0, 100, 71, 131, 70);
    let records = vec![r1, r2];

    let bam_path = temp_file("pe_merge", ".bam", &create_bam(&header, &records));
    let gtf_path = temp_file("pe_merge", ".gtf", gtf_content.as_bytes());
    let out_dir = temp_dir("pe_merge");

    let output = Command::new(env!("CARGO_BIN_EXE_rusty_telescope"))
        .args([
            "annotate",
            bam_path.to_str().unwrap(),
            // Default --min-overlap is 30; existing fixtures use 20bp matches.
            // Disable the filter here so legacy semantics are preserved.
            "--min-overlap",
            "0",
            "--gtf",
            gtf_path.to_str().unwrap(),
            "--output-dir",
            out_dir.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run binary");

    assert!(
        output.status.success(),
        "binary failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let bam_out = output_bam_path(&out_dir, &bam_path);
    // Both BAM lines should be tagged (both mates get same tags)
    let zb_tags = read_zb_tags_from_file(&bam_out);
    assert_eq!(zb_tags.len(), 2);
    assert_eq!(zb_tags[0], "GENE_A");
    assert_eq!(zb_tags[1], "GENE_A");
}

#[test]
fn test_pe_both_mates_tagged() {
    // Verify that both mates of a pair get the same ZB and ZF tags
    let gtf_content = "chr1\tsrc\texon\t100\t200\t.\t+\t.\tgene_id \"GENE_A\"";
    let header = build_header("chr1", 10000);

    let (r1, r2) = mate_pair("read1", 0, 120, 30, 160, 30);
    let records = vec![r1, r2];

    let bam_path = temp_file("pe_both_tag", ".bam", &create_bam(&header, &records));
    let gtf_path = temp_file("pe_both_tag", ".gtf", gtf_content.as_bytes());
    let out_dir = temp_dir("pe_both_tag");

    let output = Command::new(env!("CARGO_BIN_EXE_rusty_telescope"))
        .args([
            "annotate",
            bam_path.to_str().unwrap(),
            // Default --min-overlap is 30; existing fixtures use 20bp matches.
            // Disable the filter here so legacy semantics are preserved.
            "--min-overlap",
            "0",
            "--gtf",
            gtf_path.to_str().unwrap(),
            "--output-dir",
            out_dir.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run binary");

    assert!(output.status.success());
    let bam_out = output_bam_path(&out_dir, &bam_path);

    let zb_tags = read_zb_tags_from_file(&bam_out);
    assert_eq!(zb_tags, vec!["GENE_A", "GENE_A"]);

    let zf_tags = read_zf_tags_from_file(&bam_out);
    assert_eq!(
        zf_tags,
        vec![Some("GENE_A".to_string()), Some("GENE_A".to_string())]
    );
}

#[test]
fn test_pe_as_scores_summed() {
    // Two PE alignment pairs competing for the same read group.
    // Pair 1: AS=100+80=180 overlaps GENE_A
    // Pair 2: AS=150+60=210 overlaps GENE_B
    // Pair 2 should win (higher summed AS)
    let gtf_content = "\
chr1\tsrc\texon\t100\t200\t.\t+\t.\tgene_id \"GENE_A\"
chr1\tsrc\texon\t500\t600\t.\t+\t.\tgene_id \"GENE_B\"";
    let header = build_header("chr1", 10000);

    let (r1a, r1b) = mate_pair_with_as("read1", 0, 120, 30, 160, 30, 100, 80);
    let (r2a, r2b) = mate_pair_with_as("read1", 0, 520, 30, 560, 30, 150, 60);
    let records = vec![r1a, r1b, r2a, r2b];

    let bam_path = temp_file("pe_as_sum", ".bam", &create_bam(&header, &records));
    let gtf_path = temp_file("pe_as_sum", ".gtf", gtf_content.as_bytes());
    let out_dir = temp_dir("pe_as_sum");

    let output = Command::new(env!("CARGO_BIN_EXE_rusty_telescope"))
        .args([
            "annotate",
            bam_path.to_str().unwrap(),
            // Default --min-overlap is 30; existing fixtures use 20bp matches.
            // Disable the filter here so legacy semantics are preserved.
            "--min-overlap",
            "0",
            "--gtf",
            gtf_path.to_str().unwrap(),
            "--output-dir",
            out_dir.to_str().unwrap(),
            "--confidence",
            "0.5",
        ])
        .output()
        .expect("failed to run binary");

    assert!(
        output.status.success(),
        "binary failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let bam_out = output_bam_path(&out_dir, &bam_path);
    let zf_tags = read_zf_tags_from_file(&bam_out);
    // GENE_B pair has higher summed AS (210 vs 180), so should be confident
    let confident: Vec<_> = zf_tags.iter().filter_map(|t| t.as_ref()).collect();
    assert!(
        confident.iter().all(|t| t.as_str() == "GENE_B"),
        "Expected GENE_B confident, got {:?}",
        confident
    );
}

#[test]
fn test_pe_discordant_pair_annotations() {
    // Discordant pair: read1 on chr1 (GENE_A), read2 on chr2 (GENE_C)
    // Without --exclude-discordant, both annotations should be found
    let gtf_content = "\
chr1\tsrc\texon\t100\t200\t.\t+\t.\tgene_id \"GENE_A\"
chr2\tsrc\texon\t100\t200\t.\t+\t.\tgene_id \"GENE_C\"";
    let header = build_header_2refs("chr1", 10000, "chr2", 10000);

    let (r1, r2) = discordant_mate_pair("read1", 0, 120, 50, 1, 120, 50);
    let records = vec![r1, r2];

    let bam_path = temp_file("pe_discord", ".bam", &create_bam(&header, &records));
    let gtf_path = temp_file("pe_discord", ".gtf", gtf_content.as_bytes());
    let out_dir = temp_dir("pe_discord");

    let output = Command::new(env!("CARGO_BIN_EXE_rusty_telescope"))
        .args([
            "annotate",
            bam_path.to_str().unwrap(),
            // Default --min-overlap is 30; existing fixtures use 20bp matches.
            // Disable the filter here so legacy semantics are preserved.
            "--min-overlap",
            "0",
            "--gtf",
            gtf_path.to_str().unwrap(),
            "--output-dir",
            out_dir.to_str().unwrap(),
            "--confidence",
            "0.5",
        ])
        .output()
        .expect("failed to run binary");

    assert!(
        output.status.success(),
        "binary failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    // Both mates should be in output (discordant pair still processed)
    let bam_out = output_bam_path(&out_dir, &bam_path);
    let zb_tags = read_zb_tags_from_file(&bam_out);
    assert_eq!(zb_tags.len(), 2);
}

#[test]
fn test_pe_exclude_discordant() {
    // Same as above but with --exclude-discordant
    let gtf_content = "\
chr1\tsrc\texon\t100\t200\t.\t+\t.\tgene_id \"GENE_A\"
chr2\tsrc\texon\t100\t200\t.\t+\t.\tgene_id \"GENE_C\"";
    let header = build_header_2refs("chr1", 10000, "chr2", 10000);

    // One discordant pair + one concordant pair
    let (d1, d2) = discordant_mate_pair("disc_read", 0, 120, 50, 1, 120, 50);
    let (c1, c2) = mate_pair("conc_read", 0, 120, 50, 160, 50);
    let records = vec![d1, d2, c1, c2];

    let bam_path = temp_file("pe_excl_disc", ".bam", &create_bam(&header, &records));
    let gtf_path = temp_file("pe_excl_disc", ".gtf", gtf_content.as_bytes());
    let out_dir = temp_dir("pe_excl_disc");

    let output = Command::new(env!("CARGO_BIN_EXE_rusty_telescope"))
        .args([
            "annotate",
            bam_path.to_str().unwrap(),
            // Default --min-overlap is 30; existing fixtures use 20bp matches.
            // Disable the filter here so legacy semantics are preserved.
            "--min-overlap",
            "0",
            "--gtf",
            gtf_path.to_str().unwrap(),
            "--output-dir",
            out_dir.to_str().unwrap(),
            "--exclude-discordant",
        ])
        .output()
        .expect("failed to run binary");

    assert!(
        output.status.success(),
        "binary failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    // Only the concordant pair should appear (2 BAM lines)
    let bam_out = output_bam_path(&out_dir, &bam_path);
    let zb_tags = read_zb_tags_from_file(&bam_out);
    assert_eq!(zb_tags.len(), 2);
    assert!(zb_tags.iter().all(|t| t == "GENE_A"));
}

#[test]
fn test_pe_non_overlapping_mates() {
    // Two mates that don't overlap each other but both overlap the same feature.
    // Feature GENE_A: [99, 500) (GTF 100-500)
    // Read1: [99, 149) = 50bp → overlap with feature = 50
    // Read2: [399, 449) = 50bp → overlap with feature = 50
    // No intersection between reads → merged overlap = 100bp
    let gtf_content = "chr1\tsrc\texon\t100\t500\t.\t+\t.\tgene_id \"GENE_A\"";
    let header = build_header("chr1", 10000);

    let (r1, r2) = mate_pair("read1", 0, 100, 50, 400, 50);
    let records = vec![r1, r2];

    let bam_path = temp_file("pe_nonoverlap", ".bam", &create_bam(&header, &records));
    let gtf_path = temp_file("pe_nonoverlap", ".gtf", gtf_content.as_bytes());
    let out_dir = temp_dir("pe_nonoverlap");

    let output = Command::new(env!("CARGO_BIN_EXE_rusty_telescope"))
        .args([
            "annotate",
            bam_path.to_str().unwrap(),
            // Default --min-overlap is 30; existing fixtures use 20bp matches.
            // Disable the filter here so legacy semantics are preserved.
            "--min-overlap",
            "0",
            "--gtf",
            gtf_path.to_str().unwrap(),
            "--output-dir",
            out_dir.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run binary");

    assert!(output.status.success());
    let bam_out = output_bam_path(&out_dir, &bam_path);
    let zb_tags = read_zb_tags_from_file(&bam_out);
    assert_eq!(zb_tags.len(), 2);
    assert!(zb_tags.iter().all(|t| t == "GENE_A"));
}

// ============================================================
// Stranded mode: --stranded RF/FR/F/R
// ============================================================

/// Get the output strandedness TSV path derived from the input BAM path.
fn output_strandedness_path(output_dir: &std::path::Path, input_bam: &std::path::Path) -> PathBuf {
    let stem = input_bam.file_stem().unwrap().to_str().unwrap();
    output_dir.join(format!("{}_strandedness.tsv", stem))
}

/// Run `annotate` with the given extra args on a (BAM, GTF) pair, returning
/// the output BAM's ZB tags. Asserts the binary exits successfully.
fn run_annotate_extra_args(
    bam_path: &std::path::Path,
    gtf_path: &std::path::Path,
    out_dir: &std::path::Path,
    extra: &[&str],
) -> Vec<String> {
    let mut args: Vec<&str> = vec![
        "annotate",
        bam_path.to_str().unwrap(),
        "--min-overlap",
        "0",
        "--gtf",
        gtf_path.to_str().unwrap(),
        "--output-dir",
        out_dir.to_str().unwrap(),
    ];
    args.extend_from_slice(extra);
    let output = Command::new(env!("CARGO_BIN_EXE_rusty_telescope"))
        .args(&args)
        .output()
        .expect("failed to run binary");
    assert!(
        output.status.success(),
        "annotate failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    read_zb_tags_from_file(&output_bam_path(out_dir, bam_path))
}

#[test]
fn test_stranded_se_f_keeps_sense_drops_antisense() {
    // Two genes overlapping the same coordinates, opposite strands.
    let gtf_content = "\
chr1\tsrc\texon\t100\t200\t.\t+\t.\tgene_id \"PLUS\"
chr1\tsrc\texon\t100\t200\t.\t-\t.\tgene_id \"MINUS\"";

    let header = build_header("chr1", 10000);
    // SE record on the '+' strand (no 0x10 flag) overlapping both genes.
    let records = vec![named_record_strand("r1", 0, 150, 30, false)];

    let bam_path = temp_file("stranded_f", ".bam", &create_bam(&header, &records));
    let gtf_path = temp_file("stranded_f", ".gtf", gtf_content.as_bytes());
    let out_dir = temp_dir("stranded_f");

    // F-mode: read sense to transcript → keep PLUS, drop MINUS.
    let zb = run_annotate_extra_args(
        &bam_path,
        &gtf_path,
        &out_dir,
        &["--single-end", "--stranded", "f"],
    );
    assert_eq!(zb, vec!["PLUS"]);
}

#[test]
fn test_stranded_se_r_keeps_only_antisense() {
    let gtf_content = "\
chr1\tsrc\texon\t100\t200\t.\t+\t.\tgene_id \"PLUS\"
chr1\tsrc\texon\t100\t200\t.\t-\t.\tgene_id \"MINUS\"";

    let header = build_header("chr1", 10000);
    let records = vec![named_record_strand("r1", 0, 150, 30, false)];

    let bam_path = temp_file("stranded_r", ".bam", &create_bam(&header, &records));
    let gtf_path = temp_file("stranded_r", ".gtf", gtf_content.as_bytes());
    let out_dir = temp_dir("stranded_r");

    let zb = run_annotate_extra_args(
        &bam_path,
        &gtf_path,
        &out_dir,
        &["--single-end", "--stranded", "r"],
    );
    assert_eq!(zb, vec!["MINUS"]);
}

#[test]
fn test_stranded_unknown_annotation_strand_passes() {
    // GTF strand "." — should always pass regardless of mode.
    let gtf_content = "chr1\tsrc\texon\t100\t200\t.\t.\t.\tgene_id \"UNK\"";

    let header = build_header("chr1", 10000);
    let records = vec![named_record_strand("r1", 0, 150, 30, true)];

    let bam_path = temp_file("stranded_unk", ".bam", &create_bam(&header, &records));
    let gtf_path = temp_file("stranded_unk", ".gtf", gtf_content.as_bytes());
    let out_dir = temp_dir("stranded_unk");

    let zb = run_annotate_extra_args(
        &bam_path,
        &gtf_path,
        &out_dir,
        &["--single-end", "--stranded", "f"],
    );
    assert_eq!(zb, vec!["UNK"]);
}

#[test]
fn test_stranded_pe_fr_keeps_sense() {
    let gtf_content = "\
chr1\tsrc\texon\t100\t300\t.\t+\t.\tgene_id \"PLUS\"
chr1\tsrc\texon\t100\t300\t.\t-\t.\tgene_id \"MINUS\"";

    let header = build_header("chr1", 10000);
    // FR-style mate pair: read1 forward, read2 reverse.
    let (r1, r2) = mate_pair_strand("p1", 0, 150, 30, 220, 30, false, true);
    let records = vec![r1, r2];

    let bam_path = temp_file("stranded_pe_fr", ".bam", &create_bam(&header, &records));
    let gtf_path = temp_file("stranded_pe_fr", ".gtf", gtf_content.as_bytes());
    let out_dir = temp_dir("stranded_pe_fr");

    let zb = run_annotate_extra_args(&bam_path, &gtf_path, &out_dir, &["--stranded", "fr"]);
    // Both mates of the pair receive the same ZB tag.
    assert_eq!(zb, vec!["PLUS", "PLUS"]);
}

#[test]
fn test_stranded_pe_rf_keeps_antisense() {
    let gtf_content = "\
chr1\tsrc\texon\t100\t300\t.\t+\t.\tgene_id \"PLUS\"
chr1\tsrc\texon\t100\t300\t.\t-\t.\tgene_id \"MINUS\"";

    let header = build_header("chr1", 10000);
    // FR-style flags but RF interpretation: read1 still on '+' but library is
    // dUTP-stranded so read1 is antisense — read1's mapped strand '+' → unit
    // strand '+' → RF requires opposite → keeps MINUS annotation.
    let (r1, r2) = mate_pair_strand("p1", 0, 150, 30, 220, 30, false, true);
    let records = vec![r1, r2];

    let bam_path = temp_file("stranded_pe_rf", ".bam", &create_bam(&header, &records));
    let gtf_path = temp_file("stranded_pe_rf", ".gtf", gtf_content.as_bytes());
    let out_dir = temp_dir("stranded_pe_rf");

    let zb = run_annotate_extra_args(&bam_path, &gtf_path, &out_dir, &["--stranded", "rf"]);
    assert_eq!(zb, vec!["MINUS", "MINUS"]);
}

// ============================================================
// detect-strand subcommand
// ============================================================

/// Read the strandedness TSV and return (rows, recommended_mode).
/// Each row is (annotation, strand, pct_r1_same, pct_r2_same).
fn parse_strandedness_tsv(
    path: &std::path::Path,
) -> (Vec<(String, String, String, String)>, String) {
    let content = fs::read_to_string(path).unwrap();
    let mut rows = Vec::new();
    let mut recommended = String::new();
    for (i, line) in content.lines().enumerate() {
        if i == 0 {
            assert_eq!(line, "annotation\tstrand\tpct_r1_same\tpct_r2_same");
            continue;
        }
        if let Some(rest) = line.strip_prefix("# recommended_mode: ") {
            recommended = rest.to_string();
            continue;
        }
        let parts: Vec<&str> = line.splitn(4, '\t').collect();
        assert_eq!(parts.len(), 4, "row {line:?}");
        rows.push((
            parts[0].to_string(),
            parts[1].to_string(),
            parts[2].to_string(),
            parts[3].to_string(),
        ));
    }
    (rows, recommended)
}

fn run_detect_strand(bam: &std::path::Path, gtf: &std::path::Path, out: &std::path::Path) {
    let output = Command::new(env!("CARGO_BIN_EXE_rusty_telescope"))
        .args([
            "detect-strand",
            bam.to_str().unwrap(),
            "--gtf",
            gtf.to_str().unwrap(),
            "--output-dir",
            out.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run binary");
    assert!(
        output.status.success(),
        "detect-strand failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
}

#[test]
fn test_detect_strand_recommends_f() {
    let gtf_content = "chr1\tsrc\texon\t100\t300\t.\t+\t.\tgene_id \"PLUS\"";
    let header = build_header("chr1", 10000);
    // 50 SE reads, all on '+' (matches '+' annotation).
    let records: Vec<RecordBuf> = (0..50)
        .map(|i| named_record_strand(&format!("r{i}"), 0, 150 + i, 30, false))
        .collect();

    let bam_path = temp_file("detect_f", ".bam", &create_bam(&header, &records));
    let gtf_path = temp_file("detect_f", ".gtf", gtf_content.as_bytes());
    let out_dir = temp_dir("detect_f");

    run_detect_strand(&bam_path, &gtf_path, &out_dir);
    let (rows, rec) = parse_strandedness_tsv(&output_strandedness_path(&out_dir, &bam_path));
    assert_eq!(rec, "F");
    assert_eq!(rows.len(), 1);
    let (name, strand, r1, r2) = &rows[0];
    assert_eq!(name, "PLUS");
    assert_eq!(strand, "+");
    assert_eq!(r1, "100.00");
    assert_eq!(r2, "");
}

#[test]
fn test_detect_strand_recommends_r() {
    let gtf_content = "chr1\tsrc\texon\t100\t300\t.\t+\t.\tgene_id \"PLUS\"";
    let header = build_header("chr1", 10000);
    // SE reads all on '-' (antisense to '+' annotation).
    let records: Vec<RecordBuf> = (0..50)
        .map(|i| named_record_strand(&format!("r{i}"), 0, 150 + i, 30, true))
        .collect();

    let bam_path = temp_file("detect_r", ".bam", &create_bam(&header, &records));
    let gtf_path = temp_file("detect_r", ".gtf", gtf_content.as_bytes());
    let out_dir = temp_dir("detect_r");

    run_detect_strand(&bam_path, &gtf_path, &out_dir);
    let (_, rec) = parse_strandedness_tsv(&output_strandedness_path(&out_dir, &bam_path));
    assert_eq!(rec, "R");
}

#[test]
fn test_detect_strand_recommends_fr() {
    let gtf_content = "chr1\tsrc\texon\t100\t300\t.\t+\t.\tgene_id \"PLUS\"";
    let header = build_header("chr1", 10000);
    // 30 PE pairs: read1 on '+', read2 on '-' → FR library.
    let mut records: Vec<RecordBuf> = Vec::new();
    for i in 0..30 {
        let (r1, r2) = mate_pair_strand(&format!("p{i}"), 0, 150 + i, 30, 220 + i, 30, false, true);
        records.push(r1);
        records.push(r2);
    }

    let bam_path = temp_file("detect_fr", ".bam", &create_bam(&header, &records));
    let gtf_path = temp_file("detect_fr", ".gtf", gtf_content.as_bytes());
    let out_dir = temp_dir("detect_fr");

    run_detect_strand(&bam_path, &gtf_path, &out_dir);
    let (rows, rec) = parse_strandedness_tsv(&output_strandedness_path(&out_dir, &bam_path));
    assert_eq!(rec, "FR");
    let row = &rows[0];
    assert_eq!(row.0, "PLUS");
    assert_eq!(row.1, "+");
    assert_eq!(row.2, "100.00"); // pct_r1_same
    assert_eq!(row.3, "0.00"); // pct_r2_same
}

#[test]
fn test_detect_strand_recommends_rf() {
    let gtf_content = "chr1\tsrc\texon\t100\t300\t.\t+\t.\tgene_id \"PLUS\"";
    let header = build_header("chr1", 10000);
    // 30 PE pairs: read1 on '-', read2 on '+' → RF library.
    let mut records: Vec<RecordBuf> = Vec::new();
    for i in 0..30 {
        let (r1, r2) = mate_pair_strand(&format!("p{i}"), 0, 150 + i, 30, 220 + i, 30, true, false);
        records.push(r1);
        records.push(r2);
    }

    let bam_path = temp_file("detect_rf", ".bam", &create_bam(&header, &records));
    let gtf_path = temp_file("detect_rf", ".gtf", gtf_content.as_bytes());
    let out_dir = temp_dir("detect_rf");

    run_detect_strand(&bam_path, &gtf_path, &out_dir);
    let (_, rec) = parse_strandedness_tsv(&output_strandedness_path(&out_dir, &bam_path));
    assert_eq!(rec, "RF");
}

#[test]
fn test_detect_strand_recommends_unstranded() {
    let gtf_content = "chr1\tsrc\texon\t100\t300\t.\t+\t.\tgene_id \"PLUS\"";
    let header = build_header("chr1", 10000);
    // Half reads on '+', half on '-' — no strand bias.
    let mut records: Vec<RecordBuf> = Vec::new();
    for i in 0..50 {
        records.push(named_record_strand(
            &format!("r{i}"),
            0,
            150 + i,
            30,
            i % 2 == 0,
        ));
    }

    let bam_path = temp_file("detect_uns", ".bam", &create_bam(&header, &records));
    let gtf_path = temp_file("detect_uns", ".gtf", gtf_content.as_bytes());
    let out_dir = temp_dir("detect_uns");

    run_detect_strand(&bam_path, &gtf_path, &out_dir);
    let (_, rec) = parse_strandedness_tsv(&output_strandedness_path(&out_dir, &bam_path));
    assert_eq!(rec, "unstranded");
}

// ============================================================
// --skip-jaccard / --low-memory
// ============================================================

/// Reusable fixture that produces a non-empty Jaccard TSV in default mode.
/// Mirrors `test_jaccard_output` (GENE_A and GENE_B with co-occurring read groups).
fn jaccard_fixture(label: &str) -> (PathBuf, PathBuf) {
    let gtf_content = "\
chr1\tsrc\texon\t100\t200\t.\t+\t.\tgene_id \"GENE_A\"
chr1\tsrc\texon\t500\t600\t.\t+\t.\tgene_id \"GENE_B\"";
    let header = build_header("chr1", 10000);
    let records = vec![
        named_record("read1", 0, 150, 20),
        named_record("read1", 0, 550, 20),
        named_record("read2", 0, 150, 20),
    ];
    let bam_path = temp_file(label, ".bam", &create_bam(&header, &records));
    let gtf_path = temp_file(label, ".gtf", gtf_content.as_bytes());
    (bam_path, gtf_path)
}

#[test]
fn test_skip_jaccard_no_file_produced() {
    let (bam_path, gtf_path) = jaccard_fixture("skipj_alone");
    let out_dir = temp_dir("skipj_alone");

    let output = Command::new(env!("CARGO_BIN_EXE_rusty_telescope"))
        .args([
            "annotate",
            bam_path.to_str().unwrap(),
            // Default --min-overlap is 30; existing fixtures use 20bp matches.
            // Disable the filter here so legacy semantics are preserved.
            "--min-overlap",
            "0",
            "--gtf",
            gtf_path.to_str().unwrap(),
            "--output-dir",
            out_dir.to_str().unwrap(),
            "--skip-jaccard",
        ])
        .output()
        .expect("failed to run binary");
    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(!output_jaccard_path(&out_dir, &bam_path).exists());
    // Other outputs still produced
    assert!(output_summary_path(&out_dir, &bam_path).exists());
    assert!(output_bam_path(&out_dir, &bam_path).exists());
}

#[test]
fn test_low_memory_streams_unsorted_jaccard() {
    // --low-memory writes the Jaccard TSV directly while pairs are computed
    // (no Vec, no sort). The file must exist and contain the expected pair.
    let (bam_path, gtf_path) = jaccard_fixture("lowmem_stream");
    let out_dir = temp_dir("lowmem_stream");

    let output = Command::new(env!("CARGO_BIN_EXE_rusty_telescope"))
        .args([
            "annotate",
            bam_path.to_str().unwrap(),
            // Default --min-overlap is 30; existing fixtures use 20bp matches.
            // Disable the filter here so legacy semantics are preserved.
            "--min-overlap",
            "0",
            "--gtf",
            gtf_path.to_str().unwrap(),
            "--output-dir",
            out_dir.to_str().unwrap(),
            "--low-memory",
            "--min-reads",
            "1",
        ])
        .output()
        .expect("failed to run binary");
    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let jaccard_path = output_jaccard_path(&out_dir, &bam_path);
    assert!(jaccard_path.exists(), "streamed jaccard tsv should exist");
    let content = fs::read_to_string(&jaccard_path).unwrap();
    let lines: Vec<&str> = content.lines().collect();
    assert_eq!(lines[0], "annotation_a\tannotation_b\tjaccard\toverlap");
    assert!(content.contains("GENE_A\tGENE_B"));
}

#[test]
fn test_min_reads_filters_low_coverage_annotations() {
    // Default --min-reads=10 drops annotations with fewer supporting groups.
    // The fixture has only 2 groups (read1, read2), so all pairs must be
    // filtered → file is just the header.
    let (bam_path, gtf_path) = jaccard_fixture("min_reads_filter");
    let out_dir = temp_dir("min_reads_filter");

    let output = Command::new(env!("CARGO_BIN_EXE_rusty_telescope"))
        .args([
            "annotate",
            bam_path.to_str().unwrap(),
            // Default --min-overlap is 30; existing fixtures use 20bp matches.
            // Disable the filter here so legacy semantics are preserved.
            "--min-overlap",
            "0",
            "--gtf",
            gtf_path.to_str().unwrap(),
            "--output-dir",
            out_dir.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run binary");
    assert!(output.status.success());
    let jaccard_path = output_jaccard_path(&out_dir, &bam_path);
    assert!(jaccard_path.exists());
    let content = fs::read_to_string(&jaccard_path).unwrap();
    let lines: Vec<&str> = content.lines().collect();
    assert_eq!(lines, vec!["annotation_a\tannotation_b\tjaccard\toverlap"]);
}

#[test]
fn test_low_memory_conflicts_with_all_alignments() {
    let (bam_path, gtf_path) = jaccard_fixture("lowmem_conflict");
    let out_dir = temp_dir("lowmem_conflict");

    let output = Command::new(env!("CARGO_BIN_EXE_rusty_telescope"))
        .args([
            "annotate",
            bam_path.to_str().unwrap(),
            // Default --min-overlap is 30; existing fixtures use 20bp matches.
            // Disable the filter here so legacy semantics are preserved.
            "--min-overlap",
            "0",
            "--gtf",
            gtf_path.to_str().unwrap(),
            "--output-dir",
            out_dir.to_str().unwrap(),
            "--low-memory",
            "--all-alignments",
        ])
        .output()
        .expect("failed to run binary");
    assert!(!output.status.success(), "expected non-zero exit");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("--low-memory") && stderr.contains("--all-alignments"),
        "expected conflict message in stderr, got: {stderr}"
    );
}

/// Confirm the matrix/posteriors/directive refactors don't change observable
/// BAM tags or summary content. Runs the same fixture twice — default mode
/// vs. --low-memory — and asserts ZB/ZF tags and summary.tsv match. Jaccard
/// content is not compared because --low-memory writes it unsorted.
#[test]
fn test_low_memory_equivalence_to_default_outputs() {
    let (bam_a, gtf_a) = jaccard_fixture("equiv_a");
    let (bam_b, gtf_b) = jaccard_fixture("equiv_b");
    let out_a = temp_dir("equiv_a_out");
    let out_b = temp_dir("equiv_b_out");

    let run =
        |bam: &std::path::Path, gtf: &std::path::Path, dir: &std::path::Path, extra: &[&str]| {
            let mut args: Vec<&str> = vec![
                "annotate",
                bam.to_str().unwrap(),
                "--min-overlap",
                "0",
                "--gtf",
                gtf.to_str().unwrap(),
                "--output-dir",
                dir.to_str().unwrap(),
            ];
            args.extend_from_slice(extra);
            let out = Command::new(env!("CARGO_BIN_EXE_rusty_telescope"))
                .args(&args)
                .output()
                .expect("failed to run binary");
            assert!(
                out.status.success(),
                "stderr: {}",
                String::from_utf8_lossy(&out.stderr)
            );
        };

    run(&bam_a, &gtf_a, &out_a, &[]);
    run(&bam_b, &gtf_b, &out_b, &["--low-memory"]);

    let zb_a = read_zb_tags_from_file(&output_bam_path(&out_a, &bam_a));
    let zb_b = read_zb_tags_from_file(&output_bam_path(&out_b, &bam_b));
    assert_eq!(zb_a, zb_b);

    let zf_a = read_zf_tags_from_file(&output_bam_path(&out_a, &bam_a));
    let zf_b = read_zf_tags_from_file(&output_bam_path(&out_b, &bam_b));
    assert_eq!(zf_a, zf_b);

    let sum_a = fs::read_to_string(output_summary_path(&out_a, &bam_a)).unwrap();
    let sum_b = fs::read_to_string(output_summary_path(&out_b, &bam_b)).unwrap();
    assert_eq!(sum_a, sum_b);
}

#[test]
fn test_run_summary_includes_threads_used() {
    let (bam_path, gtf_path) = jaccard_fixture("threads_log");
    let out_dir = temp_dir("threads_log");

    let output = Command::new(env!("CARGO_BIN_EXE_rusty_telescope"))
        .args([
            "annotate",
            bam_path.to_str().unwrap(),
            // Default --min-overlap is 30; existing fixtures use 20bp matches.
            // Disable the filter here so legacy semantics are preserved.
            "--min-overlap",
            "0",
            "--gtf",
            gtf_path.to_str().unwrap(),
            "--output-dir",
            out_dir.to_str().unwrap(),
            "--threads",
            "1",
        ])
        .output()
        .expect("failed to run binary");
    assert!(output.status.success());

    let stem = bam_path.file_stem().unwrap().to_str().unwrap();
    let run_summary_path = out_dir.join(format!("{stem}_run_summary.txt"));
    let content = fs::read_to_string(&run_summary_path).unwrap();
    assert!(
        content.contains("Threads used:"),
        "run_summary missing 'Threads used:': {content}"
    );
    assert!(
        content.contains("EM iterations:"),
        "run_summary missing 'EM iterations:': {content}"
    );
}

// ============================================================
// --min-overlap
// ============================================================

#[test]
fn test_min_overlap_drops_short_overlaps() {
    // Read at chr1:150-170 (1-based) overlaps GENE_A (100..200) by 21 bp.
    // With --min-overlap 30 the overlap is dropped → no real annotation
    // assignment, so default-mode output is empty.
    let gtf_content = "chr1\tsrc\texon\t100\t200\t.\t+\t.\tgene_id \"GENE_A\"";
    let header = build_header("chr1", 10000);
    let records = vec![named_record("r1", 0, 150, 21)];

    let bam_path = temp_file("min_ov_drop", ".bam", &create_bam(&header, &records));
    let gtf_path = temp_file("min_ov_drop", ".gtf", gtf_content.as_bytes());
    let out_dir = temp_dir("min_ov_drop");

    let output = Command::new(env!("CARGO_BIN_EXE_rusty_telescope"))
        .args([
            "annotate",
            bam_path.to_str().unwrap(),
            "--gtf",
            gtf_path.to_str().unwrap(),
            "--output-dir",
            out_dir.to_str().unwrap(),
            "--single-end",
            "--min-overlap",
            "30",
        ])
        .output()
        .expect("failed to run binary");
    assert!(output.status.success());
    let zb = read_zb_tags_from_file(&output_bam_path(&out_dir, &bam_path));
    assert!(
        zb.is_empty(),
        "expected the sub-threshold overlap to be filtered, got {zb:?}"
    );

    // Same fixture, --min-overlap 0 → annotation appears.
    let out_dir2 = temp_dir("min_ov_keep");
    let output = Command::new(env!("CARGO_BIN_EXE_rusty_telescope"))
        .args([
            "annotate",
            bam_path.to_str().unwrap(),
            "--gtf",
            gtf_path.to_str().unwrap(),
            "--output-dir",
            out_dir2.to_str().unwrap(),
            "--single-end",
            "--min-overlap",
            "0",
        ])
        .output()
        .expect("failed to run binary");
    assert!(output.status.success());
    let zb = read_zb_tags_from_file(&output_bam_path(&out_dir2, &bam_path));
    assert_eq!(zb, vec!["GENE_A"]);
}

#[test]
fn test_min_overlap_summed_across_exons_passes() {
    // Two short same-gene exons that individually fall below threshold but
    // together exceed it: read crossing both exons must keep GENE_A under
    // --min-overlap 30.
    let gtf_content = "\
chr1\tsrc\texon\t100\t120\t.\t+\t.\tgene_id \"GENE_A\"
chr1\tsrc\texon\t150\t170\t.\t+\t.\tgene_id \"GENE_A\"";
    let header = build_header("chr1", 10000);
    // Read [100, 170] (1-based) covers both exons.
    let records = vec![named_record("r1", 0, 100, 71)];

    let bam_path = temp_file("min_ov_summed", ".bam", &create_bam(&header, &records));
    let gtf_path = temp_file("min_ov_summed", ".gtf", gtf_content.as_bytes());
    let out_dir = temp_dir("min_ov_summed");

    let output = Command::new(env!("CARGO_BIN_EXE_rusty_telescope"))
        .args([
            "annotate",
            bam_path.to_str().unwrap(),
            "--gtf",
            gtf_path.to_str().unwrap(),
            "--output-dir",
            out_dir.to_str().unwrap(),
            "--single-end",
            "--min-overlap",
            "30",
        ])
        .output()
        .expect("failed to run binary");
    assert!(output.status.success());
    let zb = read_zb_tags_from_file(&output_bam_path(&out_dir, &bam_path));
    assert_eq!(zb, vec!["GENE_A"]);
}

// ============================================================
// summary.tsv `length` column
// ============================================================

#[test]
fn test_summary_includes_length_column() {
    // GENE_A: union of two overlapping exons → 151 bp
    // GENE_B: one exon at [300, 400] (1-based inclusive) → 101 bp
    let gtf_content = "\
chr1\tsrc\texon\t100\t200\t.\t+\t.\tgene_id \"GENE_A\"
chr1\tsrc\texon\t150\t250\t.\t+\t.\tgene_id \"GENE_A\"
chr1\tsrc\texon\t300\t400\t.\t+\t.\tgene_id \"GENE_B\"";
    let header = build_header("chr1", 10000);
    // One read each so both annotations appear in the summary.
    let records = vec![
        named_record("r1", 0, 110, 50),
        named_record("r2", 0, 320, 50),
    ];

    let bam_path = temp_file("length_col", ".bam", &create_bam(&header, &records));
    let gtf_path = temp_file("length_col", ".gtf", gtf_content.as_bytes());
    let out_dir = temp_dir("length_col");

    let output = Command::new(env!("CARGO_BIN_EXE_rusty_telescope"))
        .args([
            "annotate",
            bam_path.to_str().unwrap(),
            "--gtf",
            gtf_path.to_str().unwrap(),
            "--output-dir",
            out_dir.to_str().unwrap(),
            "--single-end",
            "--min-overlap",
            "0",
        ])
        .output()
        .expect("failed to run binary");
    assert!(
        output.status.success(),
        "binary failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let summary = fs::read_to_string(output_summary_path(&out_dir, &bam_path)).unwrap();
    let lines: Vec<&str> = summary.lines().collect();
    assert_eq!(
        lines[0],
        "annotation\tfinal_count\tinitial_count\tunique_count\tlength\tprior_mass\tposterior_mass"
    );

    let mut by_name: std::collections::HashMap<&str, &str> = Default::default();
    for line in &lines[1..] {
        let parts: Vec<&str> = line.split('\t').collect();
        by_name.insert(parts[0], parts[4]);
    }
    assert_eq!(by_name.get("GENE_A").copied(), Some("151"));
    assert_eq!(by_name.get("GENE_B").copied(), Some("101"));
}

// ============================================================
// --length-correction
// ============================================================

/// Build a fixture where length-correction tilts the EM result strongly:
/// SHORT (101 bp) and LONG (~9.9kb) overlap, with a single ambiguous
/// multi-mapping read split between them at equal AS scores.
fn length_correction_fixture(label: &str) -> (PathBuf, PathBuf) {
    let gtf_content = "\
chr1\tsrc\texon\t100\t200\t.\t+\t.\tgene_id \"SHORT\"
chr1\tsrc\texon\t100\t10000\t.\t+\t.\tgene_id \"LONG\"";
    let header = build_header("chr1", 20000);
    // Two alignments under the same QNAME: one inside the SHORT exon, one
    // outside it but inside LONG. Equal AS → equal ribbonfish prior.
    let records = vec![
        named_record_with_as("r1", 0, 130, 50, Flags::SEGMENTED, 100),
        named_record_with_as("r1", 0, 5000, 50, Flags::SEGMENTED, 100),
    ];
    let bam_path = temp_file(label, ".bam", &create_bam(&header, &records));
    let gtf_path = temp_file(label, ".gtf", gtf_content.as_bytes());
    (bam_path, gtf_path)
}

#[test]
fn test_length_correction_default_per_model() {
    // ribbonfish (default --length-correction auto = ON): the EM should
    // collapse onto SHORT given its much shorter length, so r1 gets a
    // confident ZF=SHORT.
    let (bam_path, gtf_path) = length_correction_fixture("lc_ribbon");
    let out_dir = temp_dir("lc_ribbon");
    let output = Command::new(env!("CARGO_BIN_EXE_rusty_telescope"))
        .args([
            "annotate",
            bam_path.to_str().unwrap(),
            "--gtf",
            gtf_path.to_str().unwrap(),
            "--output-dir",
            out_dir.to_str().unwrap(),
            "--single-end",
            "--min-overlap",
            "0",
        ])
        .output()
        .expect("failed to run binary");
    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let zf = read_zf_tags_from_file(&output_bam_path(&out_dir, &bam_path));
    assert!(
        zf.iter().any(|t| t.as_deref() == Some("SHORT")),
        "ribbonfish should pick SHORT under length correction; got {zf:?}",
    );
    assert!(
        zf.iter().all(|t| t.as_deref() != Some("LONG")),
        "ribbonfish should not pick LONG under length correction; got {zf:?}",
    );

    // telescope (default --length-correction auto = OFF): with global AS
    // bounds equal (single AS=100 globally), the telescope prior is uniform
    // over both annotations and stays uniform through EM → neither annotation
    // crosses the default 0.9 confidence threshold, so no ZF tag.
    let out_dir2 = temp_dir("lc_telescope");
    let output = Command::new(env!("CARGO_BIN_EXE_rusty_telescope"))
        .args([
            "annotate",
            bam_path.to_str().unwrap(),
            "--gtf",
            gtf_path.to_str().unwrap(),
            "--output-dir",
            out_dir2.to_str().unwrap(),
            "--single-end",
            "--min-overlap",
            "0",
            "--model",
            "telescope",
        ])
        .output()
        .expect("failed to run binary");
    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let zf = read_zf_tags_from_file(&output_bam_path(&out_dir2, &bam_path));
    assert!(
        zf.iter().all(|t| t.is_none()),
        "telescope without length correction should be ambiguous; got {zf:?}",
    );
}

#[test]
fn test_length_correction_explicit_off_matches_no_correction() {
    // ribbonfish + --length-correction off: same uniform-priors fixture, no
    // length skew → ambiguous, no ZF tag (as in the telescope default case).
    let (bam_path, gtf_path) = length_correction_fixture("lc_off");
    let out_dir = temp_dir("lc_off");
    let output = Command::new(env!("CARGO_BIN_EXE_rusty_telescope"))
        .args([
            "annotate",
            bam_path.to_str().unwrap(),
            "--gtf",
            gtf_path.to_str().unwrap(),
            "--output-dir",
            out_dir.to_str().unwrap(),
            "--single-end",
            "--min-overlap",
            "0",
            "--length-correction",
            "off",
        ])
        .output()
        .expect("failed to run binary");
    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let zf = read_zf_tags_from_file(&output_bam_path(&out_dir, &bam_path));
    assert!(
        zf.iter().all(|t| t.is_none()),
        "explicit --length-correction off should match no-correction behavior; got {zf:?}",
    );
}

// ============================================================
// --model telescope: π/θ split, per-fragment weights, --theta-prior
// ============================================================

/// GENE_A (1000–5000), GENE_B (100–300), GENE_C (8000–9000) are pairwise
/// disjoint, so a read's annotation is decided purely by its position. C is
/// only ever used as a unique AS anchor to fix the global AS min/max.
fn disjoint_abc_gtf() -> &'static str {
    "\
chr1\tsrc\texon\t1000\t5000\t.\t+\t.\tgene_id \"GENE_A\"
chr1\tsrc\texon\t100\t300\t.\t+\t.\tgene_id \"GENE_B\"
chr1\tsrc\texon\t8000\t9000\t.\t+\t.\tgene_id \"GENE_C\""
}

/// Parse a `_summary.tsv` into annotation -> (final, initial, unique).
fn parse_summary(
    path: &std::path::Path,
) -> std::collections::HashMap<String, (usize, usize, usize)> {
    let content = fs::read_to_string(path).unwrap();
    let mut out = std::collections::HashMap::new();
    for line in content.lines().skip(1) {
        let p: Vec<&str> = line.split('\t').collect();
        out.insert(
            p[0].to_string(),
            (
                p[1].parse().unwrap(),
                p[2].parse().unwrap(),
                p[3].parse().unwrap(),
            ),
        );
    }
    out
}

/// Run `annotate --model telescope --single-end` on a fixture and return the
/// process output plus the output dir / bam path for inspection.
fn run_telescope(
    label: &str,
    gtf: &str,
    records: &[RecordBuf],
    extra: &[&str],
) -> (std::process::Output, PathBuf, PathBuf) {
    let header = build_header("chr1", 20000);
    let bam_path = temp_file(label, ".bam", &create_bam(&header, records));
    let gtf_path = temp_file(label, ".gtf", gtf.as_bytes());
    let out_dir = temp_dir(label);
    let mut args: Vec<String> = vec![
        "annotate".into(),
        bam_path.to_str().unwrap().into(),
        "--gtf".into(),
        gtf_path.to_str().unwrap().into(),
        "--output-dir".into(),
        out_dir.to_str().unwrap().into(),
        "--single-end".into(),
        "--min-overlap".into(),
        "0".into(),
        "--model".into(),
        "telescope".into(),
    ];
    args.extend(extra.iter().map(|s| s.to_string()));
    let output = Command::new(env!("CARGO_BIN_EXE_rusty_telescope"))
        .args(&args)
        .output()
        .expect("failed to run binary");
    (output, out_dir, bam_path)
}

#[test]
fn test_telescope_theta_prior_default_and_explicit_smoke() {
    // One ambiguous {A,B} read + one unique A read + a C anchor giving a
    // non-degenerate global AS range. Default and --theta-prior 0 both work.
    let records = vec![
        named_record_with_as("amb", 0, 2000, 50, Flags::SEGMENTED, 90),
        named_record_with_as("amb", 0, 150, 50, Flags::SEGMENTED, 100),
        named_record_with_as("u1", 0, 3000, 50, Flags::SEGMENTED, 80),
        named_record_with_as("c0", 0, 8500, 50, Flags::SEGMENTED, 0),
    ];
    for extra in [vec![], vec!["--theta-prior", "0"]] {
        let (out, out_dir, bam) = run_telescope("tele_smoke", disjoint_abc_gtf(), &records, &extra);
        assert!(
            out.status.success(),
            "telescope run failed (extra={extra:?}): {}",
            String::from_utf8_lossy(&out.stderr)
        );
        let zb = read_zb_tags_from_file(&output_bam_path(&out_dir, &bam));
        assert!(!zb.is_empty(), "output BAM should carry ZB tags");
        assert!(
            zb.iter().any(|z| z == "GENE_A"),
            "GENE_A should appear in ZB tags; got {zb:?}",
        );
        assert!(output_summary_path(&out_dir, &bam).exists());
    }
}

#[test]
fn test_telescope_theta_prior_regularizes_assignment() {
    // Each theta_prior regularizing read is weighted by the maximum possible
    // fragment weight, expm1(100) ≈ 2.7e43 (treated as a perfect alignment).
    // The ambiguous read's own weight is only O(1) — a wide global AS range
    // (0..10000) keeps its scaled scores ≈1.0–1.3 (weight ≈ expm1(1.3) ≈ 2.7),
    // prior ≈0.6/0.4. So even --theta-prior 1 contributes ≈2.7e43, dwarfing
    // the O(1) ambiguous mass and fully flattening θ: the confident-pick
    // boundary now sits cleanly between the only two relevant integer values,
    // 0 and 1. One ambiguous {A,B} read, no uniques, exactly 4 EM iterations
    // (at ≥5 iterations π-only reinforcement collapses regardless):
    //   --theta-prior 0 → π and θ both compound → confident collapse onto A.
    //   --theta-prior 1 → θ flattened (≈2.7e43 ≫ weighted mass) → only π
    //                     compounds, read stays sub-0.9 → no confident pick.
    let records = vec![
        named_record_with_as("c_lo", 0, 8500, 50, Flags::SEGMENTED, 0),
        named_record_with_as("c_hi", 0, 8600, 50, Flags::SEGMENTED, 10000),
        named_record_with_as("amb", 0, 2000, 50, Flags::SEGMENTED, 130),
        named_record_with_as("amb", 0, 150, 50, Flags::SEGMENTED, 100),
    ];

    // --min-mass-delta 0 disables convergence-based early stopping so EM runs
    // the full 4 iterations this closed-form analysis assumes (with the default
    // threshold the tiny total mass would converge on iteration 1).
    let (out0, dir0, bam0) = run_telescope(
        "tele_reg0",
        disjoint_abc_gtf(),
        &records,
        &["--theta-prior", "0", "--max-iter", "4", "--min-mass-delta", "0"],
    );
    let (outr, dirr, bamr) = run_telescope(
        "tele_regH",
        disjoint_abc_gtf(),
        &records,
        &["--theta-prior", "1", "--max-iter", "4", "--min-mass-delta", "0"],
    );
    assert!(out0.status.success() && outr.status.success());

    let s0 = parse_summary(&output_summary_path(&dir0, &bam0));
    let sr = parse_summary(&output_summary_path(&dirr, &bamr));
    let a0 = s0.get("GENE_A").map(|t| t.0).unwrap_or(0);
    let ar = sr.get("GENE_A").map(|t| t.0).unwrap_or(0);
    assert!(
        a0 >= 1,
        "theta-prior 0 should confidently collapse the ambiguous read onto GENE_A \
         (got A final {a0}; summary {s0:?})",
    );
    assert!(
        ar < a0,
        "theta-prior 1 should regularize away the confident collapse: GENE_A final \
         {ar} (regularized) should be < {a0} (theta-prior 0); summaries {s0:?} vs {sr:?}",
    );
}

#[test]
fn test_telescope_unique_reads_boost_pi() {
    // Ambiguous {A,B} read with equal AS on both sides (uniform telescope
    // prior). Ten unique GENE_A reads lift π_A so the ambiguous read resolves
    // confidently to GENE_A. Without the uniques it stays ambiguous (no ZF).
    let gtf = disjoint_abc_gtf();
    let mut with_uniques = vec![
        named_record_with_as("amb", 0, 2000, 50, Flags::SEGMENTED, 100),
        named_record_with_as("amb", 0, 150, 50, Flags::SEGMENTED, 100),
    ];
    for i in 0..10 {
        with_uniques.push(named_record_with_as(
            &format!("u{i}"),
            0,
            3000,
            50,
            Flags::SEGMENTED,
            100,
        ));
    }
    let (out, dir, bam) = run_telescope("tele_pi_boost", gtf, &with_uniques, &[]);
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let zf = read_zf_tags_from_file(&output_bam_path(&dir, &bam));
    assert_eq!(
        zf[0].as_deref(),
        Some("GENE_A"),
        "ambiguous read should resolve to GENE_A once uniques boost π_A; got {zf:?}",
    );

    let only_amb = vec![
        named_record_with_as("amb", 0, 2000, 50, Flags::SEGMENTED, 100),
        named_record_with_as("amb", 0, 150, 50, Flags::SEGMENTED, 100),
    ];
    let (out2, dir2, bam2) = run_telescope("tele_pi_noboost", gtf, &only_amb, &[]);
    assert!(out2.status.success());
    let zf2 = read_zf_tags_from_file(&output_bam_path(&dir2, &bam2));
    assert!(
        zf2.iter().all(|t| t.is_none()),
        "without unique support the ambiguous read should stay unresolved; got {zf2:?}",
    );
}

#[test]
fn test_telescope_weight_uses_best_alignment() {
    // Two ambiguous reads pull θ in opposite directions. The fragment with the
    // higher best-AS (huge weight) dominates the shared θ. In fixture 1 the
    // heavy fragment prefers A; in fixture 2 it prefers B. theta-prior 0 so θ
    // dynamics are not flattened. The collapse must follow the heavy fragment.
    // --min-mass-delta 0 disables convergence-based early stopping so the
    // multi-iteration collapse runs to completion (the tiny total mass would
    // otherwise converge after one iteration, before the heavy fragment wins).
    let gtf = disjoint_abc_gtf();

    // Fixture 1: P heavy & prefers A (AS A=100,B=99); Q light & prefers B
    // (AS A=11,B=12). C anchor fixes global min AS=0.
    let f1 = vec![
        named_record_with_as("c0", 0, 8500, 50, Flags::SEGMENTED, 0),
        named_record_with_as("P", 0, 2000, 50, Flags::SEGMENTED, 100),
        named_record_with_as("P", 0, 150, 50, Flags::SEGMENTED, 99),
        named_record_with_as("Q", 0, 2000, 50, Flags::SEGMENTED, 11),
        named_record_with_as("Q", 0, 150, 50, Flags::SEGMENTED, 12),
    ];
    let (o1, d1, b1) = run_telescope(
        "tele_w_a",
        gtf,
        &f1,
        &["--theta-prior", "0", "--min-mass-delta", "0"],
    );
    assert!(
        o1.status.success(),
        "{}",
        String::from_utf8_lossy(&o1.stderr)
    );
    let s1 = parse_summary(&output_summary_path(&d1, &b1));
    let a1 = s1.get("GENE_A").map(|t| t.0).unwrap_or(0);
    let b1c = s1.get("GENE_B").map(|t| t.0).unwrap_or(0);
    assert!(
        a1 > b1c,
        "heavy fragment prefers A ⇒ GENE_A should win (A={a1}, B={b1c})",
    );

    // Fixture 2: P heavy & prefers B (AS A=99,B=100); Q light & prefers A.
    let f2 = vec![
        named_record_with_as("c0", 0, 8500, 50, Flags::SEGMENTED, 0),
        named_record_with_as("P", 0, 2000, 50, Flags::SEGMENTED, 99),
        named_record_with_as("P", 0, 150, 50, Flags::SEGMENTED, 100),
        named_record_with_as("Q", 0, 2000, 50, Flags::SEGMENTED, 12),
        named_record_with_as("Q", 0, 150, 50, Flags::SEGMENTED, 11),
    ];
    let (o2, d2, b2) = run_telescope(
        "tele_w_b",
        gtf,
        &f2,
        &["--theta-prior", "0", "--min-mass-delta", "0"],
    );
    assert!(
        o2.status.success(),
        "{}",
        String::from_utf8_lossy(&o2.stderr)
    );
    let s2 = parse_summary(&output_summary_path(&d2, &b2));
    let a2 = s2.get("GENE_A").map(|t| t.0).unwrap_or(0);
    let b2c = s2.get("GENE_B").map(|t| t.0).unwrap_or(0);
    assert!(
        b2c > a2,
        "heavy fragment prefers B ⇒ GENE_B should win (A={a2}, B={b2c})",
    );
}

#[test]
fn test_telescope_max_iter_zero_equals_priors() {
    // Ambiguous {A,B}: A-side AS=100, B-side AS=0 ⇒ telescope prior ≈ [1,0].
    // With --max-iter 0 telescope_em returns priors verbatim, so the
    // confidence pick is GENE_A purely from the prior.
    let records = vec![
        named_record_with_as("amb", 0, 2000, 50, Flags::SEGMENTED, 100),
        named_record_with_as("amb", 0, 150, 50, Flags::SEGMENTED, 0),
    ];
    let (out, dir, bam) = run_telescope(
        "tele_iter0",
        disjoint_abc_gtf(),
        &records,
        &["--max-iter", "0"],
    );
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let zf = read_zf_tags_from_file(&output_bam_path(&dir, &bam));
    assert!(
        zf.iter().any(|t| t.as_deref() == Some("GENE_A")),
        "max-iter 0 should confidently pick GENE_A from the telescope prior; got {zf:?}",
    );
    assert!(
        zf.iter().all(|t| t.as_deref() != Some("GENE_B")),
        "GENE_B has prior ~0 and must never be the confident pick; got {zf:?}",
    );
}

#[test]
fn test_ribbonfish_unaffected_by_theta_prior() {
    // --theta-prior is telescope-only: passing it under the default ribbonfish
    // model must produce byte-identical ZB/ZF and summary output.
    let run = |label: &str, extra: &[&str]| {
        let (bam_path, gtf_path) = length_correction_fixture(label);
        let out_dir = temp_dir(label);
        let mut args = vec![
            "annotate",
            bam_path.to_str().unwrap(),
            "--gtf",
            gtf_path.to_str().unwrap(),
            "--output-dir",
            out_dir.to_str().unwrap(),
            "--single-end",
            "--min-overlap",
            "0",
        ];
        args.extend_from_slice(extra);
        let output = Command::new(env!("CARGO_BIN_EXE_rusty_telescope"))
            .args(&args)
            .output()
            .expect("failed to run binary");
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        let bam_out = output_bam_path(&out_dir, &bam_path);
        let zb = read_zb_tags_from_file(&bam_out);
        let zf = read_zf_tags_from_file(&bam_out);
        let summary = fs::read_to_string(output_summary_path(&out_dir, &bam_path)).unwrap();
        (zb, zf, summary)
    };
    let baseline = run("rb_tp_base", &[]);
    let with_tp = run("rb_tp_set", &["--theta-prior", "5"]);
    assert_eq!(
        baseline.0, with_tp.0,
        "ZB must be unaffected by --theta-prior"
    );
    assert_eq!(
        baseline.1, with_tp.1,
        "ZF must be unaffected by --theta-prior"
    );
    assert_eq!(
        baseline.2, with_tp.2,
        "summary must be unaffected by --theta-prior"
    );
}

#[test]
fn test_ribbonfish_length_correction_default_still_short() {
    // Guards the main.rs model-dispatch refactor: ribbonfish (default,
    // length-correction auto = ON) still collapses the ambiguous read onto the
    // much shorter SHORT annotation.
    let (bam_path, gtf_path) = length_correction_fixture("rb_lc_guard");
    let out_dir = temp_dir("rb_lc_guard");
    let output = Command::new(env!("CARGO_BIN_EXE_rusty_telescope"))
        .args([
            "annotate",
            bam_path.to_str().unwrap(),
            "--gtf",
            gtf_path.to_str().unwrap(),
            "--output-dir",
            out_dir.to_str().unwrap(),
            "--single-end",
            "--min-overlap",
            "0",
        ])
        .output()
        .expect("failed to run binary");
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let zf = read_zf_tags_from_file(&output_bam_path(&out_dir, &bam_path));
    assert!(
        zf.iter().any(|t| t.as_deref() == Some("SHORT")),
        "ribbonfish should still pick SHORT under length correction; got {zf:?}",
    );
    assert!(
        zf.iter().all(|t| t.as_deref() != Some("LONG")),
        "ribbonfish should not pick LONG; got {zf:?}",
    );
}

#[test]
fn test_telescope_length_correction_on_applies() {
    // --model telescope + explicit --length-correction on: the run succeeds
    // with no "ignored" warning, and the E-step's θ_a/L_a correction collapses
    // the ambiguous read onto the much shorter SHORT annotation — unlike the
    // telescope default (correction off), which leaves it ambiguous (no ZF).
    let (bam_path, gtf_path) = length_correction_fixture("tele_lc_on");
    let out_dir = temp_dir("tele_lc_on");
    let output = Command::new(env!("CARGO_BIN_EXE_rusty_telescope"))
        .args([
            "annotate",
            bam_path.to_str().unwrap(),
            "--gtf",
            gtf_path.to_str().unwrap(),
            "--output-dir",
            out_dir.to_str().unwrap(),
            "--single-end",
            "--min-overlap",
            "0",
            "--model",
            "telescope",
            "--length-correction",
            "on",
        ])
        .output()
        .expect("failed to run binary");
    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        !stderr.contains("length correction is ignored"),
        "telescope now supports length correction; no warning expected:\n{stderr}",
    );
    let zf = read_zf_tags_from_file(&output_bam_path(&out_dir, &bam_path));
    assert!(
        zf.iter().any(|t| t.as_deref() == Some("SHORT")),
        "telescope with length correction should pick SHORT; got {zf:?}",
    );
    assert!(
        zf.iter().all(|t| t.as_deref() != Some("LONG")),
        "telescope with length correction should not pick LONG; got {zf:?}",
    );
}

// ============================================================
// --band-no-feature: partition __no_feature__ reads into cytobands
// ============================================================

/// Two bands on chr1: p1=[0,1000), q1=[1000,2000).
const TINY_CYTOBAND: &str = "chr1\t0\t1000\tp1\tgneg\nchr1\t1000\t2000\tq1\tgneg\n";

#[test]
fn test_band_no_feature_splits_summary_and_aggregates_run_summary() {
    let gtf_content = "chr1\tsrc\texon\t100\t200\t.\t+\t.\tgene_id \"GENE_A\"";
    let header = build_header("chr1", 10000);
    // Banding only relabels no-feature reads in groups that ALSO overlap a real
    // feature. So we use mixed multimappers: each of rm1/rm2 has one alignment on
    // GENE_A (AS 10) and one intergenic alignment in band q1 (AS 100). With
    // --max-iter 0 the high-AS band is the deterministic confident pick. g1 is a
    // pure-GENE_A read to show the real-feature path still resolves normally.
    let records = vec![
        named_record_with_as("g1", 0, 150, 20, Flags::SEGMENTED, 50),
        named_record_with_as("rm1", 0, 150, 20, Flags::SEGMENTED, 10),
        named_record_with_as("rm1", 0, 1500, 20, Flags::SEGMENTED, 100),
        named_record_with_as("rm2", 0, 150, 20, Flags::SEGMENTED, 10),
        named_record_with_as("rm2", 0, 1500, 20, Flags::SEGMENTED, 100),
    ];

    let bam_path = temp_file("band_split", ".bam", &create_bam(&header, &records));
    let gtf_path = temp_file("band_split", ".gtf", gtf_content.as_bytes());
    let cyto_path = temp_file("band_split", ".cytoband.txt", TINY_CYTOBAND.as_bytes());
    let out_dir = temp_dir("band_split");

    let output = Command::new(env!("CARGO_BIN_EXE_rusty_telescope"))
        .args([
            "annotate",
            bam_path.to_str().unwrap(),
            "--min-overlap",
            "0",
            "--gtf",
            gtf_path.to_str().unwrap(),
            "--cytoband",
            cyto_path.to_str().unwrap(),
            "--band-no-feature",
            "on",
            "--max-iter",
            "0",
            "--output-dir",
            out_dir.to_str().unwrap(),
            "--single-end",
            "--include-no-feature",
        ])
        .output()
        .expect("failed to run binary");
    assert!(
        output.status.success(),
        "binary failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    // summary.tsv splits by band: q1 gets its own row reporting its real length
    // (1000), not the arbitrary __no_feature__ 1_000_000.
    let summary = fs::read_to_string(output_summary_path(&out_dir, &bam_path)).unwrap();
    let band_row = summary
        .lines()
        .find(|l| l.starts_with("__no_feature_1q1__\t"))
        .expect("summary should contain a __no_feature_1q1__ row");
    let cols: Vec<&str> = band_row.split('\t').collect();
    assert_eq!(cols[4], "1000", "band row should report the band length");
    assert_eq!(cols[1], "2", "rm1 and rm2 confidently assigned to band q1");
    // The real-feature path is unaffected: GENE_A still resolves (g1).
    let gene_row = summary
        .lines()
        .find(|l| l.starts_with("GENE_A\t"))
        .expect("summary should contain a GENE_A row");
    assert_eq!(gene_row.split('\t').collect::<Vec<_>>()[1], "1");
    // The generic single __no_feature__ bucket must not appear in band mode.
    assert!(
        !summary.lines().any(|l| l.starts_with("__no_feature__\t")),
        "band mode should not emit a generic __no_feature__ row:\n{summary}"
    );

    // run_summary.txt aggregates band picks into one "No feature" count.
    let run_summary_path = out_dir.join(format!(
        "{}_run_summary.txt",
        bam_path.file_stem().unwrap().to_str().unwrap()
    ));
    let run_summary = fs::read_to_string(run_summary_path).unwrap();
    let nf_line = run_summary
        .lines()
        .find(|l| l.starts_with("No feature:"))
        .expect("run summary should have a No feature line");
    let nf_count: usize = nf_line.split_whitespace().last().unwrap().parse().unwrap();
    assert_eq!(nf_count, 2, "No feature should aggregate both band reads");

    // With --include-no-feature, band reads are written with the band name, and
    // the confident GENE_A read (g1) is written too.
    let mut zb = read_zb_tags_from_file(&output_bam_path(&out_dir, &bam_path));
    zb.sort();
    assert_eq!(
        zb,
        vec!["GENE_A", "__no_feature_1q1__", "__no_feature_1q1__"]
    );
}

#[test]
fn test_band_no_feature_excluded_from_bam_by_default() {
    let gtf_content = "chr1\tsrc\texon\t100\t200\t.\t+\t.\tgene_id \"GENE_A\"";
    let header = build_header("chr1", 10000);
    // Mixed read whose intergenic alignment (band q1, AS 100) beats its GENE_A
    // alignment (AS 10) under --max-iter 0, so the confident pick is the band.
    let records = vec![
        named_record_with_as("rm", 0, 150, 20, Flags::SEGMENTED, 10),
        named_record_with_as("rm", 0, 1500, 20, Flags::SEGMENTED, 100),
    ];

    let bam_path = temp_file("band_excl", ".bam", &create_bam(&header, &records));
    let gtf_path = temp_file("band_excl", ".gtf", gtf_content.as_bytes());
    let cyto_path = temp_file("band_excl", ".cytoband.txt", TINY_CYTOBAND.as_bytes());
    let out_dir = temp_dir("band_excl");

    let output = Command::new(env!("CARGO_BIN_EXE_rusty_telescope"))
        .args([
            "annotate",
            bam_path.to_str().unwrap(),
            "--min-overlap",
            "0",
            "--gtf",
            gtf_path.to_str().unwrap(),
            "--cytoband",
            cyto_path.to_str().unwrap(),
            "--band-no-feature",
            "on",
            "--max-iter",
            "0",
            "--output-dir",
            out_dir.to_str().unwrap(),
            "--single-end",
        ])
        .output()
        .expect("failed to run binary");
    assert!(output.status.success());

    // Without --include-no-feature, the band-assigned confident read is gated out
    // of the BAM (same as __no_feature__), but still appears in the summary.
    let zb = read_zb_tags_from_file(&output_bam_path(&out_dir, &bam_path));
    assert!(
        zb.is_empty(),
        "band read should be excluded from BAM: {zb:?}"
    );
    let summary = fs::read_to_string(output_summary_path(&out_dir, &bam_path)).unwrap();
    assert!(
        summary
            .lines()
            .any(|l| l.starts_with("__no_feature_1q1__\t")),
        "band still appears in summary regardless of --include-no-feature"
    );
}

#[test]
fn test_band_no_feature_jaccard_at_band_level() {
    let gtf_content = "chr1\tsrc\texon\t100\t200\t.\t+\t.\tgene_id \"GENE_A\"";
    let header = build_header("chr1", 10000);
    // One multimapping read: one alignment on GENE_A, one intergenic in band q1.
    // Single-end groups them by QNAME, so GENE_A and the band co-occur.
    let records = vec![
        named_record("rm", 0, 150, 20),
        named_record("rm", 0, 1500, 20),
    ];

    let bam_path = temp_file("band_jac", ".bam", &create_bam(&header, &records));
    let gtf_path = temp_file("band_jac", ".gtf", gtf_content.as_bytes());
    let cyto_path = temp_file("band_jac", ".cytoband.txt", TINY_CYTOBAND.as_bytes());
    let out_dir = temp_dir("band_jac");

    let output = Command::new(env!("CARGO_BIN_EXE_rusty_telescope"))
        .args([
            "annotate",
            bam_path.to_str().unwrap(),
            "--min-overlap",
            "0",
            "--gtf",
            gtf_path.to_str().unwrap(),
            "--cytoband",
            cyto_path.to_str().unwrap(),
            "--band-no-feature",
            "on",
            "--output-dir",
            out_dir.to_str().unwrap(),
            "--single-end",
            "--min-reads",
            "0",
            "--similarity-threshold",
            "0",
        ])
        .output()
        .expect("failed to run binary");
    assert!(output.status.success());

    let jaccard = fs::read_to_string(output_jaccard_path(&out_dir, &bam_path)).unwrap();
    assert!(
        jaccard.contains("__no_feature_1q1__"),
        "jaccard should reference bands at band level:\n{jaccard}"
    );
}

#[test]
fn test_band_no_feature_preserves_drop_count() {
    // The invariant: which read groups are processed vs dropped is decided purely
    // by real-GTF-feature overlap, so banding must not change the retained/dropped
    // counts. g1 = pure gene; m1 = mixed (gene + intergenic); i1/i2 = all
    // intergenic (dropped either way).
    let gtf_content = "chr1\tsrc\texon\t100\t200\t.\t+\t.\tgene_id \"GENE_A\"";
    let header = build_header("chr1", 10000);
    let records = vec![
        named_record("g1", 0, 150, 20),
        named_record("m1", 0, 150, 20),
        named_record("m1", 0, 1500, 20),
        named_record("i1", 0, 1500, 20),
        named_record("i2", 0, 1600, 20),
    ];

    let bam_path = temp_file("band_inv", ".bam", &create_bam(&header, &records));
    let gtf_path = temp_file("band_inv", ".gtf", gtf_content.as_bytes());
    let cyto_path = temp_file("band_inv", ".cytoband.txt", TINY_CYTOBAND.as_bytes());

    let run = |band_mode: &str, label: &str| -> (usize, usize) {
        let out_dir = temp_dir(label);
        let output = Command::new(env!("CARGO_BIN_EXE_rusty_telescope"))
            .args([
                "annotate",
                bam_path.to_str().unwrap(),
                "--min-overlap",
                "0",
                "--gtf",
                gtf_path.to_str().unwrap(),
                "--cytoband",
                cyto_path.to_str().unwrap(),
                "--band-no-feature",
                band_mode,
                "--output-dir",
                out_dir.to_str().unwrap(),
                "--single-end",
            ])
            .output()
            .expect("failed to run binary");
        assert!(output.status.success());
        let rs = fs::read_to_string(out_dir.join(format!(
            "{}_run_summary.txt",
            bam_path.file_stem().unwrap().to_str().unwrap()
        )))
        .unwrap();
        let parse = |prefix: &str| -> usize {
            rs.lines()
                .find(|l| l.starts_with(prefix))
                .and_then(|l| l.split_whitespace().last())
                .and_then(|s| s.parse().ok())
                .unwrap_or_else(|| panic!("missing '{prefix}' in run summary:\n{rs}"))
        };
        (parse("Read groups:"), parse("Dropped"))
    };

    let on = run("on", "band_inv_on");
    let off = run("off", "band_inv_off");
    assert_eq!(
        on, off,
        "band mode must not change (read_groups, dropped) counts"
    );
    // 2 retained (g1, m1), 2 dropped (i1, i2) regardless of band mode.
    assert_eq!(on, (2, 2));
}

#[test]
fn test_unmapped_reads_dropped_and_counted() {
    let gtf_content = "chr1\tsrc\texon\t100\t200\t.\t+\t.\tgene_id \"GENE_A\"";
    let header = build_header("chr1", 10000);
    // One mapped read overlapping GENE_A, plus three unmapped reads: two with
    // RNAME '*' and one "placed unmapped" (0x4 flag, CIGAR '*', RNAME set). All
    // three are dropped by the 0x4 unmapped flag.
    let records = vec![
        named_record("r1", 0, 150, 20),
        unmapped_record("u1"),
        unmapped_record("u2"),
        placed_unmapped_record("u3", 0, 150),
    ];

    let bam_path = temp_file("unmapped", ".bam", &create_bam(&header, &records));
    let gtf_path = temp_file("unmapped", ".gtf", gtf_content.as_bytes());
    let out_dir = temp_dir("unmapped");

    let output = Command::new(env!("CARGO_BIN_EXE_rusty_telescope"))
        .args([
            "annotate",
            bam_path.to_str().unwrap(),
            "--min-overlap",
            "0",
            "--gtf",
            gtf_path.to_str().unwrap(),
            "--output-dir",
            out_dir.to_str().unwrap(),
            "--single-end",
        ])
        .output()
        .expect("failed to run binary");
    assert!(
        output.status.success(),
        "binary failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    // Unmapped reads are dropped and counted in the run summary.
    let rs = fs::read_to_string(out_dir.join(format!(
        "{}_run_summary.txt",
        bam_path.file_stem().unwrap().to_str().unwrap()
    )))
    .unwrap();
    let line = rs
        .lines()
        .find(|l| l.starts_with("Dropped (unmapped):"))
        .expect("run summary should report dropped unmapped reads");
    let count: usize = line.split_whitespace().last().unwrap().parse().unwrap();
    assert_eq!(count, 3);

    // Only the mapped read is tagged and written; the unmapped reads are gone.
    let zb = read_zb_tags_from_file(&output_bam_path(&out_dir, &bam_path));
    assert_eq!(zb, vec!["GENE_A"]);
}
