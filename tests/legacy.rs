//! Legacy (`assign`) mode against Telescope's own bundled test data.
//!
//! `telescope/data/telescope_report.tsv` is the report the Python reference
//! ships for `alignment.bam` + `annotation.gtf`; apart from the version in the
//! header line, `assign` must reproduce it byte for byte.

use std::fs;
use std::process::Command;

#[test]
fn assign_reproduces_bundled_telescope_report() {
    let root = env!("CARGO_MANIFEST_DIR");
    let out = std::env::temp_dir().join(format!("rusty_telescope_legacy_{}", std::process::id()));
    fs::create_dir_all(&out).unwrap();

    let status = Command::new(env!("CARGO_BIN_EXE_rusty_telescope"))
        .args(["assign", "--quiet", "--exp_tag", "bundled", "--outdir"])
        .arg(&out)
        .arg(format!("{root}/telescope/data/alignment.bam"))
        .arg(format!("{root}/telescope/data/annotation.gtf"))
        .status()
        .unwrap();
    assert!(status.success());

    let got = fs::read_to_string(out.join("bundled-telescope_report.tsv")).unwrap();
    let want = fs::read_to_string(format!("{root}/telescope/data/telescope_report.tsv")).unwrap();
    fs::remove_dir_all(&out).ok();

    let (got_head, got_body) = got.split_once('\n').unwrap();
    let (want_head, want_body) = want.split_once('\n').unwrap();
    assert_eq!(got_body, want_body);
    // Same run counters; only the version field differs (bundled file is 1.0.2).
    let counters = |h: &str| h.split('\t').skip(2).map(str::to_string).collect::<Vec<_>>();
    assert_eq!(counters(got_head), counters(want_head));
}

/// (flag, has XP tag) for every record of a BAM.
fn bam_records(path: &std::path::Path) -> Vec<(u16, bool)> {
    use noodles::bam;
    use noodles::sam::alignment::{RecordBuf, record::data::field::Tag};
    let mut reader = bam::io::reader::Builder.build_from_path(path).unwrap();
    let header = reader.read_header().unwrap();
    let mut rec = RecordBuf::default();
    let mut out = Vec::new();
    while reader.read_record_buf(&header, &mut rec).unwrap() != 0 {
        out.push((u16::from(rec.flags()), rec.data().get(&Tag::new(b'X', b'P')).is_some()));
    }
    out
}

fn run_updated_sam(label: &str, extra: &[&str]) -> std::path::PathBuf {
    let root = env!("CARGO_MANIFEST_DIR");
    let out = std::env::temp_dir().join(format!("rusty_telescope_{label}_{}", std::process::id()));
    fs::create_dir_all(&out).unwrap();
    let status = Command::new(env!("CARGO_BIN_EXE_rusty_telescope"))
        .args(["assign", "--quiet", "--updated_sam", "--exp_tag", "t", "--outdir"])
        .arg(&out)
        .args(extra)
        .arg(format!("{root}/telescope/data/alignment.bam"))
        .arg(format!("{root}/telescope/data/annotation.gtf"))
        .status()
        .unwrap();
    assert!(status.success());
    out
}

#[test]
fn updated_sam_legacy_keeps_every_alignment() {
    let out = run_updated_sam("usam_legacy", &["--legacy"]);
    let updated = bam_records(&out.join("t-updated.bam"));
    let other = bam_records(&out.join("t-other.bam"));
    assert!(out.join("t-tmp_tele.bam").exists());
    fs::remove_dir_all(&out).ok();
    // Every fragment in the bundled data overlaps the annotation, so all
    // 66,414 input records land in updated.bam (as with Telescope) and
    // other.bam is empty. Most are marked secondary.
    assert_eq!(updated.len(), 66414);
    assert!(other.is_empty());
    assert!(updated.iter().any(|&(flag, _)| flag & 0x100 != 0));
}

#[test]
fn updated_sam_default_keeps_only_assigned_alignments() {
    let out = run_updated_sam("usam_assigned", &[]);
    let updated = bam_records(&out.join("t-updated.bam"));
    let leftovers = out.join("t-other.bam").exists() || out.join("t-tmp_tele.bam").exists();
    fs::remove_dir_all(&out).ok();
    assert!(!leftovers);
    // One assigned alignment (two mates) per fragment at most: 1000 fragments.
    assert!(!updated.is_empty() && updated.len() <= 2000 && updated.len().is_multiple_of(2));
    assert!(updated.iter().all(|&(flag, has_xp)| flag & 0x100 == 0 && has_xp));
}
