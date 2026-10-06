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
