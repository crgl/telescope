//! `assign` against Python Telescope 1.0.4.1 (hanalysis fork) on small real-data fixtures.
//!
//! `tests/data/legacy/cases.tsv` lists the cases; `expected/` holds the reports
//! (and, for some, the updated-BAM contents) the reference produced. Every case
//! must match byte for byte. See `tests/data/legacy/README.md` for what each
//! fixture is there to exercise and how the expected files were made.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

use noodles::bam;
use noodles::sam::alignment::{
    RecordBuf, record::data::field::Tag, record_buf::data::field::Value,
};

struct Case {
    name: String,
    alignment: String,
    gtf: String,
    args: Vec<String>,
    check_updated: bool,
}

fn data_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/legacy")
}

fn cases() -> Vec<Case> {
    fs::read_to_string(data_dir().join("cases.tsv"))
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            assert_eq!(f.len(), 6, "malformed cases.tsv line: {l}");
            Case {
                name: f[0].to_string(),
                alignment: f[2].to_string(),
                gtf: f[3].to_string(),
                args: f[4].split_whitespace().filter(|a| *a != "-").map(str::to_string).collect(),
                check_updated: f[5] == "yes",
            }
        })
        .collect()
}

/// The fields `--updated_sam` is responsible for, one line per record, in the
/// layout `make_expected.sh` writes.
fn updated_table(path: &Path) -> String {
    let mut reader = bam::io::reader::Builder.build_from_path(path).unwrap();
    let header = reader.read_header().unwrap();
    let names: Vec<String> =
        header.reference_sequences().keys().map(|k| String::from_utf8_lossy(k.as_ref()).into_owned()).collect();
    let tag = |rec: &RecordBuf, t: &[u8; 2]| -> String {
        match rec.data().get(&Tag::new(t[0], t[1])) {
            None => ".".to_string(),
            Some(Value::String(s)) => String::from_utf8_lossy(s.as_ref()).into_owned(),
            Some(v) => v.as_int().map_or_else(|| "?".to_string(), |i| i.to_string()),
        }
    };
    let mut rec = RecordBuf::default();
    let mut out = String::new();
    while reader.read_record_buf(&header, &mut rec).unwrap() != 0 {
        let qname = rec.name().map_or_else(|| "*".to_string(), |n| n.to_string());
        let rname = rec.reference_sequence_id().map_or("*", |i| names[i].as_str());
        let pos = rec.alignment_start().map_or(0, usize::from);
        let mapq = rec.mapping_quality().map_or(255, u8::from);
        out.push_str(&format!(
            "{qname}\t{}\t{rname}\t{pos}\t{mapq}\t{}\t{}\t{}\t{}\t{}\n",
            u16::from(rec.flags()),
            tag(&rec, b"ZF"),
            tag(&rec, b"ZT"),
            tag(&rec, b"ZB"),
            tag(&rec, b"XP"),
            tag(&rec, b"YC"),
        ));
    }
    out
}

/// First line where two texts differ, for a readable failure.
fn first_difference(got: &str, want: &str) -> String {
    for (i, (g, w)) in got.lines().zip(want.lines()).enumerate() {
        if g != w {
            return format!("line {}:\n  got:  {g}\n  want: {w}", i + 1);
        }
    }
    format!("line counts differ: got {}, want {}", got.lines().count(), want.lines().count())
}

#[test]
fn assign_matches_reference_on_every_fixture() {
    let dir = data_dir();
    let scratch = std::env::temp_dir().join(format!("rusty_telescope_fixtures_{}", std::process::id()));
    let mut failures = Vec::new();
    let all = cases();
    assert!(all.len() >= 20, "expected the full case list, found {}", all.len());

    for case in &all {
        let out = scratch.join(&case.name);
        fs::create_dir_all(&out).unwrap();
        let mut cmd = Command::new(env!("CARGO_BIN_EXE_rusty_telescope"));
        cmd.args(["assign", "--quiet", "--exp_tag", "t", "--outdir"]).arg(&out);
        if case.check_updated {
            cmd.arg("--updated_sam");
        }
        cmd.args(&case.args).arg(dir.join(&case.alignment)).arg(dir.join(&case.gtf));
        let run = cmd.output().unwrap();
        if !run.status.success() {
            failures.push(format!("{}: exited with {}: {}", case.name, run.status, String::from_utf8_lossy(&run.stderr)));
            continue;
        }

        let got = fs::read_to_string(out.join("t-telescope_report.tsv")).unwrap();
        let want = fs::read_to_string(dir.join(format!("expected/{}.report.tsv", case.name))).unwrap();
        if got != want {
            failures.push(format!("{}: report differs at {}", case.name, first_difference(&got, &want)));
        }
        if case.check_updated {
            let got = updated_table(&out.join("t-updated.bam"));
            let want = fs::read_to_string(dir.join(format!("expected/{}.updated.tsv", case.name))).unwrap();
            if got != want {
                failures.push(format!("{}: updated BAM differs at {}", case.name, first_difference(&got, &want)));
            }
        }
    }
    fs::remove_dir_all(&scratch).ok();
    assert!(failures.is_empty(), "{} of {} cases failed:\n{}", failures.len(), all.len(), failures.join("\n"));
}
