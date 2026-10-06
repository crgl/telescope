//! Reassignment modes and the `telescope_report.tsv` writer
//! (`TelescopeLikelihood.reassign` and `Telescope.output_report`).

use std::fs::File;
use std::io::{self, BufWriter, Write};

use super::loader::RunInfo;
use super::model::{ModelFit, ScoreMatrix};
use super::numpy::{Mt19937, fmt_f2, fmt_g3, recip0, reduceat_sum};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReassignMode {
    Exclude,
    Choose,
    Average,
    Conf,
    Unique,
}

/// Per-feature totals for one reassignment mode. Telescope prints some modes
/// as integers and others with two decimals.
pub enum Column {
    Int(Vec<i64>),
    Float(Vec<f64>),
}

impl Column {
    fn cell(&self, j: usize, as_float: bool) -> String {
        match self {
            Column::Int(v) if as_float => fmt_f2(v[j] as f64),
            Column::Int(v) => v[j].to_string(),
            Column::Float(v) => fmt_f2(v[j]),
        }
    }
    fn key(&self, j: usize) -> f64 {
        match self {
            Column::Int(v) => v[j] as f64,
            Column::Float(v) => v[j],
        }
    }
}

/// Membership weights plus which entries the model dropped (empty = none).
#[derive(Clone, Copy)]
pub struct Weights<'a> {
    pub z: &'a [f64],
    pub absent: &'a [bool],
}

impl Weights<'_> {
    fn present(&self, k: usize) -> bool {
        self.absent.is_empty() || !self.absent[k]
    }
}

/// For each row, the positions of its best hits (`binmax(1)`), handed to `f`.
fn for_best_hits(m: &ScoreMatrix, w: Weights, best: &mut Vec<usize>, mut f: impl FnMut(&[usize])) {
    let z = w.z;
    for i in 0..m.n_rows {
        let r = m.row(i);
        let mut mx = f64::NEG_INFINITY;
        let mut stored = 0;
        for k in r.clone().filter(|&k| w.present(k)) {
            mx = mx.max(z[k]);
            stored += 1;
        }
        if stored < m.n_cols {
            mx = mx.max(0.0); // implicit zeros take part in a sparse max
        }
        best.clear();
        best.extend(r.filter(|&k| w.present(k) && z[k] == mx));
        f(best);
    }
}

pub fn reassign(
    m: &ScoreMatrix,
    w: Weights,
    mode: ReassignMode,
    thresh: f64,
    rng: &mut Mt19937,
) -> Column {
    let col = |k: usize| m.indices[k] as usize;
    let z = w.z;
    let mut best = Vec::new();
    match mode {
        ReassignMode::Exclude => {
            let mut out = vec![0i64; m.n_cols];
            for_best_hits(m, w, &mut best, |b| {
                if b.len() == 1 {
                    out[col(b[0])] += 1;
                }
            });
            Column::Int(out)
        }
        ReassignMode::Choose => {
            let mut out = vec![0i64; m.n_cols];
            for_best_hits(m, w, &mut best, |b| match b.len() {
                0 => {}
                1 => out[col(b[0])] += 1,
                n => out[col(b[rng.choice(n as u32) as usize])] += 1,
            });
            Column::Int(out)
        }
        ReassignMode::Average => {
            let mut out = vec![0.0; m.n_cols];
            for_best_hits(m, w, &mut best, |b| {
                let share = recip0(b.len() as f64);
                for &k in b {
                    out[col(k)] += share;
                }
            });
            Column::Float(out)
        }
        ReassignMode::Conf => {
            let mut out = vec![0.0; m.n_cols];
            let mut kept = Vec::new();
            for i in 0..m.n_rows {
                let r = m.row(i);
                kept.clear();
                kept.extend(
                    r.clone().filter(|&k| w.present(k)).map(|k| if z[k] >= thresh { z[k] } else { 0.0 }),
                );
                let recip = recip0(reduceat_sum(&kept));
                for (k, &v) in r.filter(|&k| w.present(k)).zip(&kept) {
                    out[col(k)] += v * recip;
                }
            }
            Column::Float(out)
        }
        ReassignMode::Unique => {
            let mut out = vec![0i64; m.n_cols];
            for i in 0..m.n_rows {
                if !m.is_ambiguous(i) {
                    for k in m.row(i) {
                        out[col(k)] += z[k].ceil() as u8 as i64;
                    }
                }
            }
            Column::Int(out)
        }
    }
}

/// `reassign('all', initial=True)`: every candidate with nonzero weight.
fn count_all(m: &ScoreMatrix, z: &[f64]) -> Column {
    let mut out = vec![0i64; m.n_cols];
    for (k, &x) in z.iter().enumerate() {
        if x > 0.0 {
            out[m.indices[k] as usize] += 1;
        }
    }
    Column::Int(out)
}

pub struct ReportInputs<'a> {
    pub matrix: &'a ScoreMatrix,
    pub fit: &'a ModelFit,
    pub feat_names: &'a [String],
    pub feat_lengths: &'a [u64],
    pub info: &'a RunInfo,
    pub version: &'a str,
    pub mode: ReassignMode,
    pub conf_prob: f64,
    /// `Telescope.get_random_seed()`
    pub seed: u32,
}

pub fn write_report(path: &str, r: &ReportInputs) -> io::Result<()> {
    let (m, fit) = (r.matrix, r.fit);
    // Telescope evaluates these in this order; only `choose` draws random
    // numbers, so the order fixes which draws each column sees.
    let mut rng = Mt19937::new(r.seed);
    let fin = Weights { z: &fit.z, absent: &fit.absent };
    let init = Weights { z: &fit.z_init, absent: &[] };
    let final_count = reassign(m, fin, r.mode, r.conf_prob, &mut rng);
    let final_conf = reassign(m, fin, ReassignMode::Conf, r.conf_prob, &mut rng);
    let init_aligned = count_all(m, &fit.z_init);
    let unique_count = reassign(m, fin, ReassignMode::Unique, r.conf_prob, &mut rng);
    let init_best = reassign(m, init, ReassignMode::Exclude, r.conf_prob, &mut rng);
    let init_random = reassign(m, init, ReassignMode::Choose, r.conf_prob, &mut rng);
    let init_avg = reassign(m, init, ReassignMode::Average, r.conf_prob, &mut rng);

    // Two stable descending sorts: by final_prop, then by final_count.
    let mut order: Vec<usize> = (0..m.n_cols).collect();
    order.sort_by(|&a, &b| fit.pi[b].partial_cmp(&fit.pi[a]).unwrap_or(std::cmp::Ordering::Equal));
    order.sort_by(|&a, &b| {
        final_count.key(b).partial_cmp(&final_count.key(a)).unwrap_or(std::cmp::Ordering::Equal)
    });

    let final_as_float = matches!(r.mode, ReassignMode::Average | ReassignMode::Conf);
    let i = r.info;
    let mut out = BufWriter::new(File::create(path)?);
    writeln!(
        out,
        "## RunInfo\tversion:{}\tannotated_features:{}\ttotal_fragments:{}\tpair_mapped:{}\t\
         pair_mixed:{}\tsingle_mapped:{}\tunmapped:{}\tunique:{}\tambig:{}\toverlap_unique:{}\t\
         overlap_ambig:{}",
        r.version,
        i.annotated_features,
        i.total_fragments,
        i.pair_mapped,
        i.pair_mixed,
        i.single_mapped,
        i.unmapped,
        i.unique,
        i.ambig,
        i.overlap_unique,
        i.overlap_ambig
    )?;
    writeln!(
        out,
        "transcript\ttranscript_length\tfinal_count\tfinal_conf\tfinal_prop\tinit_aligned\t\
         unique_count\tinit_best\tinit_best_random\tinit_best_avg\tinit_prop"
    )?;
    for j in order {
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            r.feat_names[j],
            r.feat_lengths[j],
            final_count.cell(j, final_as_float),
            final_conf.cell(j, true),
            fmt_g3(fit.pi[j]),
            init_aligned.cell(j, false),
            unique_count.cell(j, false),
            init_best.cell(j, false),
            init_random.cell(j, false),
            init_avg.cell(j, true),
            fmt_g3(fit.pi_init[j]),
        )?;
    }
    out.flush()
}
