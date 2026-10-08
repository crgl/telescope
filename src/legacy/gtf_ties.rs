//! `gtf-ties`: where an annotation makes `assign` choose arbitrarily.
//!
//! Wherever two or more loci cover the same bases, a read lying wholly inside
//! the shared stretch overlaps them all equally, and Telescope (hence
//! `assign`) gives it to whichever locus its interval tree yields first.
//! That order follows Python's hash-table layout rather than anything
//! biological. This command lists every such stretch, the loci competing for
//! it and the one that wins, so the effect of an annotation can be judged
//! before any reads are run.
//!
//! The winner reported is for a read lying inside the stretch, away from its
//! first base. Telescope's order also depends on the interval boundaries a
//! read spans, so a read that starts exactly on the first base can resolve
//! differently; that case is reported separately (`winner_at_first_base`).
//! A read that extends beyond the stretch overlaps the loci unequally and
//! goes to the larger overlap as usual, unless they share the extra bases.

use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

use clap::Parser;

use super::TieHashArg;
use super::annotation::{Annotation, OverlapCompat, OverlapRules};
use super::pyset;
use crate::logging::{Logger, Verbosity};

#[derive(Parser, Debug)]
pub struct GtfTiesArgs {
    /// Annotation file (GTF)
    pub gtffile: String,

    /// GTF attribute that defines a locus (note: `assign` defaults to `locus`)
    #[arg(long, default_value = "gene_id")]
    pub attribute: String,

    /// Which Python's hashing decides the winner: Telescope under Python 3.8
    /// or newer (`python38`), or under 3.7 or older (`python37`)
    #[arg(long = "tie_hash", value_enum, default_value = "python38")]
    pub tie_hash: TieHashArg,

    /// Only report competition between loci on the same strand, as happens
    /// when `assign` runs with a --stranded_mode
    #[arg(long = "same_strand")]
    pub same_strand: bool,

    /// Output directory
    #[arg(long, default_value = ".")]
    pub outdir: String,

    /// Output file prefix [default: the GTF file name without its extension]
    #[arg(long = "exp_tag")]
    pub exp_tag: Option<String>,

    /// Silence progress output
    #[arg(long)]
    pub quiet: bool,
}

#[derive(Default, Clone)]
struct LocusTally {
    shared: u64,
    won: u64,
    beats: HashMap<u32, u64>,
    loses_to: HashMap<u32, u64>,
}

pub fn run_gtf_ties(args: GtfTiesArgs) -> io::Result<()> {
    let log = Logger::new(if args.quiet { Verbosity::Quiet } else { Verbosity::Normal });
    pyset::set_pre38_tuple_hash(args.tie_hash == TieHashArg::Python37);
    let rules = OverlapRules { coords: OverlapCompat::Telescope, ties: OverlapCompat::Telescope };
    let annot = Annotation::from_gtf(&args.gtffile, &args.attribute, false, rules)?;
    log.stage(&format!("Loaded {} loci.", annot.loci.len()));

    let tag = args.exp_tag.clone().unwrap_or_else(|| {
        Path::new(&args.gtffile).file_stem().map_or_else(|| "annotation".into(), |s| s.to_string_lossy().into_owned())
    });
    let out_path = |suffix: &str| Path::new(&args.outdir).join(format!("{tag}-{suffix}"));
    let regions_path = out_path("tie_regions.tsv");
    let mut regions = BufWriter::new(File::create(&regions_path)?);
    writeln!(regions, "chrom\tstart\tend\tlength\twinner\twinner_strand\tn_competing\tlosers\twinner_at_first_base")?;

    let mut tally: Vec<LocusTally> = vec![LocusTally::default(); annot.loci.len()];
    let (mut n_regions, mut total_bp) = (0u64, 0u64);
    let mut failed: Option<io::Error> = None;
    let mut group: Vec<(u32, u8)> = Vec::new();
    let mut n_position_dependent = 0u64;
    annot.shared_regions(|chrom, begin, end, loci, at_start| {
        // One competition per strand when strands are kept apart, else one for all.
        let strands: Vec<Option<u8>> = if args.same_strand {
            let mut s: Vec<u8> = loci.iter().map(|l| l.1).collect();
            s.sort_unstable();
            s.dedup();
            s.into_iter().map(Some).collect()
        } else {
            vec![None]
        };
        for strand in strands {
            group.clear();
            group.extend(loci.iter().copied().filter(|l| strand.is_none_or(|s| l.1 == s)));
            if group.len() < 2 {
                continue;
            }
            let len = end - begin;
            let (winner, winner_strand) = group[0];
            n_regions += 1;
            total_bp += len;
            for &(l, _) in &group {
                tally[l as usize].shared += len;
            }
            tally[winner as usize].won += len;
            for &(loser, _) in &group[1..] {
                *tally[winner as usize].beats.entry(loser).or_default() += len;
                *tally[loser as usize].loses_to.entry(winner).or_default() += len;
            }
            let edge = at_start.iter().copied().find(|l| strand.is_none_or(|s| l.1 == s)).map_or(winner, |l| l.0);
            n_position_dependent += (edge != winner) as u64;
            let losers: Vec<String> =
                group[1..].iter().map(|&(l, s)| format!("{}({})", annot.loci[l as usize], s as char)).collect();
            // stored [begin, end) is GTF bases begin..=end-1
            let line = writeln!(
                regions,
                "{chrom}\t{begin}\t{}\t{len}\t{}\t{}\t{}\t{}\t{}",
                end - 1,
                annot.loci[winner as usize],
                winner_strand as char,
                group.len(),
                losers.join(","),
                annot.loci[edge as usize]
            );
            if let Err(e) = line {
                failed.get_or_insert(e);
            }
        }
    });
    if let Some(e) = failed {
        return Err(e);
    }
    regions.flush()?;

    let loci_path = out_path("tie_loci.tsv");
    let mut out = BufWriter::new(File::create(&loci_path)?);
    writeln!(out, "locus\tlength\tshared_bp\tshared_fraction\twon_bp\tlost_bp\twins_over\tloses_to")?;
    let list = |m: &HashMap<u32, u64>| -> String {
        let mut v: Vec<(&u32, &u64)> = m.iter().collect();
        v.sort_by(|a, b| b.1.cmp(a.1).then(a.0.cmp(b.0)));
        v.iter().map(|(l, bp)| format!("{}:{}", annot.loci[**l as usize], bp)).collect::<Vec<_>>().join(",")
    };
    let mut affected = 0u64;
    let mut never_win = 0u64;
    for (id, t) in tally.iter().enumerate() {
        if t.shared == 0 {
            continue;
        }
        affected += 1;
        // a locus that is shared along its whole length and never wins gets no
        // read that lies wholly inside it
        if t.won == 0 && t.shared >= annot.lengths[id] {
            never_win += 1;
        }
        writeln!(
            out,
            "{}\t{}\t{}\t{:.4}\t{}\t{}\t{}\t{}",
            annot.loci[id],
            annot.lengths[id],
            t.shared,
            t.shared as f64 / annot.lengths[id].max(1) as f64,
            t.won,
            t.shared - t.won,
            list(&t.beats),
            list(&t.loses_to)
        )?;
    }
    out.flush()?;
    log.stage(&format!(
        "{n_regions} shared regions covering {total_bp} bp; {affected} of {} loci involved; \
         {never_win} are shared along their whole length and never win",
        annot.loci.len()
    ));
    log.stage(&format!(
        "In {n_position_dependent} regions a read starting exactly on the first base goes to a different locus"
    ));
    log.stage(&format!("Wrote {} and {}", regions_path.display(), loci_path.display()));
    Ok(())
}
