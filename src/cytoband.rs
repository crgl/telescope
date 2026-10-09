use std::cmp;
use std::collections::HashMap;
use std::fs;
use std::io::{self, BufRead};

use rust_lapper::{Interval, Lapper};

use crate::intern::{Interner, Symbol};

/// Strip a leading `chr` prefix if present, else return the name unchanged.
///
/// Unlike `gtf::chr_alternate` (which *toggles* the prefix), this always
/// normalizes toward the bare form so `"chr1"` and `"1"` collapse to the same
/// key. Bands are always chr-normalized regardless of `--normalize-chr`, and
/// the bare form is also what appears in the synthetic band name.
pub fn strip_chr(name: &str) -> &str {
    name.strip_prefix("chr").unwrap_or(name)
}

/// Extract a GenBank/RefSeq contig accession from a reference name: the first
/// run of two ASCII-uppercase letters followed by ≥5 ASCII digits (the version
/// and any prefix/suffix are excluded). This is the stable identity shared by a
/// BAM's `GL000220.1` and a cytoBand's `Un_GL000220v1` /
/// `chr14_GL000220v1_random`, so unplaced/alt contigs can be matched across the
/// GenBank-accession vs UCSC-name conventions.
///
/// Returns `None` for names without such a token (e.g. `1`, `chr1`, `MT`), which
/// are matched by `strip_chr` instead and must not be touched by the alias path.
pub fn accession_token(name: &str) -> Option<String> {
    let b = name.as_bytes();
    let n = b.len();
    let mut i = 0;
    while i + 2 < n {
        if b[i].is_ascii_uppercase() && b[i + 1].is_ascii_uppercase() {
            let mut j = i + 2;
            while j < n && b[j].is_ascii_digit() {
                j += 1;
            }
            if j - (i + 2) >= 5 {
                return Some(name[i..j].to_string());
            }
        }
        i += 1;
    }
    None
}

/// A cytoBand interval suitable for rust-lapper. Coordinates are [start, stop)
/// half-open, carrying the interned Symbol of the band's synthetic name.
type BandInterval = Interval<usize, Symbol>;

/// Per-chromosome index of cytogenetic bands parsed from a UCSC `cytoBand.txt`
/// file. Chromosomes are keyed by their chr-stripped name so lookups are always
/// chr-normalized. Each band is interned (into the shared GTF interner) as
/// `__no_feature_<chr><band>__` so it participates in the matrix, priors, EM,
/// Jaccard, and output exactly like a real annotation.
pub struct BandIndex {
    lappers: HashMap<String, Lapper<usize, Symbol>>,
    /// Maps a contig accession (see `accession_token`) to its primary
    /// `strip_chr` chromosome key, so a BAM ref named with a GenBank accession
    /// (`GL000220.1`) resolves to the cytoBand entry that uses a UCSC-style name
    /// (`Un_GL000220v1`). Only populated for chroms that carry an accession.
    accession_alias: HashMap<String, String>,
    /// `[start, end)` half-open range of the interned band symbols' dense
    /// indices. Used for fast "is this symbol a band?" membership tests.
    band_range: (usize, usize),
    /// `(band symbol, length)` pairs, where length is `chromEnd - chromStart`.
    /// Consumed by `FeatureIndex::attach_bands` to populate `lengths`.
    lengths: Vec<(Symbol, u64)>,
    num_bands: usize,
}

impl BandIndex {
    /// Load a cytoBand file and intern every band into `interner`.
    pub fn from_file(path: &str, interner: &mut Interner) -> io::Result<Self> {
        let file = fs::File::open(path)?;
        let reader = io::BufReader::new(file);
        Self::ingest(reader.lines(), interner)
    }

    /// Build from cytoBand content (used in tests).
    #[cfg(test)]
    pub fn from_str(contents: &str, interner: &mut Interner) -> Self {
        let lines = contents.lines().map(|l| Ok(l.to_string()));
        Self::ingest(lines, interner).expect("from_str should not perform I/O")
    }

    fn ingest<I>(lines: I, interner: &mut Interner) -> io::Result<Self>
    where
        I: Iterator<Item = io::Result<String>>,
    {
        let range_start = interner.len();
        let mut by_chrom: HashMap<String, Vec<BandInterval>> = HashMap::new();
        let mut accession_alias: HashMap<String, String> = HashMap::new();
        let mut lengths: Vec<(Symbol, u64)> = Vec::new();
        let mut num_bands = 0usize;

        for line in lines {
            let line = line?;
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let fields: Vec<&str> = line.split('\t').collect();
            // chrom, chromStart, chromEnd, name (gieStain in col 5 is unused).
            if fields.len() < 4 {
                continue;
            }
            // cytoBand is already 0-based half-open (BED-style) — no -1 shift.
            let (Ok(start), Ok(end)) = (fields[1].parse::<usize>(), fields[2].parse::<usize>())
            else {
                continue;
            };
            if end <= start {
                continue;
            }

            let key = strip_chr(fields[0]).to_string();
            // Register the contig accession (if any) as an alias to this chrom
            // key so GenBank-accession-named reads can match. Same accession
            // never appears on two chroms, so later inserts are idempotent.
            if let Some(acc) = accession_token(fields[0]) {
                accession_alias.entry(acc).or_insert_with(|| key.clone());
            }
            // col 4 (band name) may be empty for alt/random/unplaced contigs.
            let name = format!("__no_feature_{}{}__", key, fields[3]);
            let sym = interner.intern(&name);

            by_chrom.entry(key).or_default().push(BandInterval {
                start,
                stop: end,
                val: sym,
            });
            lengths.push((sym, (end - start) as u64));
            num_bands += 1;
        }

        let lappers = by_chrom
            .into_iter()
            .map(|(chrom, intervals)| (chrom, Lapper::new(intervals)))
            .collect();
        let range_end = interner.len();

        Ok(BandIndex {
            lappers,
            accession_alias,
            band_range: (range_start, range_end),
            lengths,
            num_bands,
        })
    }

    /// Resolve a reference name to its band Lapper: try the chr-normalized key
    /// first, then fall back to the contig-accession alias.
    fn lapper_for(&self, chrom: &str) -> Option<&Lapper<usize, Symbol>> {
        if let Some(lapper) = self.lappers.get(strip_chr(chrom)) {
            return Some(lapper);
        }
        let acc = accession_token(chrom)?;
        let primary = self.accession_alias.get(&acc)?;
        self.lappers.get(primary)
    }

    /// Assign one or more genomic spans (a single read or a mate pair) to the
    /// single band with the greatest total overlap. Overlap is summed per band
    /// across the spans, then the max is taken — one band per read ("partition").
    ///
    /// Equal-overlap ties are broken deterministically by the lower band
    /// `Symbol` index (cytoBand file order). Without this, the winner would
    /// depend on `HashMap` iteration order and vary run-to-run, making
    /// band-band Jaccard output non-reproducible.
    ///
    /// Returns `None` when no span hits any band (unmapped, empty spans, or a
    /// contig absent from the cytoBand file), so the caller falls back to the
    /// generic `__no_feature__` path.
    pub fn assign_band(&self, spans: &[(&str, usize, usize)]) -> Option<(Symbol, usize)> {
        let mut overlap_by_band: HashMap<Symbol, usize> = HashMap::new();
        for &(chrom, start, stop) in spans {
            let Some(lapper) = self.lapper_for(chrom) else {
                continue;
            };
            for iv in lapper.find(start, stop) {
                let overlap = cmp::min(stop, iv.stop).saturating_sub(cmp::max(start, iv.start));
                if overlap == 0 {
                    continue;
                }
                *overlap_by_band.entry(iv.val).or_insert(0) += overlap;
            }
        }
        // Highest overlap wins; ties broken toward the lower Symbol index so the
        // result is independent of HashMap iteration order.
        overlap_by_band
            .into_iter()
            .max_by(|a, b| a.1.cmp(&b.1).then_with(|| b.0.idx().cmp(&a.0.idx())))
    }

    /// `[start, end)` dense-index range of the interned band symbols.
    pub fn band_range(&self) -> (usize, usize) {
        self.band_range
    }

    /// `(band symbol, length)` pairs for populating the shared `lengths` vector.
    pub fn lengths(&self) -> &[(Symbol, u64)] {
        &self.lengths
    }

    pub fn num_bands(&self) -> usize {
        self.num_bands
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_strip_chr() {
        assert_eq!(strip_chr("chr1"), "1");
        assert_eq!(strip_chr("chrX"), "X");
        assert_eq!(strip_chr("1"), "1");
        assert_eq!(strip_chr("X"), "X");
        assert_eq!(strip_chr("chr1_KI270706v1_random"), "1_KI270706v1_random");
    }

    #[test]
    fn test_parse_and_naming() {
        let mut interner = Interner::new();
        let bands = BandIndex::from_str("chr1\t0\t2300000\tp36.33\tgneg", &mut interner);
        assert_eq!(bands.num_bands(), 1);
        let (sym, len) = bands.lengths()[0];
        assert_eq!(interner.resolve(sym), "__no_feature_1p36.33__");
        assert_eq!(len, 2_300_000);
    }

    #[test]
    fn test_empty_band_name() {
        let mut interner = Interner::new();
        let bands = BandIndex::from_str("chr1_KI270706v1_random\t0\t175055\t\tgneg", &mut interner);
        let (sym, _) = bands.lengths()[0];
        assert_eq!(interner.resolve(sym), "__no_feature_1_KI270706v1_random__");
    }

    #[test]
    fn test_assign_band_half_open_and_normalization() {
        let mut interner = Interner::new();
        // Built from "chr1"; queries with both "chr1" and "1" must hit.
        let bands = BandIndex::from_str("chr1\t0\t100\tp1\tgneg", &mut interner);
        let p1 = interner.intern("__no_feature_1p1__");

        assert_eq!(bands.assign_band(&[("chr1", 10, 50)]), Some((p1, 40)));
        assert_eq!(bands.assign_band(&[("1", 10, 50)]), Some((p1, 40)));
        // Half-open [0,100): position 100 is outside.
        assert_eq!(bands.assign_band(&[("chr1", 100, 150)]), None);
    }

    #[test]
    fn test_assign_band_max_overlap_across_boundary() {
        let mut interner = Interner::new();
        let bands = BandIndex::from_str(
            "chr1\t0\t100\tp1\tgneg\nchr1\t100\t200\tq1\tgneg",
            &mut interner,
        );
        let q1 = interner.intern("__no_feature_1q1__");
        // Span [80,150): 20 bp in p1, 50 bp in q1 → q1 wins.
        assert_eq!(bands.assign_band(&[("chr1", 80, 150)]), Some((q1, 50)));
    }

    #[test]
    fn test_assign_band_sums_across_spans() {
        let mut interner = Interner::new();
        let bands = BandIndex::from_str("chr1\t0\t1000\tp1\tgneg", &mut interner);
        let p1 = interner.intern("__no_feature_1p1__");
        // Two mate spans both land in p1: 40 + 30 = 70 bp.
        assert_eq!(
            bands.assign_band(&[("chr1", 10, 50), ("chr1", 100, 130)]),
            Some((p1, 70))
        );
    }

    #[test]
    fn test_assign_band_tie_is_deterministic() {
        let mut interner = Interner::new();
        let bands = BandIndex::from_str(
            "chr1\t0\t1000\tp1\tgneg\nchr1\t1000\t2000\tp2\tgneg",
            &mut interner,
        );
        let p1 = interner.intern("__no_feature_1p1__");
        // [950,1050) overlaps p1 by 50 (950..1000) and p2 by 50 (1000..1050).
        // Equal overlap → deterministic tie-break picks the lower-index band (p1).
        assert_eq!(bands.assign_band(&[("chr1", 950, 1050)]), Some((p1, 50)));
    }

    #[test]
    fn test_assign_band_unknown_chrom() {
        let mut interner = Interner::new();
        let bands = BandIndex::from_str("chr1\t0\t100\tp1\tgneg", &mut interner);
        assert_eq!(bands.assign_band(&[("chr2", 10, 50)]), None);
        assert_eq!(bands.assign_band(&[]), None);
    }

    #[test]
    fn test_band_range_contiguous() {
        let mut interner = Interner::new();
        interner.intern("some_gene"); // idx 0
        let bands = BandIndex::from_str(
            "chr1\t0\t100\tp1\tgneg\nchr1\t100\t200\tq1\tgneg",
            &mut interner,
        );
        let (start, end) = bands.band_range();
        assert_eq!(start, 1);
        assert_eq!(end, 3);
    }

    #[test]
    fn test_accession_token() {
        // BAM (GenBank) and cytoBand (UCSC) forms reduce to the same accession.
        assert_eq!(accession_token("GL000220.1").as_deref(), Some("GL000220"));
        assert_eq!(
            accession_token("Un_GL000220v1").as_deref(),
            Some("GL000220")
        );
        assert_eq!(
            accession_token("chr14_GL000220v1_random").as_deref(),
            Some("GL000220")
        );
        assert_eq!(
            accession_token("22_KI270733v1_random").as_deref(),
            Some("KI270733")
        );
        assert_eq!(accession_token("KI270733.1").as_deref(), Some("KI270733"));
        // Main chroms and MT have no accession token (matched via strip_chr).
        assert_eq!(accession_token("1"), None);
        assert_eq!(accession_token("chr1"), None);
        assert_eq!(accession_token("chr14"), None);
        assert_eq!(accession_token("X"), None);
        assert_eq!(accession_token("MT"), None);
    }

    #[test]
    fn test_assign_band_matches_by_accession() {
        // cytoBand uses the UCSC-style name; BAM read uses the GenBank accession.
        let mut interner = Interner::new();
        let bands = BandIndex::from_str("Un_GL000220v1\t0\t161802\t\tgneg", &mut interner);
        let sym = interner.intern("__no_feature_Un_GL000220v1__");
        // Direct string key ("GL000220.1") misses; accession alias matches.
        assert_eq!(
            bands.assign_band(&[("GL000220.1", 10, 100)]),
            Some((sym, 90))
        );
        // The band keeps its cytoBand-derived name.
        assert_eq!(interner.resolve(sym), "__no_feature_Un_GL000220v1__");
    }

    #[test]
    fn test_accession_alias_does_not_shadow_numeric_chrom() {
        // A bare numeric chrom must still match its own band, not leak into the
        // accession path (which only fires on a strip_chr miss).
        let mut interner = Interner::new();
        let bands = BandIndex::from_str(
            "1\t0\t100\tp1\tgneg\n22_KI270733v1_random\t0\t179772\t\tgneg",
            &mut interner,
        );
        let p1 = interner.intern("__no_feature_1p1__");
        let kband = interner.intern("__no_feature_22_KI270733v1_random__");
        assert_eq!(bands.assign_band(&[("1", 10, 50)]), Some((p1, 40)));
        assert_eq!(bands.assign_band(&[("chr1", 10, 50)]), Some((p1, 40)));
        assert_eq!(
            bands.assign_band(&[("KI270733.1", 0, 100)]),
            Some((kband, 100))
        );
    }
}
