use rust_lapper::{Interval, Lapper};
use std::collections::{HashMap, HashSet};
use std::fs;
use std::io;
use std::sync::Mutex;

use crate::cytoband::BandIndex;
use crate::intern::{Interner, Symbol};
use crate::overlap::NO_FEATURE;

/// Metadata from a single GTF feature line.
#[derive(Debug, Clone, Eq, PartialEq)]
pub struct FeatureData {
    pub attributes: HashMap<String, Symbol>,
    /// Strand from GTF column 7. `Some('+')`, `Some('-')`, or `None` for `.`/unknown.
    pub strand: Option<char>,
}

/// A GTF interval suitable for rust-lapper. Coordinates are [start, stop) half-open.
pub type FeatureInterval = Interval<usize, FeatureData>;

/// Default length assigned to the synthetic `__no_feature__` annotation. Set
/// to a megabase so it competes meaningfully against real (much shorter)
/// annotations under length-correction.
pub const NO_FEATURE_LENGTH: u64 = 1_000_000;

/// Index of GTF features organized by chromosome, backed by Lapper for fast overlap queries.
pub struct FeatureIndex {
    lappers: HashMap<String, Lapper<usize, FeatureData>>,
    /// Track chromosomes we've already warned about for chr-prefix mismatch.
    warned: Mutex<HashSet<String>>,
    /// String interner for attribute values (annotation names).
    pub interner: Interner,
    /// Symbol for the synthetic `__no_feature__` annotation. Interned at startup
    /// so it has a stable index in the interner's dense space and can participate
    /// in EM exactly like a real annotation.
    pub no_feature_symbol: Symbol,
    /// Per-Symbol union-of-intervals length in base pairs, indexed by
    /// `Symbol::idx()`. Built from the `field` attribute passed to `from_file`.
    /// `__no_feature__` is set to `NO_FEATURE_LENGTH`. Symbols never seen in a
    /// matching GTF row remain at `0`. Extended with per-band lengths by
    /// `attach_bands`.
    pub lengths: Vec<u64>,
    /// Cytoband index for partitioning no-feature reads into genomic bands.
    /// `None` until `attach_bands` is called (band mode off).
    pub bands: Option<BandIndex>,
}

/// Compact, `Copy` classifier for "is this annotation a flavour of no-feature?"
/// — i.e. the generic `__no_feature__` symbol or any cytoband symbol. Threaded
/// into `compute_summary_counts` and `build_output_directives` so band picks are
/// aggregated into the run-summary's no-feature count and gated out of the BAM
/// output the same way `__no_feature__` is.
#[derive(Debug, Clone, Copy)]
pub struct NoFeatureClass {
    no_feature_symbol: Symbol,
    /// `[start, end)` dense-index range of band symbols, if band mode is on.
    band_range: Option<(usize, usize)>,
}

impl NoFeatureClass {
    /// Whether `sym` is the generic `__no_feature__` symbol or any band symbol.
    pub fn is_no_feature(&self, sym: Symbol) -> bool {
        if sym == self.no_feature_symbol {
            return true;
        }
        matches!(self.band_range, Some((start, end)) if sym.idx() >= start && sym.idx() < end)
    }

    /// The generic `__no_feature__` symbol, used as the ZB/ZF fallback.
    pub fn no_feature_symbol(&self) -> Symbol {
        self.no_feature_symbol
    }
}

/// Accumulator for `FeatureIndex` construction. Lives only during load; its
/// per-Symbol interval lists are merged into union lengths in `finish`.
struct IndexBuilder {
    by_chrom: HashMap<String, Vec<FeatureInterval>>,
    interner: Interner,
    no_feature_symbol: Symbol,
    /// (chrom, lapper_start, lapper_stop) intervals per attribute-value Symbol.
    /// Only populated for the `field` we're tagging on.
    intervals_by_symbol: HashMap<Symbol, Vec<(String, usize, usize)>>,
}

impl IndexBuilder {
    fn new() -> Self {
        let mut interner = Interner::new();
        let no_feature_symbol = interner.intern(NO_FEATURE);
        IndexBuilder {
            by_chrom: HashMap::new(),
            interner,
            no_feature_symbol,
            intervals_by_symbol: HashMap::new(),
        }
    }

    fn ingest_line(&mut self, line: &str, field: &str) {
        let Some((chrom, iv)) = parse_gtf_line(line, &mut self.interner) else {
            return;
        };
        // Track the (chrom, start, stop) interval under the field's symbol so
        // we can compute union-of-intervals length at finish time. Only the
        // single attribute we're tagging on contributes to length.
        if let Some(&sym) = iv.val.attributes.get(field) {
            self.intervals_by_symbol.entry(sym).or_default().push((
                chrom.clone(),
                iv.start,
                iv.stop,
            ));
        }
        self.by_chrom.entry(chrom).or_default().push(iv);
    }

    fn finish(self) -> FeatureIndex {
        let lappers: HashMap<String, Lapper<usize, FeatureData>> = self
            .by_chrom
            .into_iter()
            .map(|(chrom, intervals)| (chrom, Lapper::new(intervals)))
            .collect();

        // Allocate length vector sized to the interner's full width. Symbols
        // not present in `intervals_by_symbol` (e.g. attributes other than the
        // tagging field) remain at 0 — they won't be queried anyway.
        let mut lengths = vec![0u64; self.interner.len()];
        for (sym, intervals) in self.intervals_by_symbol {
            lengths[sym.idx()] = union_length(intervals);
        }
        // Synthetic __no_feature__ never appears in a GTF row; assign it a
        // fixed length so length-correction in EM treats it as a long
        // background.
        lengths[self.no_feature_symbol.idx()] = NO_FEATURE_LENGTH;

        FeatureIndex {
            lappers,
            warned: Mutex::new(HashSet::new()),
            interner: self.interner,
            no_feature_symbol: self.no_feature_symbol,
            lengths,
            bands: None,
        }
    }
}

/// Compute the total number of base pairs covered by the union of a set of
/// intervals, grouped per chromosome so cross-chrom overlaps don't merge.
fn union_length(mut intervals: Vec<(String, usize, usize)>) -> u64 {
    intervals.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
    let mut total: u64 = 0;
    let mut cur_chrom: Option<String> = None;
    let mut cur_start: usize = 0;
    let mut cur_stop: usize = 0;
    for (chrom, start, stop) in intervals {
        if cur_chrom.as_ref() == Some(&chrom) && start <= cur_stop {
            // Overlapping or adjacent on the same chrom — extend.
            if stop > cur_stop {
                cur_stop = stop;
            }
        } else {
            if cur_chrom.is_some() {
                total += (cur_stop - cur_start) as u64;
            }
            cur_chrom = Some(chrom);
            cur_start = start;
            cur_stop = stop;
        }
    }
    if cur_chrom.is_some() {
        total += (cur_stop - cur_start) as u64;
    }
    total
}

impl FeatureIndex {
    /// Load a GTF file and build the feature index.
    /// Reads line-by-line to avoid loading the entire file into memory.
    /// `field` is the GTF attribute used both for tagging and for computing
    /// per-Symbol total length (union of intervals tagged with that attribute).
    pub fn from_file(path: &str, field: &str) -> io::Result<Self> {
        use std::io::BufRead;
        let file = fs::File::open(path)?;
        let reader = io::BufReader::new(file);
        let mut builder = IndexBuilder::new();

        for line in reader.lines() {
            let line = line?;
            if line.starts_with('#') || line.is_empty() {
                continue;
            }
            builder.ingest_line(&line, field);
        }

        Ok(builder.finish())
    }

    /// Build a feature index from GTF content string (used in tests).
    #[cfg(test)]
    pub fn from_str(contents: &str, field: &str) -> Self {
        let mut builder = IndexBuilder::new();
        for line in contents.lines() {
            if line.starts_with('#') || line.is_empty() {
                continue;
            }
            builder.ingest_line(line, field);
        }
        builder.finish()
    }

    /// Find overlapping features on the given chromosome.
    /// Returns an iterator of intervals overlapping [start, stop).
    ///
    /// If `normalize_chr` is true, tries adding/removing "chr" prefix on miss.
    /// If `normalize_chr` is false but a chr-prefix variant would match, warns once to stderr.
    pub fn find(
        &self,
        chrom: &str,
        start: usize,
        stop: usize,
        normalize_chr: bool,
    ) -> Vec<&FeatureInterval> {
        // Try exact match first
        if let Some(lapper) = self.lappers.get(chrom) {
            return lapper.find(start, stop).collect();
        }

        // Try normalized variant
        let alt = chr_alternate(chrom);
        if let Some(lapper) = self.lappers.get(&alt) {
            if normalize_chr {
                return lapper.find(start, stop).collect();
            } else {
                // Warn once per chromosome
                let mut warned = self.warned.lock().unwrap();
                if warned.insert(chrom.to_string()) {
                    eprintln!(
                        "Warning: chromosome '{}' not found in GTF, but '{}' exists. \
                         Use --normalize-chr to enable automatic matching.",
                        chrom, alt
                    );
                }
            }
        }

        Vec::new()
    }

    /// Load a cytoBand file, interning each band into this index's interner and
    /// extending `lengths` with the per-band lengths. Enables band mode: future
    /// `assign_band` calls return a band symbol instead of `None`.
    ///
    /// Must be called after construction but before the parallel phases that
    /// read `interner.len()` / `lengths` (e.g. EM θ sizing).
    pub fn attach_bands(&mut self, path: &str) -> io::Result<()> {
        let bands = BandIndex::from_file(path, &mut self.interner)?;
        // The interner grew while interning bands; widen `lengths` to match and
        // fill in each band's real length (chromEnd - chromStart).
        self.lengths.resize(self.interner.len(), 0);
        for &(sym, len) in bands.lengths() {
            self.lengths[sym.idx()] = len;
        }
        self.bands = Some(bands);
        Ok(())
    }

    /// Assign a no-feature read's spans to its best-overlapping cytoband.
    /// Returns `None` when band mode is off or no span hits a band.
    pub fn assign_band(&self, spans: &[(&str, usize, usize)]) -> Option<(Symbol, usize)> {
        self.bands.as_ref()?.assign_band(spans)
    }

    /// Build the `Copy` classifier used to recognise no-feature/band symbols
    /// during summary counting and output filtering.
    pub fn no_feature_class(&self) -> NoFeatureClass {
        NoFeatureClass {
            no_feature_symbol: self.no_feature_symbol,
            band_range: self.bands.as_ref().map(|b| b.band_range()),
        }
    }

    /// Number of cytobands loaded (0 when band mode is off).
    pub fn num_bands(&self) -> usize {
        self.bands.as_ref().map_or(0, |b| b.num_bands())
    }

    /// Number of chromosomes in the index.
    pub fn num_chromosomes(&self) -> usize {
        self.lappers.len()
    }

    /// Total number of feature intervals across all chromosomes.
    pub fn num_features(&self) -> usize {
        self.lappers.values().map(|l| l.len()).sum()
    }

    /// Returns the set of chromosome names present in the index.
    #[cfg(test)]
    #[allow(dead_code)]
    pub fn chromosomes(&self) -> Vec<&String> {
        self.lappers.keys().collect()
    }
}

/// Given "chrX" return "X", given "X" return "chrX".
fn chr_alternate(name: &str) -> String {
    if let Some(stripped) = name.strip_prefix("chr") {
        stripped.to_string()
    } else {
        format!("chr{}", name)
    }
}

/// Parse a single GTF line into (chromosome, Interval).
/// Returns None if the line cannot be parsed.
fn parse_gtf_line(line: &str, interner: &mut Interner) -> Option<(String, FeatureInterval)> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 9 {
        return None;
    }

    let chrom = fields[0].to_string();
    let start: usize = fields[3].parse().ok()?;
    let end: usize = fields[4].parse().ok()?;

    // GTF is 1-based inclusive [start, end].
    // rust-lapper uses 0-based half-open [start, stop).
    // Convert: lapper_start = start - 1, lapper_stop = end
    // (since GTF start=1 means position 1, which is index 0 in 0-based)
    let lapper_start = start.checked_sub(1)?;
    let lapper_stop = end;

    let strand = parse_strand_char(fields[6]);
    let attributes = parse_attributes(fields[8], interner);

    Some((
        chrom,
        FeatureInterval {
            start: lapper_start,
            stop: lapper_stop,
            val: FeatureData { attributes, strand },
        },
    ))
}

/// Parse the strand column of a GTF row.
/// Returns `Some('+')` or `Some('-')`; `None` for `.`, empty, or unrecognized.
fn parse_strand_char(s: &str) -> Option<char> {
    match s {
        "+" => Some('+'),
        "-" => Some('-'),
        _ => None,
    }
}

/// Parse the GTF attributes column.
/// Format: `key "value"; key "value"; ...`
/// This file has doubled quotes: `key ""value""`
fn parse_attributes(attr_str: &str, interner: &mut Interner) -> HashMap<String, Symbol> {
    let mut attrs = HashMap::new();

    for pair in attr_str.split(';') {
        let pair = pair.trim();
        if pair.is_empty() {
            continue;
        }

        // Split on first whitespace to separate key from value
        if let Some(space_idx) = pair.find(|c: char| c.is_whitespace()) {
            let key = &pair[..space_idx];
            let value = pair[space_idx..].trim();

            // Strip outer quotes (handles both "value" and ""value"")
            let value = value.trim_matches('"').trim();

            attrs.insert(key.to_string(), interner.intern(value));
        }
    }

    attrs
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_attributes_doubled_quotes() {
        let mut interner = Interner::new();
        let attr = r#"gene_id ""HML2.LTR1""; transcript_id ""HML2.LTR1""; repLeft ""0"""#;
        let parsed = parse_attributes(attr, &mut interner);
        assert_eq!(
            interner.resolve(*parsed.get("gene_id").unwrap()),
            "HML2.LTR1"
        );
        assert_eq!(
            interner.resolve(*parsed.get("transcript_id").unwrap()),
            "HML2.LTR1"
        );
        assert_eq!(interner.resolve(*parsed.get("repLeft").unwrap()), "0");
    }

    #[test]
    fn test_parse_attributes_normal_quotes() {
        let mut interner = Interner::new();
        let attr = r#"gene_id "GENE1"; transcript_id "TX1""#;
        let parsed = parse_attributes(attr, &mut interner);
        assert_eq!(interner.resolve(*parsed.get("gene_id").unwrap()), "GENE1");
        assert_eq!(
            interner.resolve(*parsed.get("transcript_id").unwrap()),
            "TX1"
        );
    }

    #[test]
    fn test_parse_gtf_line() {
        let mut interner = Interner::new();
        let line =
            "chr1\tsource\texon\t100\t200\t.\t+\t.\tgene_id \"GENE1\"; transcript_id \"TX1\"";
        let (chrom, iv) = parse_gtf_line(line, &mut interner).unwrap();
        assert_eq!(chrom, "chr1");
        // GTF 1-based [100, 200] -> 0-based [99, 200)
        assert_eq!(iv.start, 99);
        assert_eq!(iv.stop, 200);
        assert_eq!(iv.val.strand, Some('+'));
        assert_eq!(
            interner.resolve(*iv.val.attributes.get("gene_id").unwrap()),
            "GENE1"
        );
    }

    #[test]
    fn test_parse_gtf_line_strand_variants() {
        let mut interner = Interner::new();
        let cases = [
            ("+", Some('+')),
            ("-", Some('-')),
            (".", None),
            ("?", None),
            ("", None),
        ];
        for (col, expected) in cases {
            let line = format!("chr1\tsrc\texon\t100\t200\t.\t{col}\t.\tgene_id \"G\"");
            let (_, iv) = parse_gtf_line(&line, &mut interner).unwrap();
            assert_eq!(iv.val.strand, expected, "strand col {col:?}");
        }
    }

    #[test]
    fn test_feature_index_find() {
        let gtf = "\
chr1\tsrc\texon\t100\t200\t.\t+\t.\tgene_id \"A\"
chr1\tsrc\texon\t300\t400\t.\t+\t.\tgene_id \"B\"
chr2\tsrc\texon\t50\t150\t.\t+\t.\tgene_id \"C\"";

        let index = FeatureIndex::from_str(gtf, "gene_id");

        // Query overlapping feature A (0-based [99, 200))
        let hits: Vec<_> = index.find("chr1", 99, 200, false);
        assert_eq!(hits.len(), 1);
        assert_eq!(
            index
                .interner
                .resolve(*hits[0].val.attributes.get("gene_id").unwrap()),
            "A"
        );

        // Query overlapping nothing on chr1
        let hits: Vec<_> = index.find("chr1", 210, 290, false);
        assert_eq!(hits.len(), 0);

        // Query on missing chromosome
        let hits: Vec<_> = index.find("chr3", 0, 1000, false);
        assert_eq!(hits.len(), 0);
    }

    #[test]
    fn test_chr_normalization() {
        let gtf = "1\tsrc\texon\t100\t200\t.\t+\t.\tgene_id \"A\"";
        let index = FeatureIndex::from_str(gtf, "gene_id");

        // Exact match works
        let hits = index.find("1", 99, 200, false);
        assert_eq!(hits.len(), 1);

        // Without normalize_chr, "chr1" won't match
        let hits = index.find("chr1", 99, 200, false);
        assert_eq!(hits.len(), 0);

        // With normalize_chr, "chr1" finds "1"
        let hits = index.find("chr1", 99, 200, true);
        assert_eq!(hits.len(), 1);
    }

    /// Test helper: find the Symbol for an annotation name in a built index.
    fn lookup_symbol(index: &FeatureIndex, name: &str) -> Symbol {
        // Walk the FeatureData attributes to find the Symbol — simpler and
        // safer than re-interning into a now-shared interner.
        for lapper in index.lappers.values() {
            for iv in lapper.iter() {
                if let Some(&sym) = iv.val.attributes.get("gene_id")
                    && index.interner.resolve(sym) == name
                {
                    return sym;
                }
            }
        }
        panic!("symbol {name:?} not found in index");
    }

    #[test]
    fn test_lengths_union_of_intervals() {
        // Two overlapping exons of GENE_A on chr1:
        //   [100, 200] (1-based inclusive) → 0-based [99, 200), length 101
        //   [150, 250] → 0-based [149, 250), length 101
        // Union: [99, 250), length 151.
        // Plus a separate non-overlapping exon for GENE_B at [300, 400] = 101 bp.
        let gtf = "\
chr1\tsrc\texon\t100\t200\t.\t+\t.\tgene_id \"GENE_A\"
chr1\tsrc\texon\t150\t250\t.\t+\t.\tgene_id \"GENE_A\"
chr1\tsrc\texon\t300\t400\t.\t+\t.\tgene_id \"GENE_B\"";
        let index = FeatureIndex::from_str(gtf, "gene_id");
        let a = lookup_symbol(&index, "GENE_A");
        let b = lookup_symbol(&index, "GENE_B");
        assert_eq!(index.lengths[a.idx()], 151);
        assert_eq!(index.lengths[b.idx()], 101);
        assert_eq!(
            index.lengths[index.no_feature_symbol.idx()],
            NO_FEATURE_LENGTH
        );
    }

    #[test]
    fn test_lengths_disjoint_chroms_summed() {
        // Same gene on two chromosomes — union per chrom, summed across.
        let gtf = "\
chr1\tsrc\texon\t100\t199\t.\t+\t.\tgene_id \"GENE_A\"
chr2\tsrc\texon\t1\t100\t.\t+\t.\tgene_id \"GENE_A\"";
        let index = FeatureIndex::from_str(gtf, "gene_id");
        let a = lookup_symbol(&index, "GENE_A");
        // chr1: [99, 199) = 100 bp; chr2: [0, 100) = 100 bp; total = 200.
        assert_eq!(index.lengths[a.idx()], 200);
    }

    #[test]
    fn test_chr_alternate() {
        assert_eq!(chr_alternate("chr1"), "1");
        assert_eq!(chr_alternate("chrX"), "X");
        assert_eq!(chr_alternate("1"), "chr1");
        assert_eq!(chr_alternate("X"), "chrX");
    }
}
