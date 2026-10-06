use std::collections::{HashMap, HashSet};
use std::io::{self, Write};

use crate::intern::{Interner, Symbol};

/// Per-cell data in the read-annotation matrix: the best alignment
/// for a given (read group, annotation) pair.
#[derive(Debug, Clone, Copy)]
pub struct MatrixCell {
    /// Index into the ReadGroup's alignments vec.
    pub alignment_idx: usize,
    /// AS tag value of the selected alignment (0 if absent).
    pub alignment_score: i64,
}

/// One row of the read-annotation matrix. Stored as a `Vec<(Symbol, MatrixCell)>`
/// sorted ascending by `Symbol` rather than a `HashMap` because rows are tiny
/// (~6 entries on real datasets) — a packed sorted Vec is both smaller (no
/// HashMap header / bucket overhead, ~3× across millions of rows) and faster
/// to iterate. Lookup is via `row_get` (binary search).
pub type MatrixRow = [(Symbol, MatrixCell)];

/// Lookup a cell by annotation symbol. `O(log n)`.
pub fn row_get(row: &MatrixRow, sym: Symbol) -> Option<&MatrixCell> {
    row.binary_search_by_key(&sym.idx(), |(s, _)| s.idx())
        .ok()
        .map(|i| &row[i].1)
}

/// Sparse matrix: rows = read groups, columns = annotations.
/// Each cell stores the best-scoring alignment for that pair.
/// Built during the selection phase, consumed by resolution, Jaccard, and output.
pub struct ReadAnnotationMatrix {
    rows: Vec<Vec<(Symbol, MatrixCell)>>,
}

impl ReadAnnotationMatrix {
    /// Caller is responsible for sorting each row by `Symbol` ascending.
    /// `select_representatives` produces sorted rows already.
    pub fn new(rows: Vec<Vec<(Symbol, MatrixCell)>>) -> Self {
        ReadAnnotationMatrix { rows }
    }

    pub fn row(&self, group_id: usize) -> &MatrixRow {
        &self.rows[group_id]
    }

    pub fn num_groups(&self) -> usize {
        self.rows.len()
    }
}

/// Sparse binary matrix: rows = read groups, columns = annotations.
/// Stores which read groups have a representative for which annotations.
/// Tracks co-occurring annotation pairs for efficient Jaccard computation.
pub struct SparseMatrix {
    annotation_to_col: HashMap<Symbol, usize>,
    col_to_annotation: Vec<Symbol>,
    columns: Vec<HashSet<usize>>,
    /// Pairs (i, j) where i < j that share at least one read group.
    co_occurring: HashSet<(usize, usize)>,
}

impl SparseMatrix {
    pub fn new() -> Self {
        SparseMatrix {
            annotation_to_col: HashMap::new(),
            col_to_annotation: Vec::new(),
            columns: Vec::new(),
            co_occurring: HashSet::new(),
        }
    }

    fn get_or_create_col(&mut self, annotation: Symbol) -> usize {
        if let Some(&col) = self.annotation_to_col.get(&annotation) {
            col
        } else {
            let col = self.col_to_annotation.len();
            self.annotation_to_col.insert(annotation, col);
            self.col_to_annotation.push(annotation);
            self.columns.push(HashSet::new());
            col
        }
    }

    /// Register all annotations for a read group at once.
    /// Tracks which annotation pairs co-occur (share a read group).
    pub fn insert_group(
        &mut self,
        read_group_id: usize,
        annotations: impl Iterator<Item = Symbol>,
    ) {
        let cols: Vec<usize> = annotations.map(|ann| self.get_or_create_col(ann)).collect();

        for &col in &cols {
            self.columns[col].insert(read_group_id);
        }

        // Record all pairwise co-occurrences
        for i in 0..cols.len() {
            for j in (i + 1)..cols.len() {
                let pair = if cols[i] < cols[j] {
                    (cols[i], cols[j])
                } else {
                    (cols[j], cols[i])
                };
                self.co_occurring.insert(pair);
            }
        }
    }

    /// Compute the Jaccard / overlap-coefficient pair for columns i, j.
    /// Returns `None` if either column has fewer than `min_reads` supporting
    /// read groups, or if the overlap coefficient is below `similarity_threshold`.
    fn pair_metrics(
        &self,
        i: usize,
        j: usize,
        similarity_threshold: f64,
        min_reads: usize,
    ) -> Option<(f64, f64)> {
        let a_size = self.columns[i].len();
        let b_size = self.columns[j].len();
        if a_size < min_reads || b_size < min_reads {
            return None;
        }
        let intersection = self.columns[i].intersection(&self.columns[j]).count();
        let union = a_size + b_size - intersection;
        let min_size = a_size.min(b_size);
        let jaccard = intersection as f64 / union as f64;
        let overlap = intersection as f64 / min_size as f64;
        if overlap < similarity_threshold {
            return None;
        }
        Some((jaccard, overlap))
    }

    /// Resolve column index → (a, b) annotation name pair, lexicographically
    /// ordered so the TSV is symmetric regardless of input column order.
    fn resolved_pair<'a>(&self, interner: &'a Interner, i: usize, j: usize) -> (&'a str, &'a str) {
        let a = interner.resolve(self.col_to_annotation[i]);
        let b = interner.resolve(self.col_to_annotation[j]);
        if a < b { (a, b) } else { (b, a) }
    }

    /// Compute pairwise Jaccard similarity between co-occurring annotation columns.
    /// Returns sparse tuples `(annotation_A, annotation_B, jaccard, overlap)`.
    ///
    /// Pairs are dropped when either annotation has fewer than `min_reads`
    /// supporting read groups, or when the overlap coefficient is below
    /// `similarity_threshold`. The remaining tuples are sorted alphabetically
    /// by resolved annotation names.
    pub fn jaccard_tuples(
        &self,
        interner: &Interner,
        similarity_threshold: f64,
        min_reads: usize,
    ) -> Vec<(String, String, f64, f64)> {
        use rayon::prelude::*;

        let pairs: Vec<(usize, usize)> = self.co_occurring.iter().copied().collect();

        let mut tuples: Vec<(String, String, f64, f64)> = pairs
            .par_iter()
            .filter_map(|&(i, j)| {
                let (jaccard, overlap) =
                    self.pair_metrics(i, j, similarity_threshold, min_reads)?;
                let (a, b) = self.resolved_pair(interner, i, j);
                Some((a.to_string(), b.to_string(), jaccard, overlap))
            })
            .collect();

        tuples.sort_by(|a, b| a.0.cmp(&b.0).then_with(|| a.1.cmp(&b.1)));
        tuples
    }

    /// Stream Jaccard pairs to a writer line-by-line — no in-memory tuples
    /// buffer, no sort. Used by `--low-memory` to avoid the O(pairs) Vec
    /// allocation. Iteration order matches `co_occurring`'s HashSet order
    /// (deterministic per process but otherwise arbitrary). The user can
    /// sort the output file post-hoc if needed.
    ///
    /// Writes the standard header line followed by zero or more rows.
    pub fn write_jaccard_streaming<W: Write>(
        &self,
        interner: &Interner,
        similarity_threshold: f64,
        min_reads: usize,
        writer: &mut W,
    ) -> io::Result<()> {
        writeln!(writer, "annotation_a\tannotation_b\tjaccard\toverlap")?;
        for &(i, j) in &self.co_occurring {
            let Some((jaccard, overlap)) = self.pair_metrics(i, j, similarity_threshold, min_reads)
            else {
                continue;
            };
            let (a, b) = self.resolved_pair(interner, i, j);
            writeln!(writer, "{a}\t{b}\t{jaccard:.2}\t{overlap:.2}")?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::intern::Interner;

    #[test]
    fn test_jaccard_known_values() {
        let mut int = Interner::new();
        let gene_a = int.intern("GENE_A");
        let gene_b = int.intern("GENE_B");

        let mut m = SparseMatrix::new();
        // read 0 → GENE_A, GENE_B
        m.insert_group(0, [gene_a, gene_b].into_iter());
        // read 1 → GENE_A only
        m.insert_group(1, [gene_a].into_iter());
        // read 2 → GENE_B only
        m.insert_group(2, [gene_b].into_iter());

        // GENE_A has reads {0, 1}, GENE_B has reads {0, 2}
        // intersection = {0}, union = {0, 1, 2}
        // Jaccard = 1/3
        let tuples = m.jaccard_tuples(&int, 0.0, 0);
        assert_eq!(tuples.len(), 1);
        assert_eq!(tuples[0].0, "GENE_A");
        assert_eq!(tuples[0].1, "GENE_B");
        assert!((tuples[0].2 - 1.0 / 3.0).abs() < 1e-10);
        assert!((tuples[0].3 - 0.5).abs() < 1e-10); // overlap = 1/2
    }

    #[test]
    fn test_jaccard_identical_sets() {
        let mut int = Interner::new();
        let a = int.intern("A");
        let b = int.intern("B");

        let mut m = SparseMatrix::new();
        m.insert_group(0, [a, b].into_iter());
        m.insert_group(1, [a, b].into_iter());

        let tuples = m.jaccard_tuples(&int, 0.0, 0);
        assert_eq!(tuples.len(), 1);
        assert!((tuples[0].2 - 1.0).abs() < 1e-10);
        assert!((tuples[0].3 - 1.0).abs() < 1e-10); // overlap = 1.0
    }

    #[test]
    fn test_jaccard_disjoint() {
        let mut int = Interner::new();
        let a = int.intern("A");
        let b = int.intern("B");

        let mut m = SparseMatrix::new();
        m.insert_group(0, [a].into_iter());
        m.insert_group(1, [b].into_iter());

        let tuples = m.jaccard_tuples(&int, 0.0, 0);
        assert!(tuples.is_empty()); // no co-occurrence
    }

    #[test]
    fn test_single_annotation() {
        let mut int = Interner::new();
        let a = int.intern("A");

        let mut m = SparseMatrix::new();
        m.insert_group(0, [a].into_iter());
        m.insert_group(1, [a].into_iter());

        let tuples = m.jaccard_tuples(&int, 0.0, 0);
        assert!(tuples.is_empty()); // only one annotation, no pairs
    }

    #[test]
    fn test_empty_matrix() {
        let int = Interner::new();
        let m = SparseMatrix::new();
        let tuples = m.jaccard_tuples(&int, 0.0, 0);
        assert!(tuples.is_empty());
    }

    #[test]
    fn test_min_reads_filters_low_coverage() {
        // GENE_A has only 1 supporting read; GENE_B has 5; their pair should
        // be dropped when min_reads >= 2.
        let mut int = Interner::new();
        let a = int.intern("A");
        let b = int.intern("B");
        let mut m = SparseMatrix::new();
        m.insert_group(0, [a, b].into_iter());
        for g in 1..5 {
            m.insert_group(g, [b].into_iter());
        }
        // min_reads = 0 → keep
        assert_eq!(m.jaccard_tuples(&int, 0.0, 0).len(), 1);
        // min_reads = 2 → A has only 1 supporting group → drop
        assert_eq!(m.jaccard_tuples(&int, 0.0, 2).len(), 0);
    }

    #[test]
    fn test_streaming_jaccard_matches_collected_content() {
        // Both code paths should emit the same set of pairs (modulo order).
        let mut int = Interner::new();
        let a = int.intern("A");
        let b = int.intern("B");
        let c = int.intern("C");
        let mut m = SparseMatrix::new();
        m.insert_group(0, [a, b].into_iter());
        m.insert_group(1, [a, b, c].into_iter());
        m.insert_group(2, [b, c].into_iter());

        let collected = m.jaccard_tuples(&int, 0.0, 0);

        let mut buf = Vec::new();
        m.write_jaccard_streaming(&int, 0.0, 0, &mut buf).unwrap();
        let streamed = String::from_utf8(buf).unwrap();
        let mut streamed_lines: Vec<&str> = streamed.lines().skip(1).collect(); // skip header
        streamed_lines.sort_unstable();

        let mut collected_lines: Vec<String> = collected
            .iter()
            .map(|(a, b, j, o)| format!("{a}\t{b}\t{j:.2}\t{o:.2}"))
            .collect();
        collected_lines.sort_unstable();

        assert_eq!(
            streamed_lines
                .iter()
                .map(|s| s.to_string())
                .collect::<Vec<_>>(),
            collected_lines
        );
    }
}
