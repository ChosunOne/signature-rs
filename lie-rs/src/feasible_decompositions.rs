/// A compact, CSR-style storage for coefficient pairs.
///
/// Only pairs `(i, j)` with `i != j` and `deg(i) + deg(j) <= max_degree`
/// are stored in ascending `j` order per row.
#[derive(Clone)]
pub(crate) struct FeasibleDecompositions<U> {
    /// Row `i` occupies entries `row_offsets[i] .. row_offsets[i + 1]`.
    row_offsets: Vec<u32>,
    /// Per entry (in row-major order): the paired basis index `j`.
    pair_j: Vec<u32>,
    /// Entry `e`'s decomposition occupies `decomp_offsets[e] .. decomp_offsets[e + 1]`
    /// of the two flat arrays.
    decomp_offsets: Vec<u32>,
    /// Flat decomposition: basis indices of the decomposition terms.
    decomp_indices: Vec<u32>,
    /// Flat decomposition: parallel non-zero coefficients
    decomp_coefficients: Vec<U>,
}

impl<U> FeasibleDecompositions<U> {
    /// Number of rows
    pub(crate) fn num_rows(&self) -> usize {
        self.row_offsets.len() - 1
    }

    /// Total number of stored feasible pairs
    pub(crate) fn len(&self) -> usize {
        self.pair_j.len()
    }

    /// Iterate row `i`'s entries as `(j, decomposition indices, decomposition coefficients)`
    /// in ascending `j` order.
    pub(crate) fn row(&self, i: usize) -> impl Iterator<Item = (usize, &[u32], &[U])> + '_ {
        let start = self.row_offsets[i] as usize;
        let end = self.row_offsets[i + 1] as usize;
        let pair_j = &self.pair_j[start..end];
        let offsets = &self.decomp_offsets[start..=end];
        let indices = &self.decomp_indices;
        let coefficients = &self.decomp_coefficients;
        pair_j.iter().enumerate().map(move |(k, &j)| {
            let from = offsets[k] as usize;
            let to = offsets[k + 1] as usize;
            (j as usize, &indices[from..to], &coefficients[from..to])
        })
    }

    /// The stored decomposition for the feasible pair `(i, j)`, if any.
    pub(crate) fn get(&self, i: usize, j: usize) -> Option<(&[u32], &[U])> {
        let start = self.row_offsets[i] as usize;
        let end = self.row_offsets[i + 1] as usize;
        let row = &self.pair_j[start..end];
        let pos = start + row.binary_search(&(j as u32)).ok()?;
        let from = self.decomp_offsets[pos] as usize;
        let to = self.decomp_offsets[pos + 1] as usize;
        Some((
            &self.decomp_indices[from..to],
            &self.decomp_coefficients[from..to],
        ))
    }
}

/// Incremental builder used during `LieSeries::new`: pairs arrive
/// row-major (`i` ascending, `j` ascending within a row) as the
/// structure constants are computed.
pub(crate) struct Builder<U> {
    rows: Vec<Vec<(u32, Vec<u32>, Vec<U>)>>,
}

impl<U: Clone> Builder<U> {
    pub(crate) fn new(num_rows: usize) -> Self {
        Self {
            rows: (0..num_rows).map(|_| Vec::new()).collect(),
        }
    }

    pub(crate) fn push(&mut self, i: usize, j: usize, indices: &[usize], coefficients: &[U]) {
        self.rows[i].push((
            u32::try_from(j).expect("basis index exceeds u32"),
            indices
                .iter()
                .map(|&x| u32::try_from(x).expect("basis index exceeds u32"))
                .collect(),
            coefficients.to_vec(),
        ));
    }

    pub(crate) fn finish(self) -> FeasibleDecompositions<U> {
        let num_rows = self.rows.len();
        let entries = self.rows.iter().map(|r| r.len()).sum::<usize>();
        let total = self
            .rows
            .iter()
            .flat_map(|r| r.iter())
            .map(|(_, idx, _)| idx.len())
            .sum::<usize>();

        let mut row_offsets = Vec::with_capacity(num_rows + 1);
        let mut pair_j = Vec::with_capacity(entries);
        let mut decomp_offsets = Vec::with_capacity(entries + 1);
        let mut decomp_indices = Vec::with_capacity(total);
        let mut decomp_coefficients = Vec::with_capacity(total);
        decomp_offsets.push(0u32);

        let mut entry_cursor = 0u32;
        let mut decomp_cursor = 0u32;
        for row in &self.rows {
            row_offsets.push(entry_cursor);
            for (j, indices, coefficients) in row {
                pair_j.push(*j);
                decomp_indices.extend_from_slice(indices);
                decomp_coefficients.extend_from_slice(coefficients);
                decomp_cursor += indices.len() as u32;
                decomp_offsets.push(decomp_cursor);
                entry_cursor += 1;
            }
        }
        row_offsets.push(entry_cursor);

        FeasibleDecompositions {
            row_offsets,
            pair_j,
            decomp_offsets,
            decomp_indices,
            decomp_coefficients,
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn builder_round_trip() {
        let mut b = Builder::<f64>::new(3);
        b.push(0, 1, &[2, 0], &[1.5, -2.0]);
        b.push(0, 2, &[1], &[0.25]);
        b.push(2, 0, &[], &[]);
        let t = b.finish();

        assert_eq!(t.num_rows(), 3);
        assert_eq!(t.len(), 3);

        let row0 = t.row(0).collect::<Vec<_>>();
        assert_eq!(row0.len(), 2);
        assert_eq!(row0[0].0, 1);
        assert_eq!(row0[0].1, &[2u32, 0]);
        assert_eq!(row0[0].2, &[1.5, -2.0]);
        assert_eq!(row0[1].0, 2);
        assert_eq!(row0[1].1, &[1u32]);

        assert_eq!(t.row(1).count(), 0);

        let row2 = t.row(2).collect::<Vec<_>>();
        assert_eq!(row2.len(), 1);
        assert_eq!(row2[0].0, 0);
        assert_eq!(row2[0].1, &[]);
        assert_eq!(row2[0].2, &[]);
    }

    #[test]
    fn empty_table() {
        let t = Builder::<f64>::new(4).finish();
        assert_eq!(t.num_rows(), 4);
        assert_eq!(t.len(), 0);
        for i in 0..4 {
            assert_eq!(t.row(i).count(), 0);
        }
    }
}
