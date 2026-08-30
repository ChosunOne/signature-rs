use std::cmp::Ordering;

/// A canonical pair `(i, j)` with `i < j`. Its decomposition slice is
/// `decomp_start ..` up to the next flat entry's `decomp_start`.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Entry {
    pub(crate) i: u16,
    pub(crate) j: u16,
    pub(crate) decomp_start: u32,
}

/// A contiguous run of canonical pairs with `deg(i) = p`, `deg(j) = q`.
#[derive(Clone, Copy, Debug)]
pub(crate) struct BlockMeta {
    /// Degree of the smaller index of every pair in the block.
    pub(crate) p: u8,
    /// Degree of the larger index of every pair in the block.
    pub(crate) q: u8,
    /// Bracket degree `p + q`: every decomposition of the block
    /// scatters into the result's degree-`target` slice only.
    pub(crate) target: u8,
    /// Entry range `entries[start .. end]`, sorted by `(i, j)`.
    /// The decomposition slice of entry `k` is
    /// `entries[k].decomp_start .. entries[k + 1].decomp_start`
    pub(crate) start: u32,
    pub(crate) end: u32,
}

/// A compact, CSR-style storage for coefficient pairs.
///
/// Only *canonical* pairs `(i, j)` with `i != j` and `deg(i) + deg(j) <= max_degree`
/// are stored in ascending `j` order per row.
#[derive(Clone)]
pub(crate) struct FeasibleDecompositions<U> {
    /// Per basis word: its Lyndon degree. The basis is degree-grouped
    /// so this array is non-decreasing.
    index_degrees: Vec<u8>,
    /// For each degree `d`, the basis index where degree-`d` words
    /// start. `degree_range(target)` slices the result vector into
    /// the per-degree write regions used by the commutation kernel.
    degree_starts: Vec<u32>,
    /// Degree blocks, sorted by `(target, p, q)`; entry ranges are
    /// `entries[start .. end]`.
    blocks: Vec<BlockMeta>,
    /// Per block contiguous entries, sorted by `(i, j)` within a
    /// block.
    entries: Vec<Entry>,
    /// Flat decomposition: basis indices *relative to*
    /// `degree_starts[target of the owning block]`.
    decomp_indices: Vec<u32>,
    /// Flat decomposition: basis indices (public `get` path)
    decomp_indices_abs: Vec<u32>,
    /// Flat decomposition: parallel non-zero coefficients
    decomp_coefficients: Vec<U>,
    num_entries: usize,
}

impl<U> FeasibleDecompositions<U> {
    /// Total number of stored feasible pairs
    pub(crate) fn len(&self) -> usize {
        self.num_entries
    }

    /// The maximum basis degree
    pub(crate) fn max_degree(&self) -> usize {
        self.degree_starts.len() - 2
    }

    #[inline]
    pub(crate) fn degree_of(&self, i: usize) -> usize {
        self.index_degrees[i] as usize
    }

    #[inline]
    pub(crate) fn blocks(&self) -> &[BlockMeta] {
        &self.blocks
    }

    /// A block's entry slice, extended by one: `slice[k + 1]` is the
    /// flat successor of `slice[k]`, so entry `k`'s decomposition
    /// range is `slice[k].decomp_start .. slice[k + 1].decomp_start`.
    /// Zip the span with itself shifted one to iterate branch-free.
    #[inline]
    pub(crate) fn entry_span(&self, block: &BlockMeta) -> &[Entry] {
        &self.entries[block.start as usize..=block.end as usize]
    }

    /// The flat relative decomposition index array.
    #[inline]
    pub(crate) fn decomp_indices_rel(&self) -> &[u32] {
        &self.decomp_indices
    }

    /// The flat decomposition coefficient array
    #[inline]
    pub(crate) fn decomp_coeffs(&self) -> &[U] {
        &self.decomp_coefficients
    }

    #[inline]
    pub(crate) fn degree_start(&self, degree: usize) -> usize {
        self.degree_starts[degree] as usize
    }

    /// The result-vector range `[start, end)` holding the degree-`target`
    /// basis words: the disjoint write region of every block with this
    /// target.
    #[inline]
    pub(crate) fn degree_range(&self, target: u8) -> (usize, usize) {
        let t = target as usize;
        (
            self.degree_starts[t] as usize,
            self.degree_starts[t + 1] as usize,
        )
    }

    /// The stored decomposition for the feasible pair `(i, j)`, if any.
    pub(crate) fn get(&self, i: usize, j: usize) -> Option<(&[u32], &[U], bool)> {
        let (min, max, swapped) = match i.cmp(&j) {
            Ordering::Less => (i, j, false),
            Ordering::Greater => (j, i, true),
            Ordering::Equal => return None,
        };

        let (min, max) = (u16::try_from(min).ok()?, u16::try_from(max).ok()?);

        let p = self.index_degrees[min as usize];
        let q = self.index_degrees[max as usize];
        let key = (p as u16 + q as u16, p, q);
        let pos = self
            .blocks
            .partition_point(|b| (b.target as u16, b.p, b.q) < key);
        let block = self.blocks.get(pos)?;
        if (block.p, block.q) != (p, q) {
            return None; // pair is infeasible
        }
        let entries = &self.entries[block.start as usize..block.end as usize];
        let pos = entries.partition_point(|e| (e.i, e.j) < (min, max));
        let entry = entries.get(pos)?;
        if (entry.i, entry.j) != (min, max) {
            return None;
        }

        let from = entry.decomp_start as usize;
        let next = self.entries[block.start as usize + pos + 1].decomp_start as usize;
        Some((
            &self.decomp_indices_abs[from..next],
            &self.decomp_coefficients[from..next],
            swapped,
        ))
    }
}

/// Incremental builder used during `LieSeries::new`: canonical pairs arrive
/// row-major (`i` ascending, `j` ascending within a row) as the
/// structure constants are computed.
pub(crate) struct Builder<U> {
    rows: Vec<Vec<(u16, Vec<u32>, Vec<U>)>>,
    degrees: Vec<u8>,
}

impl<U: Clone> Builder<U> {
    /// `degrees[i]` is the Lyndon degree of basis word `i`; the
    /// basis must be degree-grouped.
    pub(crate) fn new(degrees: &[u8]) -> Self {
        Self {
            rows: (0..degrees.len()).map(|_| Vec::new()).collect(),
            degrees: degrees.to_vec(),
        }
    }

    /// Records the decomposition of the canonical pair `(i, j)`
    /// with `i < j`. Callers must push pairs in row-major order
    /// (`i` ascending, `j > i` ascending within a row) and must
    /// not push infeasible pairs.
    ///
    /// The decomposition must only hit degree-`deg(i) + deg(j)`
    /// basis words.
    pub(crate) fn push(&mut self, i: usize, j: usize, indices: &[usize], coefficients: &[U]) {
        self.rows[i].push((
            u16::try_from(j).expect("basis index exceeds u16"),
            indices
                .iter()
                .map(|&x| u32::try_from(x).expect("basis index exceeds u32"))
                .collect(),
            coefficients.to_vec(),
        ));
    }

    pub(crate) fn finish(self) -> FeasibleDecompositions<U> {
        let degrees = &self.degrees;
        assert!(
            degrees.len() <= u16::MAX as usize,
            "basis size exceeds u16 indices"
        );

        let max_degree: usize = degrees.last().copied().unwrap_or(0) as usize;
        let mut degree_starts: Vec<u32> = vec![0u32; max_degree + 2];
        let mut d: usize = 0usize;
        for (index, &g) in degrees.iter().enumerate() {
            while d < g as usize {
                d += 1;
                degree_starts[d] = index as u32;
            }
        }
        degree_starts[max_degree + 1] = degrees.len() as u32;

        // Flatten and sort all entries by (target, p, q, i, j): blocks in
        // kernel iteration order, entries sorted within their block.
        let mut flat = Vec::new();
        for (i, row) in self.rows.iter().enumerate() {
            for &(ref j, ref indices, ref coefficients) in row {
                let (p, q) = (degrees[i], degrees[*j as usize]);
                debug_assert!(
                    p <= q,
                    "non-degree-grouped basis: i < j but deg(i) > deg(j)"
                );
                flat.push((
                    p + q,
                    p,
                    q,
                    i as u16,
                    *j,
                    indices.clone(),
                    coefficients.clone(),
                ));
            }
        }
        flat.sort_unstable_by(|a, b| (a.0, a.1, a.2, a.3, a.4).cmp(&(b.0, b.1, b.2, b.3, b.4)));

        let total: usize = flat
            .iter()
            .map(|&(_, _, _, _, _, ref idx, _)| idx.len())
            .sum();
        let mut blocks = Vec::<BlockMeta>::new();
        let mut entries = Vec::with_capacity(flat.len());
        let mut decomp_indices = Vec::with_capacity(total);
        let mut decomp_indices_abs = Vec::with_capacity(total);
        let mut decomp_coefficients = Vec::with_capacity(total);

        for (target, p, q, i, j, indices, coefficients) in flat {
            debug_assert!(
                target as usize <= max_degree,
                "infeasible pair ({i}, {j}) pushed: target degree {target} > {max_degree}"
            );
            let target_start: u32 = degree_starts[target as usize];
            if blocks
                .last()
                .is_none_or(|b| (b.target, b.p, b.q) != (target, p, q))
            {
                blocks.push(BlockMeta {
                    p,
                    q,
                    target,
                    start: entries.len() as u32,
                    end: 0,
                });
            }
            debug_assert!(
                indices.iter().all(|&x| degrees[x as usize] == target),
                "decomposition of pair ({i}, {j}) is not degree-homogeneous"
            );
            entries.push(Entry {
                i,
                j,
                decomp_start: decomp_indices.len() as u32,
            });
            for &x in &indices {
                debug_assert!(
                    x >= target_start,
                    "decomposition index below target-degree slice start"
                );
                decomp_indices.push(x - target_start);
            }
            decomp_indices_abs.extend_from_slice(&indices);
            decomp_coefficients.extend_from_slice(&coefficients);
        }
        // Close the blocks' entry ranges (the sentinel entry appended below
        // bounds the last real entry's decomposition slice).
        let mut end: u32 = entries.len() as u32;
        for block in blocks.iter_mut().rev() {
            block.end = end;
            end = block.start;
        }
        let total: u32 = decomp_indices.len() as u32;

        let num_entries: usize = entries.len();
        entries.push(Entry {
            i: 0,
            j: 0,
            decomp_start: total,
        });

        FeasibleDecompositions {
            index_degrees: degrees.clone(),
            degree_starts,
            blocks,
            entries,
            decomp_indices,
            decomp_indices_abs,
            decomp_coefficients,
            num_entries,
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    impl<U> FeasibleDecompositions<U> {
        pub(crate) fn block_iter<'a>(
            &'a self,
            block: &'a BlockMeta,
        ) -> impl Iterator<Item = (&'a Entry, &'a [u32], &'a [U])> + 'a {
            let span = self.entry_span(block);
            span[..span.len() - 1]
                .iter()
                .zip(span[1..].iter())
                .map(move |(e, next)| {
                    let from = e.decomp_start as usize;
                    let to = next.decomp_start as usize;
                    (
                        e,
                        &self.decomp_indices[from..to],
                        &self.decomp_coefficients[from..to],
                    )
                })
        }

        pub(crate) fn iter_entries(
            &self,
        ) -> impl Iterator<Item = (usize, usize, &[u32], &[U])> + '_ {
            self.blocks.iter().flat_map(move |block| {
                let span = self.entry_span(block);
                span[..span.len() - 1]
                    .iter()
                    .zip(span[1..].iter())
                    .map(move |(e, next)| {
                        let from = e.decomp_start as usize;
                        let to = next.decomp_start as usize;
                        (
                            e.i as usize,
                            e.j as usize,
                            &self.decomp_indices_abs[from..to],
                            &self.decomp_coefficients[from..to],
                        )
                    })
            })
        }
    }

    #[test]
    fn builder_round_trip() {
        let mut b = Builder::<f64>::new(&[1, 1, 2, 2, 3, 3]);
        b.push(0, 1, &[3, 2], &[1.5, -2.0]);
        b.push(0, 2, &[4], &[0.25]);
        b.push(1, 2, &[], &[]);
        let t = b.finish();

        assert_eq!(t.len(), 3);
        assert_eq!(t.max_degree(), 3);
        assert_eq!(t.degree_range(1), (0, 2));
        assert_eq!(t.degree_range(2), (2, 4));
        assert_eq!(t.degree_range(3), (4, 6));

        let blocks = t
            .blocks
            .iter()
            .map(|b| (b.p, b.q, b.target))
            .collect::<Vec<_>>();
        assert_eq!(blocks, [(1, 1, 2), (1, 2, 3)]);
        let block12 = t.blocks().iter().find(|b| (b.p, b.q) == (1, 2)).unwrap();
        let entries12 = t
            .block_iter(block12)
            .map(|(e, indices, coefficients)| (e.i, e.j, indices, coefficients))
            .collect::<Vec<_>>();
        assert_eq!(entries12.len(), 2);
        assert_eq!((entries12[0].0, entries12[0].1), (0, 2));
        assert_eq!(entries12[0].2, &[0u32]); // relative to degree_start(3) = 4
        assert_eq!(entries12[0].3, &[0.25]);
        assert_eq!((entries12[1].0, entries12[1].1), (1, 2));
        assert_eq!(entries12[1].2, &[] as &[u32]);
        assert_eq!(entries12[1].3, &[] as &[f64]);

        assert_eq!(t.get(0, 1), Some((&[3u32, 2][..], &[1.5, -2.0][..], false)));
        assert_eq!(t.get(1, 0), Some((&[3u32, 2][..], &[1.5, -2.0][..], true)));
        assert_eq!(t.get(2, 1), Some((&[] as &[u32], &[] as &[f64], true)));
        assert_eq!(t.get(0, 0), None);
        assert_eq!(t.get(1, 1), None);
        assert_eq!(t.get(0, 3), None);
    }

    #[test]
    fn empty_table() {
        let t = Builder::<f64>::new(&[1, 1, 2]).finish();
        assert_eq!(t.len(), 0);
        assert_eq!(t.blocks().len(), 0);
        assert_eq!(t.get(0, 1), None);
    }

    #[test]
    fn relative_indices_match_absolute() {
        let mut b: Builder<f64> = Builder::<f64>::new(&[1, 1, 2, 2, 3, 3]);
        // target 2: only indices 2, 3 are degree-2 words.
        b.push(0, 1, &[2, 3], &[1.0, 2.0]);
        // target 3: only indices 4, 5 are degree-3 words.
        b.push(0, 2, &[4, 5], &[-1.0, 1.0]);
        let t: FeasibleDecompositions<f64> = b.finish();
        for block in t.blocks() {
            let (start, _) = t.degree_range(block.target);
            for (e, rel, _) in t.block_iter(block) {
                let (abs, _, swapped) = t.get(e.i as usize, e.j as usize).expect("stored");
                assert!(!swapped);
                for (rel, abs) in rel.iter().zip(abs) {
                    assert_eq!(start + *rel as usize, *abs as usize);
                }
            }
        }
    }
}
