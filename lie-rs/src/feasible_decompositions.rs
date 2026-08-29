//! Feasible-pair decomposition table: a compact, CSR-style replacement for
//! the O(D²) `commutator_basis_map` / `commutator_basis_map_coefficients`
//! tables.
//!
//! Only *canonical* pairs `(i, j)` with `i < j` and
//! `deg(i) + deg(j) <= max_degree` are stored. Entries are grouped into
//! *anagram units*: all pairs whose brackets land in the same letter
//! multiset are contiguous, sorted by `(target degree, content class)`.
//! Two different units never scatter onto the same basis word, so units
//! are conflict-free write regions for a parallel kernel.
//!
//! Decomposition indices are stored relative to the start of the
//! target-degree slice, with an absolute copy for [`FeasibleDecompositions::get`].

/// A canonical pair `(i, j)` with `i < j`. Its decomposition slice is
/// `decomp_start ..` up to the next flat entry's `decomp_start`.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Entry {
    pub(crate) i: u32,
    pub(crate) j: u32,
    pub(crate) decomp_start: u32,
}

/// A contiguous run of canonical pairs whose brackets land in the same
/// letter multiset: the unit's decompositions only hit degree-`target`
/// words of that content, so two units never write the same basis word.
#[derive(Clone, Debug)]
pub(crate) struct UnitMeta {
    /// Bracket degree `p + q` (equal to `|gamma|`).
    pub(crate) target: u8,
    /// The unit's letter multiset (counts per alphabet letter). Two units
    /// with the same `target` differ in `gamma` and vice versa.
    pub(crate) gamma: Vec<u8>,
    /// Bit `p` set iff the unit contains pairs whose smaller index has
    /// degree `p`.
    pub(crate) p_mask: [u64; 2],
    /// Entry range `entries[start .. end]`, sorted by `(p, q, i, j)`.
    /// The decomposition slice of entry `k` is
    /// `entries[k].decomp_start .. entries[k + 1].decomp_start`.
    pub(crate) start: u32,
    pub(crate) end: u32,
}

/// Coefficient type is generic; decomposition indices are stored as `u32`.
#[derive(Clone)]
pub(crate) struct FeasibleDecompositions<U> {
    /// Per basis word: its Lyndon degree. The basis is degree-grouped, so
    /// this array is non-decreasing.
    index_degrees: Vec<u8>,
    /// Per basis word: its letter content (counts per alphabet letter).
    index_contents: Vec<Vec<u8>>,
    /// For each degree `d`, the basis index where degree-`d` words start.
    /// `degree_range(target)` slices the result vector into the per-degree
    /// write regions used by the commutation kernel.
    degree_starts: Vec<u32>,
    /// Anagram units, sorted by `(target, content class)`; entry ranges
    /// are `entries[start .. end]`.
    units: Vec<UnitMeta>,
    /// Per-block contiguous entries, sorted by `(i, j)` within a block.
    entries: Vec<Entry>,
    /// Flat decomposition: basis indices relative to the owning unit's
    /// target-degree slice (kernel scatter path).
    decomp_indices: Vec<u32>,
    /// Flat decomposition: absolute basis indices (`get` path).
    decomp_indices_abs: Vec<u32>,
    /// Flat decomposition: parallel non-zero coefficients.
    decomp_coefficients: Vec<U>,
    /// Number of real entries (the trailing sentinel entry is excluded).
    num_entries: usize,
}

/// Incremental builder used during `LieSeries::new`: canonical pairs arrive
/// row-major (`i` ascending, `j > i` ascending within a row) as the structure
/// constants are computed.
pub(crate) struct Builder<U> {
    /// Per row `i`: `(j, decomposition indices, decomposition coefficients)`.
    rows: Vec<Vec<(u32, Vec<u32>, Vec<U>)>>,
    /// Letter content of each basis word (counts per alphabet letter), in
    /// basis order.
    contents: Vec<Vec<u8>>,
    /// Degrees of the basis words (content sums), in basis order
    /// (non-decreasing).
    degrees: Vec<u8>,
}

impl<U: Clone> Builder<U> {
    /// `contents[i]` is the letter multiset of basis word `i`; the basis
    /// must be degree-grouped (degrees non-decreasing along basis order).
    pub(crate) fn new(contents: &[Vec<u8>]) -> Self {
        let degrees: Vec<u8> = contents
            .iter()
            .map(|c| c.iter().copied().sum::<u8>())
            .collect();
        Self {
            rows: (0..contents.len()).map(|_| Vec::new()).collect(),
            contents: contents.to_vec(),
            degrees,
        }
    }

    /// Records the decomposition of the canonical pair `(i, j)`. Pairs
    /// must arrive row-major; decompositions must only hit
    /// degree-`deg(i) + deg(j)` basis words.
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
        let degrees = &self.degrees;
        assert!(
            degrees.len() <= u32::MAX as usize,
            "basis size exceeds u32 indices"
        );

        // Degree starts: for each degree, the first basis index carrying it.
        // Absent degrees point at the next present degree (their slice is
        // empty). The basis is degree-grouped, so a single forward scan
        // suffices.
        let max_degree = degrees.last().copied().unwrap_or(0) as usize;
        let mut degree_starts = vec![0u32; max_degree + 2];
        let mut d = 0usize;
        for (index, &g) in degrees.iter().enumerate() {
            while d < g as usize {
                d += 1;
                degree_starts[d] = index as u32;
            }
        }
        degree_starts[max_degree + 1] = degrees.len() as u32;

        // Flatten and sort by (target, content class, p, q, i, j): anagram
        // units in kernel iteration order, entries sorted within their
        // unit.
        let contents = &self.contents;
        let mut flat: Vec<(u8, Vec<u8>, u8, u8, u32, u32, Vec<u32>, Vec<U>)> = Vec::new();
        for (i, row) in self.rows.iter().enumerate() {
            for (j, indices, coefficients) in row {
                let (p, q) = (degrees[i], degrees[*j as usize]);
                debug_assert!(p <= q, "non-degree-grouped basis: i < j but deg(i) > deg(j)");
                let mut gamma = contents[i].clone();
                for k in 0..gamma.len() {
                    gamma[k] += contents[*j as usize][k];
                }
                flat.push((
                    p + q,
                    gamma,
                    p,
                    q,
                    i as u32,
                    *j,
                    indices.clone(),
                    coefficients.clone(),
                ));
            }
        }
        flat.sort_unstable_by(|a, b| {
            (a.0, &a.1, a.2, a.3, a.4, a.5).cmp(&(b.0, &b.1, b.2, b.3, b.4, b.5))
        });

        let total: usize = flat.iter().map(|(_, _, _, _, _, _, idx, _)| idx.len()).sum();
        let mut units: Vec<UnitMeta> = Vec::new();
        let mut entries = Vec::with_capacity(flat.len());
        let mut decomp_indices = Vec::with_capacity(total);
        let mut decomp_indices_abs = Vec::with_capacity(total);
        let mut decomp_coefficients = Vec::with_capacity(total);

        for (target, gamma, p, _q, i, j, indices, coefficients) in flat {
            debug_assert!(
                target as usize <= max_degree,
                "infeasible pair ({i}, {j}) pushed: target degree {target} > {max_degree}"
            );
            let target_start = degree_starts[target as usize];
            if units
                .last()
                .is_none_or(|u| (u.target, &u.gamma) != (target, &gamma))
            {
                units.push(UnitMeta {
                    target,
                    gamma: gamma.clone(),
                    p_mask: [0, 0],
                    start: entries.len() as u32,
                    end: 0,
                });
            }
            let unit = units.last_mut().unwrap();
            let mask_word = &mut unit.p_mask[(p / 64) as usize];
            *mask_word |= 1u64 << (p % 64);
            debug_assert!(
                indices.iter().all(|&x| contents[x as usize] == gamma),
                "decomposition of pair ({i}, {j}) is not content-homogeneous"
            );
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
        // Close the units' entry ranges (the sentinel entry appended below
        // bounds the last real entry's decomposition slice).
        let mut end = entries.len() as u32;
        for unit in units.iter_mut().rev() {
            unit.end = end;
            end = unit.start;
        }
        let total = decomp_indices.len() as u32;

        let num_entries = entries.len();
        entries.push(Entry {
            i: 0,
            j: 0,
            decomp_start: total,
        });

        FeasibleDecompositions {
            index_degrees: degrees.clone(),
            index_contents: contents.clone(),
            degree_starts,
            units,
            entries,
            decomp_indices,
            decomp_indices_abs,
            decomp_coefficients,
            num_entries,
        }
    }
}

impl<U> FeasibleDecompositions<U> {
    /// Total number of stored feasible pairs.
    pub(crate) fn len(&self) -> usize {
        self.num_entries
    }

    /// The maximum basis degree.
    #[cfg_attr(not(test), allow(dead_code))]
    pub(crate) fn max_degree(&self) -> usize {
        self.degree_starts.len() - 2
    }

    /// The Lyndon degree of basis word `i` (O(1)).
    #[inline]
    pub(crate) fn degree_of(&self, i: usize) -> usize {
        self.index_degrees[i] as usize
    }

    /// Anagram units in kernel iteration order.
    #[inline]
    pub(crate) fn units(&self) -> &[UnitMeta] {
        &self.units
    }

    /// The flat entry array in storage order.
    #[inline]
    pub(crate) fn entries(&self) -> &[Entry] {
        &self.entries
    }

    /// A unit's entry slice, extended by one: `slice[k + 1]` is the flat
    /// successor of `slice[k]`, so entry `k`'s decomposition range is
    /// `slice[k].decomp_start .. slice[k + 1].decomp_start`.
    #[cfg_attr(not(test), allow(dead_code))]
    #[inline]
    pub(crate) fn entry_span(&self, unit: &UnitMeta) -> &[Entry] {
        &self.entries[unit.start as usize..=unit.end as usize]
    }

    /// The flat relative decomposition index array (kernel scatter path).
    #[inline]
    pub(crate) fn decomp_indices_rel(&self) -> &[u32] {
        &self.decomp_indices
    }

    /// The flat decomposition coefficient array.
    #[inline]
    pub(crate) fn decomp_coeffs(&self) -> &[U] {
        &self.decomp_coefficients
    }

    /// A unit's entries as `(entry, decomposition indices, decomposition
    /// coefficients)`. Decomposition indices are relative to the start of
    /// the unit's target-degree slice.
    #[cfg(test)]
    pub(crate) fn unit_iter<'a>(
        &'a self,
        unit: &'a UnitMeta,
    ) -> impl Iterator<Item = (&'a Entry, &'a [u32], &'a [U])> + 'a {
        let span = self.entry_span(unit);
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

    /// The basis index where degree-`degree` words start (see
    /// [`FeasibleDecompositions::degree_range`]).
    #[inline]
    pub(crate) fn degree_start(&self, degree: usize) -> usize {
        self.degree_starts[degree] as usize
    }

    /// The result-vector range `[start, end)` holding the degree-`target`
    /// basis words.
    #[inline]
    pub(crate) fn degree_range(&self, target: u8) -> (usize, usize) {
        let t = target as usize;
        (
            self.degree_starts[t] as usize,
            self.degree_starts[t + 1] as usize,
        )
    }

    /// All entries in storage order as
    /// `(i, j, absolute decomposition indices, decomposition coefficients)`.
    #[cfg(test)]
    pub(crate) fn iter_entries(&self) -> impl Iterator<Item = (usize, usize, &[u32], &[U])> + '_ {
        self.units.iter().flat_map(move |unit| {
            let span = self.entry_span(unit);
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

    /// The decomposition of the bracket `[basis[i], basis[j]]`, if stored.
    ///
    /// Querying with `j < i` returns the canonical data with
    /// `swapped = true`; the caller must negate the coefficients. `i == j`
    /// is never stored.
    pub(crate) fn get(&self, i: usize, j: usize) -> Option<(&[u32], &[U], bool)> {
        let (min, max, swapped) = match i.cmp(&j) {
            std::cmp::Ordering::Less => (i, j, false),
            std::cmp::Ordering::Greater => (j, i, true),
            std::cmp::Ordering::Equal => return None,
        };
        let (min, max) = (
            u32::try_from(min).ok()?,
            u32::try_from(max).ok()?,
        );
        let (c_min, c_max) = (
            &self.index_contents[min as usize],
            &self.index_contents[max as usize],
        );
        let mut gamma = c_min.clone();
        for (g, c) in gamma.iter_mut().zip(c_max) {
            *g += c;
        }
        let target = gamma.iter().copied().sum::<u8>();
        // Units are sorted by (target, content class) and the class
        // determines the target, so the unit is found with one binary
        // search. The target is widened: arbitrary (infeasible) queries may
        // overflow `u8`.
        let pos = self
            .units
            .partition_point(|u| (u.target as u16, &u.gamma) < (target as u16, &gamma));
        let unit = self.units.get(pos)?;
        if unit.gamma != gamma {
            return None; // no unit for this content pair: infeasible pair
        }
        let entries = &self.entries[unit.start as usize..unit.end as usize];
        let pos = entries.partition_point(|e| (e.i, e.j) < (min, max));
        let entry = entries.get(pos)?;
        if (entry.i, entry.j) != (min, max) {
            return None;
        }
        let from = entry.decomp_start as usize;
        // The flat successor bounds the slice (next unit's first entry or
        // the trailing sentinel; both carry the right `decomp_start`).
        let next = self.entries[unit.start as usize + pos + 1].decomp_start as usize;
        Some((
            &self.decomp_indices_abs[from..next],
            &self.decomp_coefficients[from..next],
            swapped,
        ))
    }
}

#[cfg(test)]
mod test {
    use super::*;

    /// Five fake basis words over a 2-letter alphabet: `a`, `b` (degree 1),
    /// `ab` (degree 2), `aab`, `abb` (degree 3). Degree-`t` decompositions
    /// may only hit degree-`t` words whose content is the multiset union of
    /// the pair's contents — e.g. `[a, b]` -> `ab`, `[a, ab]` -> `aab`.
    #[test]
    fn builder_round_trip() {
        let contents = vec![
            vec![1, 0], // 0: a
            vec![0, 1], // 1: b
            vec![1, 1], // 2: ab
            vec![2, 1], // 3: aab
            vec![1, 2], // 4: abb
        ];
        let mut b = Builder::<f64>::new(&contents);
        // Row-major pushes of canonical pairs, including an empty
        // decomposition and a row with no feasible pairs.
        // (0,1): a+b, target 2, class (1,1,0) -> word 2 (ab)
        b.push(0, 1, &[2], &[1.5]);
        // (0,2): a+ab, target 3, class (2,1,0) -> word 3 (aab)
        b.push(0, 2, &[3], &[0.25]);
        // (1,2): b+ab, target 3, class (1,2,0): word 4 (abb) exists but the
        // pair's decomposition is empty here; the unit still exists.
        b.push(1, 2, &[], &[]);
        let t = b.finish();

        assert_eq!(t.len(), 3);
        assert_eq!(t.max_degree(), 3);
        assert_eq!(t.degree_range(1), (0, 2));
        assert_eq!(t.degree_range(2), (2, 3));
        assert_eq!(t.degree_range(3), (3, 5));

        // Units sorted by (target, content class), entries sorted by (i, j).
        let units: Vec<_> = t
            .units()
            .iter()
            .map(|u| (u.target, u.gamma.clone()))
            .collect();
        assert_eq!(
            units,
            [(2, vec![1, 1]), (3, vec![1, 2]), (3, vec![2, 1])]
        );
        let unit12: &UnitMeta = t
            .units()
            .iter()
            .find(|u| u.target == 3 && u.gamma[0] == 1)
            .unwrap();
        let entries12: Vec<_> = t
            .unit_iter(unit12)
            .map(|(e, indices, coefficients)| (e.i, e.j, indices, coefficients))
            .collect();
        assert_eq!(entries12.len(), 1);
        assert_eq!((entries12[0].0, entries12[0].1), (1, 2));
        assert_eq!(entries12[0].2, &[] as &[u32]);
        assert_eq!(entries12[0].3, &[] as &[f64]);

        // Canonical lookups, including the mirrored (negated) query.
        assert_eq!(t.get(0, 1), Some((&[2u32][..], &[1.5][..], false)));
        assert_eq!(t.get(1, 0), Some((&[2u32][..], &[1.5][..], true)));
        assert_eq!(t.get(2, 1), Some((&[] as &[u32], &[] as &[f64], true)));
        assert_eq!(t.get(0, 0), None);
        assert_eq!(t.get(1, 1), None);
        assert_eq!(t.get(0, 3), None); // feasible degrees, but never pushed
    }

    #[test]
    fn empty_table() {
        let t = Builder::<f64>::new(&[vec![1, 0], vec![0, 1], vec![2, 0]]).finish();
        assert_eq!(t.len(), 0);
        assert_eq!(t.max_degree(), 2);
        assert_eq!(t.units().len(), 0);
        assert_eq!(t.get(0, 1), None);
    }

    /// The kernel scatters with indices relative to the target-degree slice;
    /// the absolute copy used by `get` must reconstruct the same indices.
    #[test]
    fn relative_indices_match_absolute() {
        let contents = vec![
            vec![1, 0], // 0: a
            vec![0, 1], // 1: b
            vec![1, 1], // 2: ab   (degree 2, class (1,1,0))
            vec![2, 1], // 3: aab  (degree 3, class (2,1,0))
        ];
        let mut b = Builder::<f64>::new(&contents);
        // target 2, class (1,1,0): word 2 is the only degree-2 word of that
        // content.
        b.push(0, 1, &[2], &[1.0]);
        // target 3, class (2,1,0): word 3 is the only degree-3 word of that
        // content.
        b.push(0, 2, &[3], &[-1.0]);
        let t = b.finish();

        for unit in t.units() {
            let (start, _) = t.degree_range(unit.target);
            for (e, rel, _) in t.unit_iter(unit) {
                let (abs, _, swapped) = t.get(e.i as usize, e.j as usize).expect("stored");
                assert!(!swapped);
                for (rel, abs) in rel.iter().zip(abs) {
                    assert_eq!(start + *rel as usize, *abs as usize);
                }
            }
        }
    }
}