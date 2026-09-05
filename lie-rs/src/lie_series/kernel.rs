use super::*;
use super::gating::KernelGating;

/// Receiver for the absolute basis indices a kernel call scatters onto.
/// The collecting implementation deduplicates through a dirty bitset.
pub(crate) trait ScatterSink {
    #[inline(always)]
    fn scatter(&mut self, _index: usize) {}
}

pub(crate) struct NoSink;

impl ScatterSink for NoSink {}

pub(crate) struct CollectSink<'a> {
    dirty: &'a mut [u64],
    out: &'a mut Vec<usize>,
    cutoff: usize,
}

impl ScatterSink for CollectSink<'_> {
    #[inline(always)]
    fn scatter(&mut self, index: usize) {
        if index >= self.cutoff {
            return;
        }
        let word = &mut self.dirty[index / 64];
        let bit = 1u64 << (index % 64);
        if *word & bit == 0 {
            *word |= bit;
            self.out.push(index);
        }
    }
}

impl<
    T: Clone + Ord + Generator + Hash + Eq,
    U: Clone
        + Default
        + One
        + Zero
        + Eq
        + MulAssign
        + Neg<Output = U>
        + Hash
        + AddAssign
        + 'static,
> LieSeries<T, U>
{
    /// Indices of non-zero coefficients that the commutation kernel will
    /// actually use: non-zero values on basis elements below `max_degree`.
    pub fn nonzero_coefficient_indices(&self, coefficients: &[U]) -> Vec<usize> {
        let cutoff = self
            .feasible_decompositions
            .degree_start(self.max_degree)
            .min(coefficients.len());
        coefficients[..cutoff]
            .iter()
            .enumerate()
            .filter(|(_, c)| !c.is_zero())
            .map(|x| x.0)
            .collect()
    }

    /// Total number of stored feasible decomposition pairs of this series'
    /// table — the static per-fold sweep-work upper bound (a steady fold's
    /// gating resolves most of the table's entries once operand supports
    /// are dense).
    pub fn feasible_table_len(&self) -> usize {
        self.feasible_decompositions.len()
    }

    /// Allocless [`Self::nonzero_coefficient_indices`] equality test:
    /// `true` iff `coefficients`' kernel-visible nonzero support equals
    /// `support` (both understood over `[0, cutoff)`, `support` sorted
    /// ascending — exactly what the reference call returns). Walks the
    /// coefficients once, matching nonzero positions against `support`
    /// in order; early-exits on the first mismatch. Same result as
    /// `self.nonzero_coefficient_indices(coefficients) == support`
    /// without the per-candidate `Vec` allocation.
    pub fn has_support(&self, coefficients: &[U], support: &[usize]) -> bool {
        let cutoff = self
            .feasible_decompositions
            .degree_start(self.max_degree)
            .min(coefficients.len());
        let mut remaining = support.iter().copied();
        let mut expected = remaining.next();
        for (i, c) in coefficients[..cutoff].iter().enumerate() {
            if !c.is_zero() {
                match expected {
                    Some(j) if j == i => expected = remaining.next(),
                    _ => return false,
                }
            }
        }
        expected.is_none()
    }

    #[cfg_attr(
        feature = "tracing",
        tracing::instrument(
            level = "debug",
            name = "lie_series_commutator_coefficients",
            skip_all,
            fields(
                basis_len = a_series.basis.len(),
                max_degree = a_series.max_degree
            )
        )
    )]
    pub fn commutator_coefficients(
        a_series: &LieSeries<T, U>,
        a_coefficients: &[U],
        b_coefficients: &[U],
        result_coefficients: &mut [U],
    ) where
        U: Send + Sync,
    {
        let a_nonzero = a_series.nonzero_coefficient_indices(a_coefficients);
        let b_nonzero = a_series.nonzero_coefficient_indices(b_coefficients);

        LieSeries::commutator_coefficients_with_nonzero(
            a_series,
            a_coefficients,
            &a_nonzero,
            b_coefficients,
            &b_nonzero,
            result_coefficients,
        );
    }

    pub fn commutator_coefficients_with_nonzero_collecting(
        a_series: &LieSeries<T, U>,
        a_coefficients: &[U],
        a_nonzero: &[usize],
        b_coefficients: &[U],
        b_nonzero: &[usize],
        result_coefficients: &mut [U],
        dirty: &mut [u64],
        targets: &mut Vec<usize>,
    ) {
        for word in dirty.iter_mut() {
            *word = 0;
        }
        targets.clear();
        let cutoff = a_series
            .feasible_decompositions
            .degree_start(a_series.max_degree);
        let gating = Self::kernel_prologue(a_series, a_nonzero, b_nonzero);
        let table = &a_series.feasible_decompositions;
        // The gating's word classes carry public positions: one single-writer
        // store per word, straight into the public result buffer (the old
        // class-scratch + permutation detour densified a scatter of repeated
        // `+=`s — one store per word needs no densifying). The sink reports
        // public basis indices directly.
        Self::sweep_words_serial(
            a_series,
            a_coefficients,
            b_coefficients,
            &gating,
            result_coefficients,
            table.entries(),
            table.decomp_indices_rel(),
            &mut CollectSink {
                dirty,
                out: targets,
                cutoff,
            },
        );
    }

    /// Class-contiguous variant of
    /// [`Self::commutator_coefficients_with_nonzero_collecting`]: operands
    /// are class-ordered, support lists class-indexed, and the sink reports
    /// class-indexed positions. No scratch, no permutation.
    #[allow(clippy::too_many_arguments)]
    /// Handle to the compiled feasible-decomposition table (basis-derived
    /// integer data — no `T`). `#[doc(hidden)]` kernel-plumbing accessor
    /// for the DAG crate's level-parallel collecting rebuild, which must
    /// capture only Send+Sync, T-free state (see `class_collect_kernel`).
    #[doc(hidden)]
    pub fn feasible_decompositions_handle(&self) -> &Arc<FeasibleDecompositions<U>> {
        &self.feasible_decompositions
    }

    pub(crate) fn commutator_coefficients_class_with_nonzero_collecting(
        a_series: &LieSeries<T, U>,
        order: &ClassOrder,
        a_coefficients: &[U],
        a_nonzero_cls: &[usize],
        b_coefficients: &[U],
        b_nonzero_cls: &[usize],
        result_coefficients: &mut [U],
        dirty: &mut [u64],
        targets: &mut Vec<usize>,
    ) {
        Self::class_collect_kernel(
            &a_series.feasible_decompositions,
            a_series.basis.len(),
            a_series.max_degree,
            order,
            a_coefficients,
            a_nonzero_cls,
            b_coefficients,
            b_nonzero_cls,
            result_coefficients,
            dirty,
            targets,
        );
    }

    /// T-free core of
    /// [`Self::commutator_coefficients_class_with_nonzero_collecting`]:
    /// everything the class collect kernel touches is basis-derived
    /// integer data (the feasible table, the class order) plus U values —
    /// the letter type never enters the kernel. This is what lets the
    /// collecting rebuild run level-parallel without shipping `T` across
    /// threads (see `CommutatorDag::ensure_lists_steady`).
    #[doc(hidden)]
    pub fn class_collect_kernel(
        table: &FeasibleDecompositions<U>,
        basis_len: usize,
        max_degree: usize,
        order: &ClassOrder,
        a_coefficients: &[U],
        a_nonzero_cls: &[usize],
        b_coefficients: &[U],
        b_nonzero_cls: &[usize],
        result_coefficients: &mut [U],
        dirty: &mut [u64],
        targets: &mut Vec<usize>,
    ) {
        for word in dirty.iter_mut() {
            *word = 0;
        }
        targets.clear();
        let cutoff = table.degree_start(max_degree);
        let gating = Self::kernel_prologue_cached_class_core(
            table,
            basis_len,
            a_nonzero_cls,
            b_nonzero_cls,
            order,
            &mut GatingCache::default(),
        );
        Self::sweep_words_serial_core(
            table.decomp_coeffs(),
            a_coefficients,
            b_coefficients,
            &gating,
            result_coefficients,
            order.entries_cls(),
            order.decomp_cls(),
            &mut CollectSink {
                dirty,
                out: targets,
                cutoff,
            },
        );
    }

    pub fn commutator_coefficients_with_nonzero(
        a_series: &LieSeries<T, U>,
        a_coefficients: &[U],
        a_nonzero: &[usize],
        b_coefficients: &[U],
        b_nonzero: &[usize],
        result_coefficients: &mut [U],
    ) where
        U: Send + Sync,
    {
        let mut job = KernelJob {
            a: a_coefficients.as_ptr(),
            a_len: a_coefficients.len(),
            a_nonzero,
            b: b_coefficients.as_ptr(),
            b_len: b_coefficients.len(),
            b_nonzero,
            result: result_coefficients.as_mut_ptr(),
            result_len: result_coefficients.len(),
            a_shift: IDENTITY_SHIFTS.as_ptr(),
            b_shift: IDENTITY_SHIFTS.as_ptr(),
            r_shift: IDENTITY_SHIFTS.as_ptr(),
        };
        commutator_coefficients_batch(a_series, std::slice::from_mut(&mut job));
        // `job` moved into the slice; its fields are all references/raw
        // pointers — nothing to drop.
    }
    pub(super) fn sweep_words_serial<S: ScatterSink>(
        a_series: &LieSeries<T, U>,
        a_coefficients: &[U],
        b_coefficients: &[U],
        gating: &KernelGating,
        result_coefficients: &mut [U],
        entries: &[Entry],
        decomp_rels: &[u32],
        sink: &mut S,
    ) {
        Self::sweep_words_serial_core(
            a_series.feasible_decompositions.decomp_coeffs(),
            a_coefficients,
            b_coefficients,
            gating,
            result_coefficients,
            entries,
            decomp_rels,
            sink,
        );
    }

    /// T-free core of [`Self::sweep_words_serial`]: everything the sweep
    /// touches is basis-derived integer tables plus U values — the letter
    /// type never enters the kernel (this is what lets the collecting
    /// rebuild run level-parallel; see `CommutatorDag::ensure_lists_steady`).
    fn sweep_words_serial_core<S: ScatterSink>(
        decomp_coefficients: &[U],
        a_coefficients: &[U],
        b_coefficients: &[U],
        gating: &KernelGating,
        result_coefficients: &mut [U],
        entries: &[Entry],
        decomp_rels: &[u32],
        sink: &mut S,
    ) {
        // Single-phase per-unit sweep: each unit's active entries are
        // visited in table order, each term computed exactly once and its
        // row contributions added straight into the result buffer. A
        // word's adds happen inside its one owning unit, in table-entry
        // order — exactly the serial scatter's per-word `+=` sequence —
        // and callers' buffers are accumulated into (their current value
        // starts the sum). Degree-slice starts are identical in both
        // working layouts, so `rs + rel` addresses the working-space
        // result directly.
        for u in gating.units.iter() {
            for tp in u.ticket_start as usize..u.ticket_end as usize {
                let ticket = gating.tickets[tp];
                let e = (ticket & TICKET_INDEX_MASK) as usize;
                let entry = entries[e];
                let p_active = ticket & TICKET_P_ACTIVE != 0;
                let q_active = ticket & TICKET_Q_ACTIVE != 0;
                let (i, j) = (entry.i as usize, entry.j as usize);
                // Orientation (`a = j, b = i` only) negates: `[basis[j],
                // basis[i]]` is the negation of the stored decomposition.
                // `raw_mul`/`raw_add_assign` skip the float wrappers'
                // per-op NaN checks (raw-float fast path, see `raw_mul`'s
                // NaN policy).
                let term = if p_active {
                    let mut t = raw_mul(&a_coefficients[i], &b_coefficients[j]);
                    if q_active {
                        raw_add_assign(&mut t, &-raw_mul(&a_coefficients[j], &b_coefficients[i]));
                    }
                    t
                } else {
                    -raw_mul(&a_coefficients[j], &b_coefficients[i])
                };
                let from = entry.decomp_start as usize;
                let to = entries[e + 1].decomp_start as usize;
                let base = u.rs as usize;
                for (k, &rel) in decomp_rels[from..to].iter().enumerate() {
                    let pos = base + rel as usize;
                    raw_add_assign(
                        &mut result_coefficients[pos],
                        &raw_mul(&decomp_coefficients[from + k], &term),
                    );
                }
            }
        }
        // All units are complete: report every active position (globally
        // ascending — the same sequence the per-word fan-in reported). No
        // unit ever writes another's words.
        for &pos in gating.unit_words.iter() {
            sink.scatter(pos as usize);
        }
    }

    pub fn commutator_assign(&mut self, other: &Self)
    where
        U: Send + Sync,
    {
        let mut coefficients = vec![U::default(); self.coefficients.len()];
        LieSeries::commutator_coefficients(
            self,
            &self.coefficients,
            &other.coefficients,
            &mut coefficients,
        );
        self.coefficients = coefficients;
    }
}

impl<
    T: Clone + Ord + Generator + Hash,
    U: Clone
        + Default
        + One
        + Zero
        + Eq
        + MulAssign
        + Neg<Output = U>
        + Hash
        + AddAssign
        + Send
        + Sync
        + 'static,
> Commutator<&Self> for LieSeries<T, U>
{
    type Output = Self;

    /// Calculates the lie bracket `[A, B]` for a lie series for terms within the commutator basis.
    fn commutator(&self, other: &Self) -> Self::Output {
        let mut coefficients = vec![U::default(); self.coefficients.len()];
        LieSeries::<T, U>::commutator_coefficients(
            self,
            // other,
            &self.coefficients,
            &other.coefficients,
            &mut coefficients,
        );
        self.with_coefficients(coefficients)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ordered_float::NotNan;
    use rstest::rstest;

    #[rstest]
    #[case(2, 2, vec![1, 2, 3], vec![4, 5, 6], vec![0, 0, -3])]
    #[case(2, 2, vec![3, 2, 1], vec![1, 2, 3], vec![0, 0, 4])]
    #[case(2, 3, vec![1, 2, 3, 4, 5], vec![6, 7, 8, 9, 10], vec![0, 0, -5, -10, 5])]
    #[case(3, 3,
        vec![1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        vec![5, 3, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        vec![0, 0, 0, -7, -14, -7, 0, 0, 0, 0, 0, 0, 0, 0])]
    #[case(3, 3,
        vec![1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        vec![0, 0, 0, -7, -14, -7, 0, 0, 0, 0, 0, 0, 0, 0],
        vec![0, 0, 0, 0, 0, 0, -7, -14, 14, 14, 49, 42, -14, 21])]
    #[case(3, 4,
        vec![
            1, 2, 3, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0
        ],
        vec![
            5, 3, 1, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0
        ],
        vec![
            0, 0, 0, -7, -14, -7, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0
        ],
    )]
    #[case(3, 4,
        vec![
            1, 2, 3, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0
        ],
        vec![
            0, 0, 0, -7, -14, -7, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0
        ],
        vec![
            0, 0, 0, 0, 0, 0, -7, -14,
            14, 14, 49, 42, -14, 21, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0,
        ],
    )]
    fn test_lie_series_commutation(
        #[case] num_generators: usize,
        #[case] basis_depth: usize,
        #[case] a_coefficients: Vec<i128>,
        #[case] b_coefficients: Vec<i128>,
        #[case] expected_coefficients: Vec<i128>,
    ) {
        use lyndon_rs::lyndon::{LyndonBasis, Sort};
        use num_rational::Ratio;

        let basis = LyndonBasis::<u8>::new(num_generators, Sort::Lexicographical)
            .generate_basis(basis_depth);
        let a_coefficients = a_coefficients
            .into_iter()
            .map(Ratio::<i128>::from_integer)
            .collect::<Vec<_>>();
        let a = LieSeries::new(basis.clone(), a_coefficients);

        let b_coefficients = b_coefficients
            .into_iter()
            .map(Ratio::<i128>::from_integer)
            .collect::<Vec<_>>();
        let b = LieSeries::new(basis.clone(), b_coefficients);
        let expected_coefficients = expected_coefficients
            .into_iter()
            .map(Ratio::<i128>::from_integer)
            .collect::<Vec<_>>();

        let series = comm![a, b];
        assert_eq!(series.coefficients.len(), expected_coefficients.len());
        dbg!(&series.coefficients);
        for (i, c) in series.coefficients.iter().enumerate() {
            assert_eq!(
                *c, expected_coefficients[i],
                "{i}: {c:?} != {:?}",
                expected_coefficients[i]
            );
        }
    }

    /// The raw-float fast path (`NotNan<f64>` / `f64` dispatch) and the
    /// generic path (`Ratio<i128>`) must agree. Integer-valued coefficients
    /// keep every intermediate exactly representable in `f64` (magnitudes
    /// stay far below 2^53), so the comparison is exact. (Ported from the
    /// original raw-float fast path change.)
    #[test]
    fn raw_float_kernel_matches_rationals() {
        use lyndon_rs::lyndon::{LyndonBasis, Sort};
        use num_rational::Ratio;
        use num_traits::ToPrimitive;

        for (d, m) in [(2usize, 6usize), (3, 5)] {
            let words = LyndonBasis::<u8>::new(d, Sort::Lexicographical).generate_basis(m);
            let coeffs = |salt: usize| {
                (0..words.len())
                    .map(|i| ((i * 7 + salt * 13) % 21) as i128 - 10)
                    .collect::<Vec<_>>()
            };
            let (a_int, b_int) = (coeffs(1), coeffs(2));

            // Raw-float path.
            let a_f = LieSeries::<u8, NotNan<f64>>::new(
                words.clone(),
                a_int
                    .iter()
                    .map(|&x| NotNan::new(x as f64).unwrap())
                    .collect::<Vec<_>>(),
            );
            let b_f = LieSeries::<u8, NotNan<f64>>::new(
                words.clone(),
                b_int
                    .iter()
                    .map(|&x| NotNan::new(x as f64).unwrap())
                    .collect::<Vec<_>>(),
            );
            let ab_f = a_f.commutator(&b_f);

            // Generic exact path.
            let a_r = LieSeries::<u8, Ratio<i128>>::new(
                words.clone(),
                a_int.iter().map(|&x| Ratio::from_integer(x)).collect(),
            );
            let b_r = LieSeries::<u8, Ratio<i128>>::new(
                words.clone(),
                b_int.iter().map(|&x| Ratio::from_integer(x)).collect(),
            );
            let ab_r = a_r.commutator(&b_r);

            for (x, y) in ab_f.coefficients.iter().zip(&ab_r.coefficients) {
                assert_eq!(
                    x.into_inner(),
                    y.to_f64().unwrap(),
                    "d={d} m={m}: raw float vs exact rationals"
                );
            }
        }
    }

    /// The fused kernel evaluates both bracket orientations of every canonical
    /// pair, so `[A, A]` must vanish exactly — every pair contributes
    /// `c * (a_min * a_max - a_max * a_min) = 0`.
    #[test]
    fn commutator_of_series_with_itself_is_zero() {
        use lyndon_rs::lyndon::{LyndonBasis, Sort};

        for (d, m) in [(2usize, 6usize), (3, 5)] {
            let words: Vec<LyndonWord<u8>> =
                LyndonBasis::<u8>::new(d, Sort::Lexicographical).generate_basis(m);
            let coefficients: Vec<_> = (0..words.len())
                .map(|i: usize| NotNan::new(((i % 7 + 1) as f64) / 3.0 - 1.1).unwrap())
                .collect();
            let series: LieSeries<u8, NotNan<f64>> = LieSeries::new(words, coefficients);
            let result: LieSeries<u8, NotNan<f64>> = series.commutator(&series);
            assert!(
                result.coefficients.iter().all(|c| c.is_zero()),
                "non-zero coefficient in [A, A] for d={d}, m={m}"
            );
        }
    }

    /// ADVERSARIAL (write-access division): the kernel sweep must equal an
    /// INDEPENDENT per-word fan-in reference — built by walking the table's
    /// entries in table order with presence resolved straight from the
    /// support lists — bit for bit, across coefficient types (raw-float
    /// fast path and exact rationals), layouts (public direct and
    /// class-contiguous), thread counts, and adversarial inputs (all-zero,
    /// single-hot, planted exact cancellations). Distinct failure: this
    /// guards arithmetic and per-word accumulation order, not gating
    /// structure (the transposition test in `gating`).
    #[test]
    fn write_class_sweep_matches_independent_fanin_reference() {
        use lyndon_rs::lyndon::{LyndonBasis, Sort};
        use std::collections::BTreeMap;

        use super::super::ClassOrderedCommutation;

        // The reference: result[w] += c * term over the table's entries in
        // table order, term = (p_active ? a_i*b_j : 0) - (q_active ?
        // a_j*b_i : 0), entries with neither orientation skipped. Same add
        // sequence per word as the kernel's write classes, derived without
        // touching the gating.
        fn reference<U>(
            series: &LieSeries<u8, U>,
            a: &[U],
            a_present: &[bool],
            b: &[U],
            b_present: &[bool],
        ) -> Vec<U>
        where
            U: Clone + Default + Mul<Output = U> + AddAssign + Neg<Output = U> + PartialEq,
        {
            let table = &series.feasible_decompositions;
            let entries = table.entries();
            let coeffs = table.decomp_coeffs();
            let mut acc: BTreeMap<usize, U> = BTreeMap::new();
            for (ei, entry) in entries[..entries.len() - 1].iter().enumerate() {
                let (i, j) = (entry.i as usize, entry.j as usize);
                let p_active = a_present[i] && b_present[j];
                let q_active = a_present[j] && b_present[i];
                if !p_active && !q_active {
                    continue;
                }
                let term = if p_active {
                    let mut t = a[i].clone() * b[j].clone();
                    if q_active {
                        t += -(a[j].clone() * b[i].clone());
                    }
                    t
                } else {
                    -(a[j].clone() * b[i].clone())
                };
                let from = entry.decomp_start as usize;
                let to = entries[ei + 1].decomp_start as usize;
                let rs = table.degree_start(
                    table.degree_of(i) + table.degree_of(j),
                );
                for dp in from..to {
                    // Absolute decomp position -> the row's target: the
                    // relative array is stored per entry in row order.
                    let rel = table.decomp_indices_rel()[dp];
                    let w = rs + rel as usize;
                    let contrib = coeffs[dp].clone() * term.clone();
                    *acc.entry(w).or_insert_with(U::default) += contrib;
                }
            }
            let mut out = vec![U::default(); series.coefficients.len()];
            for (w, v) in acc {
                out[w] = v;
            }
            out
        }

        for (d, m) in [(2usize, 4usize), (3usize, 6usize)] {
            let words: Vec<LyndonWord<u8>> =
                LyndonBasis::<u8>::new(d, Sort::Lexicographical).generate_basis(m);
            let len = words.len();
            let table_deg = |i: usize| words[i].len();
            let _ = table_deg;

            // Adversarial value sets: zeros (support holes), signed values
            // that plant exact cancellations between entries feeding the
            // same word, and asymmetric magnitudes.
            let mk = |f: &dyn Fn(usize) -> f64| {
                (0..len)
                    .map(|i| NotNan::new(f(i)).unwrap())
                    .collect::<Vec<_>>()
            };
            let a_sets = [
                mk(&|_| 0.0),
                mk(&|i| if i == 0 { 1.0 } else { 0.0 }),
                mk(&|i| ((i % 5 + 1) as f64) / 3.0 - 1.0),
                mk(&|i| if i % 2 == 0 { 1.0 } else { -1.0 }),
            ];
            let b_sets = [
                mk(&|_| 0.0),
                mk(&|i| if i == 1 { 1.0 } else { 0.0 }),
                mk(&|i| ((i * 3 + 2) % 7) as f64 / 2.0 - 1.5),
                mk(&|i| if i % 3 == 0 { 2.0 } else { -0.5 }),
            ];

            let cutoff_len = len; // full lists; the kernel filters by value
            for a_coefficients in &a_sets {
                for b_coefficients in &b_sets {
                    let a: LieSeries<u8, NotNan<f64>> =
                        LieSeries::new(words.clone(), a_coefficients.clone());
                    let b: LieSeries<u8, NotNan<f64>> =
                        LieSeries::new(words.clone(), b_coefficients.clone());
                    let a_nonzero = a.nonzero_coefficient_indices(a_coefficients);
                    let b_nonzero = b.nonzero_coefficient_indices(b_coefficients);
                    let ap = {
                        let mut p = vec![false; cutoff_len];
                        for &i in &a_nonzero {
                            p[i] = true;
                        }
                        p
                    };
                    let bp = {
                        let mut p = vec![false; cutoff_len];
                        for &i in &b_nonzero {
                            p[i] = true;
                        }
                        p
                    };
                    let expected =
                        reference(&a, a_coefficients, &ap, b_coefficients, &bp);

                    for threads in [1usize, 4usize] {
                        let pool = rayon::ThreadPoolBuilder::new()
                            .num_threads(threads)
                            .build()
                            .expect("pool");
                        pool.install(|| {
                            // Public direct kernel.
                            let mut out =
                                vec![NotNan::<f64>::default(); len];
                            LieSeries::commutator_coefficients_with_nonzero(
                                &a,
                                a_coefficients,
                                &a_nonzero,
                                b_coefficients,
                                &b_nonzero,
                                &mut out,
                            );
                            assert_eq!(
                                out, expected,
                                "public kernel differs (d={d}, m={m}, t={threads}, \
                                 a={:?}, b={:?})",
                                a_coefficients.iter().map(|x| x.into_inner()).collect::<Vec<_>>(),
                                b_coefficients.iter().map(|x| x.into_inner()).collect::<Vec<_>>()
                            );

                            // Class-contiguous kernel: operands class-
                            // ordered, result mapped back through `inv`.
                            let order = a.class_order();
                            let a_cls = a.class_coefficients(&order, a_coefficients);
                            let b_cls = b.class_coefficients(&order, b_coefficients);
                            let a_nz_cls: Vec<usize> =
                                a_nonzero.iter().map(|&i| order.inv()[i] as usize).collect();
                            let b_nz_cls: Vec<usize> =
                                b_nonzero.iter().map(|&i| order.inv()[i] as usize).collect();
                            let mut out_cls =
                                vec![NotNan::<f64>::default(); len];
                            a.class_commutation(
                                &order,
                                &a_cls,
                                &a_nz_cls,
                                &b_cls,
                                &b_nz_cls,
                                &mut out_cls,
                                &mut GatingCache::default(),
                            );
                            let public = a.public_coefficients(&order, &out_cls);
                            assert_eq!(
                                public, expected,
                                "class kernel differs (d={d}, m={m}, t={threads})"
                            );
                        });
                    }
                }
            }
        }
    }

    /// ADVERSARIAL (write-access division): a word class stays OWNED under
    /// total cancellation — `a == b` makes every term exactly zero, yet
    /// every reachable word is still swept (one single-writer store of the
    /// exact zero) and still reported by the collecting variant: the
    /// reported targets are a WRITE-set superset, never a value-filtered
    /// set. Distinct failure: guards target ownership and the superset
    /// property, not arithmetic (the fan-in test above).
    #[test]
    fn write_class_canceled_word_stays_owned_and_reported() {
        use lyndon_rs::lyndon::{LyndonBasis, Sort};

        for (d, m) in [(2usize, 5usize), (3usize, 4usize)] {
            let words: Vec<LyndonWord<u8>> =
                LyndonBasis::<u8>::new(d, Sort::Lexicographical).generate_basis(m);
            let len = words.len();
            let coeffs: Vec<_> = (0..len)
                .map(|i| NotNan::new(((i % 9 + 1) as f64) / 4.0 - 1.1).unwrap())
                .collect();
            let a: LieSeries<u8, NotNan<f64>> =
                LieSeries::new(words.clone(), coeffs.clone());
            let b: LieSeries<u8, NotNan<f64>> = LieSeries::new(words, coeffs);
            let full: Vec<usize> = (0..len).collect();

            // [A, A] = 0 exactly (float products are commutative), but
            // every word with a feasible decomposition is still written.
            let mut out = vec![NotNan::<f64>::default(); len];
            LieSeries::commutator_coefficients_with_nonzero(
                &a,
                &a.coefficients,
                &full,
                &b.coefficients,
                &full,
                &mut out,
            );
            for (w, v) in out.iter().enumerate() {
                assert_eq!(*v, NotNan::<f64>::default(), "word {w} not exactly zero");
            }

            // The collecting variant must still report every written word:
            // value-filtered targets would be unsound for list reuse (a
            // canceled target can go nonzero again when values change).
            let mut dirty = vec![0u64; len / 64 + 1];
            let mut targets = Vec::new();
            LieSeries::commutator_coefficients_with_nonzero_collecting(
                &a,
                &a.coefficients,
                &full,
                &b.coefficients,
                &full,
                &mut out,
                &mut dirty,
                &mut targets,
            );
            // Independent expectation: every word targeted by any active
            // contribution of any entry (full supports -> every entry).
            let table = &a.feasible_decompositions;
            let entries = table.entries();
            let decomp = table.decomp_indices_rel();
            let degrees = table.index_degrees();
            let mut expected: Vec<usize> = Vec::new();
            for (ei, entry) in entries[..entries.len() - 1].iter().enumerate() {
                let (i, j) = (entry.i as usize, entry.j as usize);
                let from = entry.decomp_start as usize;
                let to = entries[ei + 1].decomp_start as usize;
                let rs = table.degree_start(degrees[i] as usize + degrees[j] as usize);
                for &rel in &decomp[from..to] {
                    expected.push(rs + rel as usize);
                }
            }
            expected.sort_unstable();
            expected.dedup();
            // The sink reports kernel-VISIBLE positions only: the degree-
            // `max_degree` tail lies above the support cutoff (nothing can
            // consume those values in a truncated BCH fold), so both sides
            // filter at `degree_start(max_degree)`.
            let cutoff = table.degree_start(table.max_degree());
            expected.retain(|&w| w < cutoff);
            targets.sort_unstable();
            assert_eq!(
                targets, expected,
                "collecting targets lost canceled words (write-set superset violated)"
            );
        }
    }
}
