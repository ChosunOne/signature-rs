use super::*;

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
