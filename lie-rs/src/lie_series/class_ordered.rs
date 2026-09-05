use super::*;

/// Class-contiguous working mode for the commutation kernel.
///
/// The class ordering depends only on the basis: every series over the
/// same basis rearranges coefficients identically, so one [`ClassOrder`]
/// handle — cheap to clone, `Arc` internally — amortizes the O(basis)
/// planning across operand series, across every kernel call of a fold,
/// and across batches of folds. Obtain the handle once via
/// [`Self::class_order`], keep all intermediate work in class order, and
/// pay for exactly one [`Self::public_coefficients`] epilogue at the end.
pub trait ClassOrderedCommutation<T, U> {
    /// The basis' class-contiguous ordering: created on first request
    /// (then cached on the series' table and shared), or prebuilt at
    /// series construction when a degree slice outgrows L1.
    fn class_order(&self) -> Arc<ClassOrder>;

    /// Rearranges public-order coefficients into class-contiguous order
    /// (the fold's operand conversion, once per fold).
    fn class_coefficients(&self, order: &ClassOrder, coefficients: &[U]) -> Vec<U>;

    /// The final epilogue: class-contiguous coefficients back to public
    /// basis order. Call once, after the last class-ordered kernel call.
    fn public_coefficients(&self, order: &ClassOrder, class_coefficients: &[U]) -> Vec<U>;

    /// Commutator entirely in class-contiguous order: `a`/`b` are
    /// class-ordered coefficient slices, the support lists class-indexed,
    /// and `result` receives the class-ordered sum directly — no scratch,
    /// no permutation. `order` must describe this series' basis.
    fn class_commutation(
        &self,
        order: &ClassOrder,
        a: &[U],
        a_nonzero: &[usize],
        b: &[U],
        b_nonzero: &[usize],
        result: &mut [U],
        cache: &mut GatingCache,
    );

    /// Batch form of [`Self::class_commutation`] for folds: jobs carry
    /// class-ordered operands, class-indexed support lists, and
    /// class-ordered result buffers.
    fn class_commutation_batch(
        &self,
        order: &ClassOrder,
        jobs: &mut [KernelJob<'_, U>],
        cache: &mut GatingCache,
    );

    /// Whole-fold form: `levels[l]` holds the jobs of dependency level `l`
    /// (operand slices may alias results of strictly earlier levels). The
    /// engine plans class-contiguous balanced packs for every level up front
    /// and runs the fold as one parallel pack-slot walk with level counters
    /// — no per-level dispatches.
    fn class_commutation_batch_fold(
        &self,
        order: &ClassOrder,
        levels: &mut [Vec<KernelJob<'_, U>>],
        cache: &mut GatingCache,
    ) {
        for level in levels.iter_mut() {
            if level.len() == 1 {
                self.class_commutation(
                    order,
                    // SAFETY: single-job serial path (exclusive).
                    unsafe { std::slice::from_raw_parts(level[0].a, level[0].a_len) },
                    level[0].a_nonzero,
                    unsafe { std::slice::from_raw_parts(level[0].b, level[0].b_len) },
                    level[0].b_nonzero,
                    unsafe { std::slice::from_raw_parts_mut(level[0].result, level[0].result_len) },
                    cache,
                );
            } else {
                self.class_commutation_batch(order, level, cache);
            }
        }
    }

    /// Collecting variant used by folds to (re)derive each result's
    /// scatter-target superset: the reported indices are class-indexed.
    #[allow(clippy::too_many_arguments)]
    fn class_commutation_with_nonzero_collecting(
        &self,
        order: &ClassOrder,
        a: &[U],
        a_nonzero: &[usize],
        b: &[U],
        b_nonzero: &[usize],
        result: &mut [U],
        dirty: &mut [u64],
        targets: &mut Vec<usize>,
    );
}

impl<T, U> ClassOrderedCommutation<T, U> for LieSeries<T, U>
where
    T: Clone + Ord + Generator + Hash,
    U: Clone
        + Default
        + One
        + Zero
        + Eq
        + MulAssign
        + Neg<Output = U>
        + Hash
        + Mul<Output = U>
        + AddAssign
        + Send
        + Sync
        + 'static,
{
    fn class_order(&self) -> Arc<ClassOrder> {
        self.feasible_decompositions.ensure_class_order()
    }

    fn class_coefficients(&self, order: &ClassOrder, coefficients: &[U]) -> Vec<U> {
        let inv = order.inv();
        let mut out = vec![U::default(); coefficients.len()];
        for (w, c) in coefficients.iter().enumerate() {
            out[inv[w] as usize] = c.clone();
        }
        out
    }

    fn public_coefficients(&self, order: &ClassOrder, class_coefficients: &[U]) -> Vec<U> {
        let inv = order.inv();
        let mut out = vec![U::default(); class_coefficients.len()];
        for (k, &src) in inv.iter().enumerate() {
            out[k] = class_coefficients[src as usize].clone();
        }
        out
    }

    fn class_commutation(
        &self,
        order: &ClassOrder,
        a: &[U],
        a_nonzero: &[usize],
        b: &[U],
        b_nonzero: &[usize],
        result: &mut [U],
        cache: &mut GatingCache,
    ) {
        let mut job = KernelJob {
            a: a.as_ptr(),
            a_len: a.len(),
            a_nonzero,
            b: b.as_ptr(),
            b_len: b.len(),
            b_nonzero,
            result: result.as_mut_ptr(),
            result_len: result.len(),
            a_shift: IDENTITY_SHIFTS.as_ptr(),
            b_shift: IDENTITY_SHIFTS.as_ptr(),
            r_shift: IDENTITY_SHIFTS.as_ptr(),
        };
        commutator_coefficients_class_batch_with_cache(
            self,
            order,
            std::slice::from_mut(&mut job),
            cache,
        );
    }

    fn class_commutation_batch(
        &self,
        order: &ClassOrder,
        jobs: &mut [KernelJob<'_, U>],
        cache: &mut GatingCache,
    ) {
        commutator_coefficients_class_batch_with_cache(self, order, jobs, cache);
    }

    fn class_commutation_with_nonzero_collecting(
        &self,
        order: &ClassOrder,
        a: &[U],
        a_nonzero: &[usize],
        b: &[U],
        b_nonzero: &[usize],
        result: &mut [U],
        dirty: &mut [u64],
        targets: &mut Vec<usize>,
    ) {
        LieSeries::commutator_coefficients_class_with_nonzero_collecting(
            self, order, a, a_nonzero, b, b_nonzero, result, dirty, targets,
        );
    }

    fn class_commutation_batch_fold(
        &self,
        order: &ClassOrder,
        levels: &mut [Vec<KernelJob<'_, U>>],
        cache: &mut GatingCache,
    ) {
        // The one-dispatch pack-slot walk is the default: its packs are
        // claimed dynamically, so a waiter never blocks on runnable fold
        // work and the walk stays live even when the pool is shared with
        // outer parallelism. Two cases route elsewhere: with a single
        // worker nothing can run in parallel, so the walk's planning (task
        // lists, packs, claim counters) would be pure overhead — the
        // per-level dispatch's serial branch sweeps the identical order
        // with no planning; and `FOLD_ENGINE=0` keeps the per-level
        // dispatch for A/B comparison.
        static PER_LEVEL_DISPATCH: OnceLock<bool> = OnceLock::new();
        let per_level = *PER_LEVEL_DISPATCH
            .get_or_init(|| std::env::var_os("FOLD_ENGINE").is_some_and(|v| v == *"0"));
        if per_level || rayon::current_num_threads() == 1 {
            for level in levels.iter_mut() {
                self.class_commutation_batch(order, level, cache);
            }
        } else {
            commutator_coefficients_class_fold_with_cache(self, order, levels, cache);
        }
    }
}
