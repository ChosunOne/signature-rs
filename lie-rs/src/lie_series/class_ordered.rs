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

#[cfg(test)]
mod tests {
    use super::*;
    use ordered_float::NotNan;

    /// The trait's internal working mode (`class_commutation`) must agree
    /// bit-for-bit with the direct kernel after one
    /// `public_coefficients` epilogue, at every thread count; the
    /// collecting variant must report the class-indexed image of the
    /// direct layout's first-touch sequence.
    #[test]
    fn class_commutation_round_trip_is_bit_identical() {
        use lyndon_rs::lyndon::{LyndonBasis, Sort};

        for (d, m, force) in [(2usize, 4usize, true), (3usize, 10usize, false)] {
            let words: Vec<LyndonWord<u8>> =
                LyndonBasis::<u8>::new(d, Sort::Lexicographical).generate_basis(m);
            let a_coefficients: Vec<_> = (0..words.len())
                .map(|i: usize| NotNan::new(((i % 11 + 1) as f64) / 7.0 - 0.9).unwrap())
                .collect();
            let b_coefficients: Vec<_> = (0..words.len())
                .map(|i: usize| NotNan::new(((i * 5 + 3) % 13) as f64 / 6.0 - 1.3).unwrap())
                .collect();
            let mut direct = LieSeries::new(words.clone(), a_coefficients.clone());
            let class = LieSeries::new(words, b_coefficients.clone());
            if force {
                Arc::make_mut(&mut direct.feasible_decompositions).clear_class_order();
            }
            let order = class.class_order();
            assert_eq!(order.inv().len(), class.coefficients.len());

            let a_cls = class.class_coefficients(&order, &a_coefficients);
            let b_cls = class.class_coefficients(&order, &b_coefficients);
            let inv_of = |i: usize| order.inv()[i] as usize;
            let a_nonzero = direct.nonzero_coefficient_indices(&a_coefficients);
            let b_nonzero = direct.nonzero_coefficient_indices(&b_coefficients);
            let a_nz_cls: Vec<usize> = a_nonzero.iter().copied().map(inv_of).collect();
            let b_nz_cls: Vec<usize> = b_nonzero.iter().copied().map(inv_of).collect();

            let run = |threads: usize, series: &LieSeries<u8, NotNan<f64>>| {
                let pool = rayon::ThreadPoolBuilder::new()
                    .num_threads(threads)
                    .build()
                    .expect("thread pool");
                let mut out = vec![NotNan::<f64>::default(); series.coefficients.len()];
                pool.install(|| {
                    LieSeries::commutator_coefficients(
                        series,
                        &a_coefficients,
                        &b_coefficients,
                        &mut out,
                    )
                });
                out
            };

            let reference = run(1, &direct);
            for threads in [1usize, 2, 4, 8] {
                let pool = rayon::ThreadPoolBuilder::new()
                    .num_threads(threads)
                    .build()
                    .expect("thread pool");
                let mut result = vec![NotNan::<f64>::default(); class.coefficients.len()];
                pool.install(|| {
                    class.class_commutation(
                        &order,
                        &a_cls,
                        &a_nz_cls,
                        &b_cls,
                        &b_nz_cls,
                        &mut result,
                        &mut GatingCache::default(),
                    )
                });
                let public = class.public_coefficients(&order, &result);
                assert_eq!(
                    reference, public,
                    "class round trip differs for d={d}, m={m}, threads={threads}"
                );
            }

            // Collecting variant: class-indexed first-touch sequence.
            let mut direct_result = vec![NotNan::<f64>::default(); direct.coefficients.len()];
            let mut direct_dirty = vec![0u64; direct.coefficients.len() / 64 + 1];
            let mut direct_targets = Vec::new();
            LieSeries::commutator_coefficients_with_nonzero_collecting(
                &direct,
                &a_coefficients,
                &a_nonzero,
                &b_coefficients,
                &b_nonzero,
                &mut direct_result,
                &mut direct_dirty,
                &mut direct_targets,
            );
            let mut class_result = vec![NotNan::<f64>::default(); class.coefficients.len()];
            let mut class_dirty = vec![0u64; class.coefficients.len() / 64 + 1];
            let mut class_targets = Vec::new();
            class.class_commutation_with_nonzero_collecting(
                &order,
                &a_cls,
                &a_nz_cls,
                &b_cls,
                &b_nz_cls,
                &mut class_result,
                &mut class_dirty,
                &mut class_targets,
            );
            for (k, &src) in order.inv().iter().enumerate() {
                assert_eq!(
                    direct_result[k], class_result[src as usize],
                    "collecting result differs for d={d}, m={m} at public {k}"
                );
            }
            // Class targets are internal positions: map them back with
            // `perm` (internal -> public) before comparing. Under the
            // write-access division each variant emits its targets sorted
            // ascending in ITS OWN layout's position order (one store per
            // word class, classes in position order), so the public
            // sequences differ by the layout permutation — the invariant is
            // the same target SET, each position reported exactly once
            // (per-word accumulation order is guarded by the bit-identical
            // result comparison above).
            let perm = order.perm();
            let mut relabeled: Vec<usize> = class_targets
                .iter()
                .copied()
                .map(|p| perm[p] as usize)
                .collect();
            relabeled.sort_unstable();
            let mut direct_sorted = direct_targets.clone();
            direct_sorted.sort_unstable();
            assert_eq!(direct_sorted, relabeled, "collecting targets differ");
        }
    }

    /// The class-contiguous scatter layout must produce bit-identical
    /// results to the direct layout, at every thread count, on both kernel
    /// entry points: the layout is a pure relabeling of write addresses and
    /// never reorders the per-word accumulation sequence. `(2, 4)` forces
    /// the layout on a small table (fast build); `(3, 10)` qualifies
    /// through the real slice-size threshold (its degree-10 slice exceeds
    /// 4096 words).
    #[test]
    fn commutator_is_bit_identical_across_scatter_layouts() {
        use lyndon_rs::lyndon::{LyndonBasis, Sort};

        for (d, m, force) in [(2usize, 4usize, true), (3usize, 10usize, false)] {
            let words: Vec<LyndonWord<u8>> =
                LyndonBasis::<u8>::new(d, Sort::Lexicographical).generate_basis(m);
            let a_coefficients: Vec<_> = (0..words.len())
                .map(|i: usize| NotNan::new(((i % 11 + 1) as f64) / 7.0 - 0.9).unwrap())
                .collect();
            let b_coefficients: Vec<_> = (0..words.len())
                .map(|i: usize| NotNan::new(((i * 5 + 3) % 13) as f64 / 6.0 - 1.3).unwrap())
                .collect();
            let mut direct = LieSeries::new(words.clone(), a_coefficients.clone());
            let mut class = LieSeries::new(words, b_coefficients.clone());
            if force {
                Arc::make_mut(&mut class.feasible_decompositions).force_class_order();
            } else {
                // (3, 10) auto-qualifies through the threshold: strip the
                // reference's layout so the pair compares direct vs class
                // on the same table.
                Arc::make_mut(&mut direct.feasible_decompositions).clear_class_order();
            }
            assert!(
                direct
                    .feasible_decompositions
                    .cached_class_order()
                    .is_none(),
                "expected the reference series to run the direct layout"
            );
            assert!(
                class.feasible_decompositions.cached_class_order().is_some(),
                "expected the layout series to carry a class order"
            );
            let len = direct.coefficients.len();

            let run = |series: &LieSeries<u8, NotNan<f64>>, threads: usize| {
                let pool = rayon::ThreadPoolBuilder::new()
                    .num_threads(threads)
                    .build()
                    .expect("thread pool");
                let mut out = vec![NotNan::<f64>::default(); len];
                pool.install(|| {
                    LieSeries::commutator_coefficients(
                        series,
                        &a_coefficients,
                        &b_coefficients,
                        &mut out,
                    )
                });
                out
            };

            let reference = run(&direct, 1);
            for threads in [1usize, 2, 4, 8] {
                let out = run(&class, threads);
                assert_eq!(
                    reference, out,
                    "class layout differs from direct for d={d}, m={m}, threads={threads}"
                );
            }

            // Collecting entry point: result, nonzero set, and first-touch
            // order must all match the direct layout exactly.
            let collect = |series: &LieSeries<u8, NotNan<f64>>| {
                let a_nonzero = series.nonzero_coefficient_indices(&a_coefficients);
                let b_nonzero = series.nonzero_coefficient_indices(&b_coefficients);
                let mut result = vec![NotNan::<f64>::default(); len];
                let mut dirty = vec![0u64; len / 64 + 1];
                let mut targets = Vec::new();
                LieSeries::commutator_coefficients_with_nonzero_collecting(
                    series,
                    &a_coefficients,
                    &a_nonzero,
                    &b_coefficients,
                    &b_nonzero,
                    &mut result,
                    &mut dirty,
                    &mut targets,
                );
                (result, targets)
            };
            let (direct_result, direct_targets) = collect(&direct);
            let (class_result, class_targets) = collect(&class);
            assert_eq!(direct_result, class_result, "collecting result differs");
            assert_eq!(direct_targets, class_targets, "collecting targets differ");
        }
    }
}
