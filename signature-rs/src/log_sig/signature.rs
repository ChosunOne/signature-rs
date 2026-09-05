use super::*;
use super::engine::{
    SeriesTemplate, TOURNAMENT_LEAF_CHUNK, TOURNAMENT_MAX_LEAVES, TOURNAMENT_MIN_DISPLACEMENTS,
    audit_coefficients_no_nan, cohort_merge_group, fold_batch_sequential, fold_one_displacement,
    leaf_group_engine, return_pooled_dag, scalar_merge_fold, take_pooled_dag,
};

/// Represents a computed log signature of a path.
///
/// A log signature captures the essential geometric and algebraic features of a path
/// through its representation as a series in the free Lie algebra. This structure
/// contains both the computed coefficients and the Baker-Campbell-Hausdorff series
/// needed for concatenation operations.
pub struct LogSignature<T, U> {
    /// The main series containing the log signature coefficients.
    pub series: LieSeries<T, U>,
    /// The BCH series used for concatenating log signatures.
    pub bch_series: LieSeries<u8, U>,
    /// Compiled evaluation plan for [`LogSignature::concatenate_assign`]
    pub(super) dag: CommutatorDag<U>,
}

impl<T: Debug, U: Debug> Debug for LogSignature<T, U> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("LogSignature")
            .field("series", &self.series)
            .field("bch_series", &self.bch_series)
            .field("dag", &self.dag)
            .finish()
    }
}

impl<T: Clone, U: Clone> Clone for LogSignature<T, U> {
    fn clone(&self) -> Self {
        Self {
            series: self.series.clone(),
            bch_series: self.bch_series.clone(),
            dag: self.dag.clone_shallow(),
        }
    }
}

impl<T, U> Index<usize> for LogSignature<T, U> {
    type Output = U;

    fn index(&self, index: usize) -> &Self::Output {
        &self.series[index]
    }
}

impl<
    T: Clone + Ord + Generator + Hash,
    U: Clone + One + Zero + Eq + MulAssign + Neg<Output = U> + Hash,
> Index<LyndonWord<T>> for LogSignature<T, U>
{
    type Output = U;

    fn index(&self, index: LyndonWord<T>) -> &Self::Output {
        &self.series[index]
    }
}

impl<
    T: Clone + Ord + Generator + Hash,
    U: Clone + One + Zero + Eq + MulAssign + Neg<Output = U> + Hash,
> Index<&LyndonWord<T>> for LogSignature<T, U>
{
    type Output = U;

    fn index(&self, index: &LyndonWord<T>) -> &Self::Output {
        &self.series[index]
    }
}

impl<T, U> IndexMut<usize> for LogSignature<T, U> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.series[index]
    }
}

impl<
    T: Clone + Ord + Generator + Hash,
    U: Clone + One + Zero + Eq + MulAssign + Neg<Output = U> + Hash,
> IndexMut<LyndonWord<T>> for LogSignature<T, U>
{
    fn index_mut(&mut self, index: LyndonWord<T>) -> &mut Self::Output {
        &mut self.series[index]
    }
}

impl<
    T: Clone + Ord + Generator + Hash,
    U: Clone + One + Zero + Eq + MulAssign + Neg<Output = U> + Hash,
> IndexMut<&LyndonWord<T>> for LogSignature<T, U>
{
    fn index_mut(&mut self, index: &LyndonWord<T>) -> &mut Self::Output {
        &mut self.series[index]
    }
}

impl<
    T: Clone + Ord + Generator + Hash,
    U: Clone
        + Mul<Output = U>
        + MulAssign
        + AddAssign
        + Hash
        + Eq
        + Default
        + One
        + Zero
        + Neg<Output = U>
        + 'static,
> LogSignature<T, U>
{
    /// Concatenates two log signatures using the Baker-Campbell-Hausdorff formula.
    ///
    /// This operation computes the log signature of the path formed by concatenating
    /// the paths represented by `self` and `rhs`. The result captures the geometry
    /// of the combined path.
    #[must_use]
    pub fn concatenate(&self, rhs: &Self) -> Self
    where
        U: Send + Sync,
    {
        let mut concatenated_log_sig = self.clone();

        concatenated_log_sig.concatenate_assign(rhs);

        concatenated_log_sig
    }

    /// In-place concatenation of log signatures.
    ///
    /// This is the mutable version of `concatenate` that modifies `self` instead
    /// of creating a new log signature. More memory-efficient for chaining operations.
    pub fn concatenate_assign(&mut self, rhs: &Self)
    where
        U: Send + Sync,
    {
        self.concatenate_assign_coefficients(&rhs.series.coefficients);
    }

    pub(crate) fn concatenate_assign_coefficients(&mut self, rhs_coefficients: &[U])
    where
        U: Send + Sync,
    {
        let original_coefficients = self.series.coefficients.clone();

        let a_nonzero = self
            .series
            .nonzero_coefficient_indices(&original_coefficients);
        let b_nonzero = self.series.nonzero_coefficient_indices(rhs_coefficients);

        self.dag.evaluate(
            &self.series,
            &original_coefficients,
            &a_nonzero,
            rhs_coefficients,
            &b_nonzero,
        );

        // The DAG accumulated every term in class-contiguous order; this
        // single call applies the BCH weights and the one epilogue back to
        // public basis order.
        self.dag
            .accumulate_terms(&mut self.series.coefficients, rhs_coefficients);
        self.audit_no_nan();
    }

    /// The commutation kernel's raw-float fast path does not check NaN per
    /// operation (see `lie_rs::raw_mul`'s NaN policy); the only way NaN can
    /// arise is coefficient overflow producing infinities whose combination
    /// cancels. Audit the accumulator once per fold step so that failure
    /// stays loud instead of silently persisting a NaN through the `NotNan`
    /// invariant.
    fn audit_no_nan(&self) {
        audit_coefficients_no_nan(&self.series.coefficients);
    }

    /// Folds every displacement in `rhss` into `self`, in order.
    ///
    /// Three strategies, selected per call:
    ///
    /// - **Cohort tournament** (cohort-capable coefficient types — raw
    ///   repr-transparent floats on AVX2 CPUs with the kill switch off,
    ///   see `cohort_capable` — and at least 2 displacements): the
    ///   balanced tournament below, with EVERY round's folds routed
    ///   through the SIMD-across-folds cohort engine so a fold costs a
    ///   fraction of a scalar walk — 2-4 folds share one plan walk and
    ///   one kernel pass at every level of the tree, from the leaf
    ///   round's very first fold (shared-support groups run their whole
    ///   chunks as one dispatch) through every merge round that has a
    ///   pair to fold. Only truly lone folds (the last merge of an odd
    ///   leaf count) stay scalar. The tree shape depends ONLY on
    ///   `rhss.len()` — never on the thread count — so a given input
    ///   always reduces along the same tree and results are
    ///   reproducible at any pool size.
    /// - **Scalar tournament** (at least `TOURNAMENT_MIN_DISPLACEMENTS`
    ///   displacements and `rayon::current_num_threads() > 1`, for
    ///   non-capable types such as `Ratio<i128>`): the same balanced
    ///   tree with every round on the scalar per-lane engine.
    /// - **Sequential** (everything else): the folds run one-by-one
    ///   until the accumulator's support is the full basis and the node
    ///   lists are steady for it; the remaining folds — which must share
    ///   one displacement support — then run as ONE continuous batch
    ///   dispatch (see `CommutatorDag::fold_batch`), which keeps the
    ///   accumulator class-ordered and reuses one plan and one job table
    ///   for the whole batch, turning the per-fold gather/accumulate
    ///   phases into parallel in-graph stages. This path is bit-identical
    ///   to folding each displacement with
    ///   `Self::concatenate_assign_coefficients`.
    ///
    /// Rounding caveat: f32/f64 accumulation is not associative, so any
    /// tournament's reassociated tree is NOT bit-identical to the
    /// sequential path — observed worst-case shifts are a few ulps
    /// (~1e-13 absolute on coefficients of magnitude ~1e2 in adversarial
    /// test batches). The caveat now covers ALL cohort-capable batch
    /// sizes (the cohort tournament engages from 2 displacements up);
    /// within one strategy results are bit-stable: the sequential path is
    /// bit-identical to folding each displacement with
    /// `Self::concatenate_assign_coefficients`, and every tournament's
    /// tree shape depends only on `rhss.len()`, so a given input reduces
    /// along the same tree at any pool size. Exact coefficient types
    /// (e.g. `Ratio<i128>`) are insensitive to association order and
    /// produce identical results on every path — and they keep the
    /// historical engagement gate (the cohort engine is float-only), so
    /// small exact-type batches still run the bit-identical sequential
    /// path.
    ///
    /// No nested-recursion guard is needed: the tournament's leaves call
    /// the free `fold_batch_sequential` engine directly — the engine
    /// contains no driver-selection logic, so it cannot re-enter this
    /// method — and neither it nor `concatenate_assign_coefficients`
    /// re-enters this driver.
    pub fn concatenate_batch_coefficients(&mut self, rhss: &[&[U]])
    where
        U: Send + Sync,
    {
        // The tournament-vs-sequential choice is a TREE decision and must
        // not depend on the cohort kill switch: the switch swaps the
        // engines INSIDE the tournament (cohort vs scalar leaves/merges),
        // never the association tree, so the kill-switch oracle runs
        // compare the same tree's two engines bit for bit. Hence the
        // type/CPU capability ([`cohort_type_capable`]) here, with the
        // switch read inside [`Self::tournament_reduce`].
        if rhss.len() >= 2 && cohort_type_capable::<U>() {
            let reduced = self.tournament_reduce(rhss);
            self.concatenate_assign_coefficients(&reduced);
            return;
        }
        if rhss.len() >= TOURNAMENT_MIN_DISPLACEMENTS && rayon::current_num_threads() > 1 {
            let reduced = self.tournament_reduce(rhss);
            self.concatenate_assign_coefficients(&reduced);
            return;
        }
        self.concatenate_batch_sequential(rhss);
    }

    /// The strictly sequential fold chain — the historical body of
    /// [`Self::concatenate_batch_coefficients`], kept verbatim as the
    /// low-displacement path and as the bit-identical reference the
    /// tournament's rounding caveat is stated against. The walk itself
    /// lives on the free function [`fold_batch_sequential`] (a bare
    /// (`series`, `dag`) pair) so the tournament's leaves can run the
    /// identical engine — warm-up folds plus one steady batch dispatch —
    /// on their private accumulators.
    pub(super) fn concatenate_batch_sequential(&mut self, rhss: &[&[U]])
    where
        U: Send + Sync,
    {
        fold_batch_sequential(&mut self.series, &mut self.dag, rhss);
    }

    /// Balanced-binary-tournament reduction of `rhss` to the coefficient
    /// array of the concatenated sub-path.
    ///
    /// Leaves are contiguous chunks of at least [`TOURNAMENT_LEAF_CHUNK`]
    /// displacements — adaptively
    /// `max(TOURNAMENT_LEAF_CHUNK, rhss.len().div_ceil(TOURNAMENT_MAX_LEAVES))`,
    /// a pure function of `rhss.len()` — each folded from a zeroed
    /// accumulator with a private DAG through the SEQUENTIAL batch driver
    /// ([`fold_batch_sequential`], the exact engine of
    /// [`Self::concatenate_batch_sequential`]). A leaf therefore pays the
    /// same warm-up + steady-batch treatment as the public sequential path:
    /// its support-growth folds and node-list rebuilds run once at the
    /// chunk's start and then amortize over the whole chunk, instead of
    /// being paid per displacement (the failure mode of the historical
    /// fixed 8-displacement leaf, whose chunks never outgrew the batch
    /// warm-up). Within a leaf the association order is bit-identical to
    /// the sequential engine's on the same chunk, so the tournament's
    /// reassociation dust enters ONLY at merge boundaries, and the leaf
    /// shape matches the historical fixed-8 tree for batches of up to
    /// `TOURNAMENT_MAX_LEAVES * TOURNAMENT_LEAF_CHUNK` displacements.
    ///
    /// Each merge round folds the right chunk result into the left one —
    /// adjacent halves only — so the reduction is exactly the balanced
    /// binary tree over the input order: chunk boundaries are positional
    /// (`leaf * chunk`), a short tail chunk and an odd pass-through are
    /// determined by `rhss.len()` alone, and no part of the shape is chosen
    /// by the scheduler (the chunk formula involves no thread count, so
    /// results are bit-reproducible across pool sizes). A raw
    /// single-displacement chunk needs no special case: it is folded from
    /// a zeroed full-length accumulator like any other leaf, which keeps
    /// every leaf result in the same full-length representation regardless
    /// of the displacement slices' length.
    ///
    /// Recursion: leaves call [`fold_batch_sequential`] — the engine
    /// itself, which contains no driver-selection logic — so the
    /// tournament cannot re-engage inside a leaf by construction; the
    /// selection exists only in
    /// [`Self::concatenate_batch_coefficients`].
    ///
    /// DAG pooling: one `Mutex<Vec<CommutatorDag<U>>>` serves leaves and
    /// merges alike. A pooled dag is a `clone_shallow` copy — it shares the
    /// immutable compiled plan (`Arc<DagStructure>`) and carries private
    /// scratch. `evaluate` re-derives the node lists through a collecting
    /// pass whenever a fold's atom supports differ from the dag's stored
    /// ones, so a pooled dag is immediately usable for arbitrary supports:
    /// reuse can only save that rebuild, never skip it. The per-dag
    /// `GatingCache` is keyed by support fingerprints and valid for this
    /// series' (`Arc`-shared) decomposition table, which every tournament
    /// accumulator shares with `self.series` — a stale entry can therefore
    /// never be hit under different supports. Pooling bounds the deep
    /// copies of node lists/dirty bitsets at one per concurrently running
    /// fold instead of one per fold.
    pub(super) fn tournament_reduce(&self, rhss: &[&[U]]) -> Vec<U>
    where
        U: Send + Sync,
    {
        use rayon::prelude::*;

        // Template accumulator: shares the basis and tables with
        // `self.series` (Arc bumps only — see `SeriesTemplate`'s safety
        // contract); every task installs its own zeroed or chunk-result
        // coefficients before folding. The closure captures the wrapper
        // through this shared reference (precise capture would otherwise
        // reach past the wrapper to `&LieSeries`, whose `Send`/`Sync`
        // depend on `T`).
        let template = SeriesTemplate(self.series.clone());
        let template = &template;
        // Captured by the tasks instead of `self` itself: a closure using
        // `self.dag` captures the whole `&LogSignature` (reference derefs
        // defeat precise capture), which drags the `T`-bearing series into
        // the auto-trait requirements.
        let source_dag = &self.dag;
        let basis_len = self.series.coefficients.len();
        // One private-dag pool for leaves and merges alike; it grows to at
        // most the peak number of concurrently running folds and is dropped
        // with the tournament.
        let dag_pool: Mutex<Vec<CommutatorDag<U>>> = Mutex::new(Vec::new());

        // The cohort capability: raw repr-transparent float coefficients +
        // AVX2 + kill switch (see `cohort_capable`). False keeps every
        // round on the scalar paths below — the pre-cohort tournament,
        // bit for bit.
        let cohort_on = cohort_capable::<U>();

        // Adaptive leaf chunk: large enough that each leaf's sequential
        // warm-up amortizes over its chunk, capped at ~32 leaves to keep
        // in-flight oversubscription modest — see `TOURNAMENT_MAX_LEAVES`.
        // Cohort-capable types shrink the chunk below the tournament floor
        // for small batches so even 2-3-displacement batches form 2-4-lane
        // leaf groups — the SIMD-across-folds engine then covers EVERY fold
        // of the batch (the whole point of the cohort routing), not just
        // long batches' merge rounds. Both formulas are pure functions of
        // `rhss.len()` — never of the thread count — so the tree shape (and
        // therefore results, bit for bit) stays reproducible at any pool
        // size within each regime.
        let chunk = if rhss.len() < TOURNAMENT_MIN_DISPLACEMENTS {
            rhss.len().div_ceil(4).max(1)
        } else {
            TOURNAMENT_LEAF_CHUNK.max(rhss.len().div_ceil(TOURNAMENT_MAX_LEAVES))
        };

        // Cohort engagement is unconditional (capability-gated only, see
        // `cohort_capable`): the post-F2 single-phase per-unit sweep made the
        // 4-lane SoA engine the LIGHT path — one shared plan walk per group
        // and SIMD-4 math — so it wins or ties at every measured pool width,
        // including the regime the former wide-pool gate protected. The old
        // gate (`threads <= 16 || terms >= 256 || table >= 1000 || groups >=
        // threads`, consts COHORT_WIDE_POOL_*) dated to the pre-F2 two-phase
        // engine whose stage chains serialized at width; re-measured after
        // F2 (prof_concat e2e driver, 3 interleaved reps) it only misgated:
        // forcing the cohort past the old gate at 32 workers measured
        // 8x3 −29% (1.50→1.09ms), 2x8 −18% (4.87→3.99ms), 3x8 −30%
        // (scalar counterfactual 63→48ms), 2x12 −46% (322→216ms), 4x6 −8%;
        // no cell favored the scalar fallback. The lone-pair merge task
        // below stays scalar by construction (a 1-lane cohort is the scalar
        // engine).

        // A leaf folded through the sequential batch engine (the scalar
        // leaf task body — also the partial group's path).
        // Leaf round: chunks fold their disjoint displacement ranges
        // independently (the tree's bottom level), through the steady-chunk
        // leaf engine — one reachable-support plan per leaf, then the
        // leaf's whole chunk as ONE batch dispatch. Groups of two-to-four
        // leaves with support-uniform displacements run as ONE SoA cohort
        // dispatch (up to 4 folds per plan walk from the very first fold);
        // every other group — and every group when the cohort engine is
        // unavailable — folds its lanes through the scalar batch engine
        // under the same plans, bit-identically.
        let leaf_count = rhss.len().div_ceil(chunk);
        let cohort_leaf = cohort_on;
        let mut level: Vec<Vec<U>> = (0..leaf_count.div_ceil(4))
            .into_par_iter()
            .flat_map_iter(|g| {
                let first = g * 4;
                let last = (first + 4).min(leaf_count);
                leaf_group_engine(
                    template,
                    source_dag,
                    &dag_pool,
                    basis_len,
                    rhss,
                    chunk,
                    first,
                    last - first,
                    cohort_leaf,
                )
            })
            .collect();

        // Merge rounds: pairwise folds over adjacent results. A single
        // left-over result passes through to the next round; the clone
        // (at most one basis-length copy per round) keeps the parallel map
        // borrow-free and changes nothing about the tree shape.
        while level.len() > 1 {
            let passthrough = if level.len() % 2 == 1 {
                level.pop()
            } else {
                None
            };
            let n_pairs = level.len() / 2;
            level = if cohort_on {
                (0..n_pairs.div_ceil(4))
                    .into_par_iter()
                    .flat_map_iter(|g| {
                        let first = g * 4;
                        let last = (first + 4).min(n_pairs);
                        if last - first >= 2 {
                            cohort_merge_group(
                                template,
                                source_dag,
                                &dag_pool,
                                &level,
                                &(first..last).collect::<Vec<usize>>(),
                            )
                        } else {
                            // A lone pair: today's scalar merge task.
                            let p = first;
                            vec![scalar_merge_fold(
                                template,
                                source_dag,
                                &dag_pool,
                                &level[2 * p],
                                &level[2 * p + 1],
                            )]
                        }
                    })
                    .collect()
            } else {
                level
                    .par_chunks(2)
                    .map(|pair| match pair {
                        [left, right] => {
                            let mut series = template.0.clone();
                            series.coefficients = left.clone();
                            let mut dag = take_pooled_dag(&dag_pool, source_dag);
                            fold_one_displacement(&mut series, &mut dag, right);
                            return_pooled_dag(&dag_pool, dag);
                            series.coefficients
                        }
                        [single] => single.clone(),
                        _ => unreachable!("par_chunks(2) yields chunks of length 1 or 2"),
                    })
                    .collect()
            };
            if let Some(single) = passthrough {
                level.push(single);
            }
        }

        level.pop().expect("at least one leaf chunk")
    }

    /// Folds every log signature in `rhss` into `self`, in order — the
    /// batch form of [`Self::concatenate_assign`]. Same batching behavior
    /// and guarantees as [`Self::concatenate_batch_coefficients`], including
    /// the adaptive tournament reduction and its rounding caveat for
    /// non-exact coefficient types.
    pub fn concatenate_batch(&mut self, rhss: &[Self])
    where
        U: Send + Sync,
    {
        let coeffs: Vec<&[U]> = rhss
            .iter()
            .map(|r| r.series.coefficients.as_slice())
            .collect();
        self.concatenate_batch_coefficients(&coeffs);
    }
}
