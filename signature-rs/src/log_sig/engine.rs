use super::*;

/// True exactly for a NaN payload. `x != x` holds for no value except NaN,
/// which keeps the check available for every coefficient type without adding
/// a float-trait bound.
#[inline]
#[allow(clippy::eq_op)]
fn is_nan_value<U: PartialEq>(c: &U) -> bool {
    c != c
}

/// Minimum number of displacement slices for
/// [`LogSignature::concatenate_batch_coefficients`] to prefer the parallel
/// tournament reduction over the sequential fold chain for types the
/// cohort engine cannot carry (exact rationals and other non-raw-float
/// coefficients). Below this the per-fold/batch sequential path is cheaper
/// than the tournament's tree overhead (extra DAG warm-ups, pool
/// coordination, one final merge fold). Cohort-capable types ignore this
/// floor — their tournament engages from 2 displacements up, because the
/// SIMD-across-folds engine reuses one plan walk across 2-4 folds and
/// amortizes the tree overhead away.
pub(super) const TOURNAMENT_MIN_DISPLACEMENTS: usize = 16;

/// Minimum displacements per tournament leaf chunk. Leaves fold a
/// contiguous chunk of at least this many displacements through the
/// SEQUENTIAL batch driver (`fold_batch_sequential`), so a leaf's
/// support-growth warm-up folds and node-list rebuilds amortize over its
/// chunk exactly as they do for the sequential engine, instead of being
/// paid once per displacement. The floor (rather than a larger minimum)
/// keeps the leaf shape — and with it every fold's association order —
/// identical to the historical fixed 8-displacement-chunk tree for
/// batches of up to
/// `TOURNAMENT_MAX_LEAVES * TOURNAMENT_LEAF_CHUNK` displacements, which
/// is the regime the bit-level tournament tests pin.
pub(super) const TOURNAMENT_LEAF_CHUNK: usize = 8;

/// Upper bound on the tournament's leaf count: the adaptive chunk is
/// `rhss.len().div_ceil(TOURNAMENT_MAX_LEAVES)` raised to
/// [`TOURNAMENT_LEAF_CHUNK`], so a batch of `n` displacements runs at
/// most `TOURNAMENT_MAX_LEAVES` leaves. Each leaf now folds through the
/// warm-up + steady-batch sequential engine, so fewer-but-larger leaves
/// strictly amortize better: at 1000 displacements (the 2x12/3x8 e2e
/// path length) this is 32 warm-batched leaves instead of ~125
/// 8-displacement leaves that each paid cold support-growth rebuilds,
/// plus 31 merge folds instead of 124. ~32 in-flight leaves also keeps
/// oversubscription modest on the 16–32-thread pools the e2e workload
/// targets. The chunk is a pure function of `rhss.len()` — never of the
/// thread count — so the tree shape (and therefore results, bit for
/// bit) stays reproducible at any pool size.
pub(super) const TOURNAMENT_MAX_LEAVES: usize = 32;

/// NaN audit for an arbitrary coefficient slice — the
/// `LogSignature::audit_no_nan` check factored off `&self` so the
/// tournament's locally owned accumulators (plain `LieSeries` values, not
/// `LogSignature`s) can run the same per-fold audit.
#[inline]
pub(super) fn audit_coefficients_no_nan<U: PartialEq>(coefficients: &[U]) {
    if coefficients.iter().any(|c| is_nan_value(c)) {
        panic!("log-signature coefficients overflowed to NaN");
    }
}

/// One BCH fold of `rhs_coefficients` into `series` using `dag`: the exact
/// call sequence of `LogSignature::concatenate_assign_coefficients`,
/// factored out for the tournament runner, whose leaf and merge folds own
/// their accumulator/dag pair locally instead of going through a
/// `LogSignature`. Keeping the two in lockstep is what makes every
/// tournament fold semantically identical to the public per-fold path:
/// same support computation (with the degree cutoff), same DAG evaluation
/// (including the collecting rebuild whenever the atom supports differ from
/// the dag's stored ones), same term accumulation, same NaN audit. The
/// tournament's merge folds call this directly; a leaf's fold chain calls
/// it through `fold_batch_sequential`'s warm-up phase.
pub(super) fn fold_one_displacement<T, U>(
    series: &mut LieSeries<T, U>,
    dag: &mut CommutatorDag<U>,
    rhs_coefficients: &[U],
) where
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
{
    let original_coefficients = series.coefficients.clone();

    let a_nonzero = series
        .nonzero_coefficient_indices(&original_coefficients);
    let b_nonzero = series.nonzero_coefficient_indices(rhs_coefficients);

    dag.evaluate(
        series,
        &original_coefficients,
        &a_nonzero,
        rhs_coefficients,
        &b_nonzero,
    );

    // The DAG accumulated every term in class-contiguous order; this
    // single call applies the BCH weights and the one epilogue back to
    // public basis order.
    dag.accumulate_terms(&mut series.coefficients, rhs_coefficients);
    audit_coefficients_no_nan(&series.coefficients);
}

/// The eligibility check of [`try_fold_batch_rest`], factored out for the
/// cohort driver (which must test several lanes' readiness without folding
/// any of them): `Some((a_nonzero, b_nonzero))` when the whole `rest` can
/// run as ONE batch dispatch on the (`dag`, `acc`) pair, `None` otherwise.
/// `acc` is the accumulator's current coefficient slice (the caller's
/// series' coefficients — passed explicitly so cohort lanes can check
/// against a locally owned accumulator without building a series).
fn batch_rest_eligible<T, U>(
    series: &LieSeries<T, U>,
    acc: &[U],
    dag: &CommutatorDag<U>,
    rest: &[&[U]],
) -> Option<(Vec<usize>, Vec<usize>)>
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
        + AddAssign
        + Send
        + Sync
        + 'static,
{
    if rest.is_empty() {
        return None;
    }
    let a_nonzero = series.nonzero_coefficient_indices(acc);
    // All displacements must share one support — compared allocless (the
    // per-candidate `nonzero_coefficient_indices` Vec this replaces
    // allocated once per displacement, dominating the batch attempt's cost
    // on long batches).
    let b_nonzero = series.nonzero_coefficient_indices(rest[0]);
    if rest[1..].iter().any(|r| !series.has_support(r, &b_nonzero)) {
        return None;
    }
    // Eligibility: the node lists are the built fixed point for these
    // supports and the accumulator's support already equals the reachable
    // set — level-0 gating then uses masks that cover every position any
    // later fold can touch, which stays sound even if mid-batch
    // accumulation cancels values to zero (the node lists are
    // scatter-target supersets).
    if !dag.batch_eligible(&a_nonzero, &b_nonzero) {
        return None;
    }
    Some((a_nonzero, b_nonzero))
}

/// Attempts to batch ALL of `rhss` in one dispatch on the (`series`,
/// `dag`) pair; `false` means the batch contract does not hold yet
/// (sparse accumulator, displacements with differing supports, or lists
/// not steady) and the caller folds one displacement per-fold instead.
///
/// The body of `LogSignature`'s former `try_fold_batch_rest` method,
/// factored off `&self` onto a bare (`series`, `dag`) pair so the
/// tournament's leaf accumulators (plain `LieSeries` values with pooled
/// DAGs, not `LogSignature`s) can run the exact sequential batch driver
/// the public path runs. Keeping them in lockstep is what makes a leaf's
/// association order bit-identical to the sequential engine's on the
/// same chunk.
fn try_fold_batch_rest<T, U>(
    series: &mut LieSeries<T, U>,
    dag: &mut CommutatorDag<U>,
    rhss: &[&[U]],
) -> bool
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
        + AddAssign
        + Send
        + Sync
        + 'static,
{
    if rhss.is_empty() {
        return true;
    }
    let Some((a_nonzero, b_nonzero)) =
        batch_rest_eligible(series, &series.coefficients, dag, rhss)
    else {
        return false;
    };
    dag.fold_batch(series, rhss, &a_nonzero, &b_nonzero);
    audit_coefficients_no_nan(&series.coefficients);
    true
}

/// The strictly sequential fold chain on a bare (`series`, `dag`) pair:
/// folds one-by-one until the accumulator's support is the full basis
/// and the node lists are steady for it, then dispatches the remaining
/// folds as ONE continuous batch.
///
/// This is the historical body of
/// `LogSignature::concatenate_batch_sequential`, factored off `&self`
/// for the tournament's leaves: each leaf folds its displacement chunk
/// from a zeroed accumulator through this exact engine, so the leaf
/// gets the same warm-up + steady-batch treatment (and the same
/// association order) as the public sequential path, with the chunk's
/// support-growth rebuilds amortized over the whole chunk instead of
/// paid per displacement.
pub(super) fn fold_batch_sequential<T, U>(
    series: &mut LieSeries<T, U>,
    dag: &mut CommutatorDag<U>,
    rhss: &[&[U]],
) where
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
{
    let mut i = 0usize;
    while i < rhss.len() {
        // Batch the remaining folds when eligible; otherwise fold one
        // — growing support or rebuilding lists exactly as the
        // per-fold path does — and re-check.
        if !try_fold_batch_rest(series, dag, &rhss[i..]) {
            fold_one_displacement(series, dag, rhss[i]);
            i += 1;
        } else {
            return;
        }
    }
}

/// A `LieSeries` template shared read-only into rayon tasks regardless of
/// `T`'s auto traits (the public batch driver carries no `T: Send + Sync`
/// bound, so `T`-bearing references cannot cross into tasks directly).
///
/// # Safety contract
///
/// The `unsafe Send`/`Sync` impls rely on the tournament never touching
/// `T` data off the calling thread:
///
/// 1. `LieSeries::clone` — the only operation tasks perform on the
///    template — is three `Arc` refcount bumps plus a `Vec<U>` copy; it
///    never reads `T` data.
/// 2. The fold path (`nonzero_coefficient_indices`, `evaluate`,
///    `accumulate_terms`) reads only the `U` coefficients, `basis.len()`
///    and the `U`/`u8` decomposition tables; the `T`-bearing basis and
///    commutator-term arrays are never dereferenced, let alone mutated.
/// 3. The template outlives the whole parallel region (rayon joins before
///    `tournament_reduce` returns) and `self.series` outlives the
///    template, so a task's series clone dropped on a worker thread never
///    releases the last `Arc` reference of a `T`-bearing allocation — no
///    `T` destructor ever runs off-thread.
pub(super) struct SeriesTemplate<U, T>(pub(super) LieSeries<T, U>);

// SAFETY: see the struct-level safety contract.
unsafe impl<U: Send + Sync, T> Send for SeriesTemplate<U, T> {}
// SAFETY: see the struct-level safety contract.
unsafe impl<U: Send + Sync, T> Sync for SeriesTemplate<U, T> {}

/// Pops a pooled DAG, falling back to a fresh shallow clone of `source`
/// (the caller's compiled plan) when the pool is empty.
pub(super) fn take_pooled_dag<U>(
    pool: &Mutex<Vec<CommutatorDag<U>>>,
    source: &CommutatorDag<U>,
) -> CommutatorDag<U> {
    pool.lock()
        .unwrap_or_else(|poisoned| poisoned.into_inner())
        .pop()
        .unwrap_or_else(|| source.clone_shallow())
}

/// Returns a finished DAG to the pool for the next fold.
pub(super) fn return_pooled_dag<U>(pool: &Mutex<Vec<CommutatorDag<U>>>, dag: CommutatorDag<U>) {
    pool.lock()
        .unwrap_or_else(|poisoned| poisoned.into_inner())
        .push(dag);
}

/// True when every displacement slice in the leaf group's chunks has
/// kernel-visible support ⊆ `b_nonzero` (the shared-plan condition of
/// `leaf_group_engine`). Subset — not equality — is the sound condition:
/// the shared plan's masks are derived from the LARGEST support, and a
/// smaller-support slice gathers exact zeros into every unshared mask
/// position, contributing exact-zero adds (bit-identical to the scalar
/// engine skipping them — the sign of zero is unobservable on
/// accumulators).
fn leaf_group_supports_subset<T, U>(
    template: &SeriesTemplate<U, T>,
    rhss: &[&[U]],
    chunk: usize,
    first_leaf: usize,
    group_leaves: usize,
    b_nonzero: &[usize],
) -> bool
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
        + AddAssign
        + Send
        + Sync
        + 'static,
{
    (0..group_leaves).all(|l| {
        let start = (first_leaf + l) * chunk;
        let end = (start + chunk).min(rhss.len());
        rhss[start..end].iter().all(|r| {
            template
                .0
                .nonzero_coefficient_indices(r)
                .iter()
                .all(|&i| b_nonzero.binary_search(&i).is_ok())
        })
    })
}

/// Derives a leaf's reachable-support plan: the fixed point
/// `a* = union(lists(a*, b)) ∪ b` of the DAG's scatter-target sets, in
/// public basis order. See `cohort_leaf_group` for why the per-node lists
/// must be rebuilt for exactly these atom supports before any shared-plan
/// dispatch (the interior jobs' operand gating reads them verbatim;
/// lists built for smaller atom supports miss the DAG's acc-side units
/// and their targets, dropping real contributions once the accumulators
/// grow past the displacement support). The iteration is monotone and
/// bounded by the basis: two passes in practice — the first build already
/// reaches the full word closure (every closure word has a pure-B
/// decomposition chain), the second repairs the per-node target sets for
/// the acc-side units.
pub(super) fn leaf_steady_plan<T, U>(
    template: &SeriesTemplate<U, T>,
    dag: &mut CommutatorDag<U>,
    zero_acc: &[U],
    first_slice: &[U],
    b: &[usize],
) -> Vec<usize>
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
        + AddAssign
        + Send
        + Sync
        + 'static,
{
    let mut a_plan = b.to_vec();
    loop {
        dag.ensure_lists_steady(&template.0, zero_acc, &a_plan, first_slice, b);
        let order = template.0.class_order();
        let perm = order.perm();
        let mut next: Vec<usize> = dag
            .nonzeros
            .iter()
            .flat_map(|l| l.iter().map(|&p| perm[p] as usize))
            .collect();
        next.extend_from_slice(b);
        next.sort_unstable();
        next.dedup();
        let converged = next == a_plan;
        a_plan = next;
        if converged {
            return a_plan;
        }
    }
}

/// The group's steady plan, via the process-wide leaf plan cache: a cache
/// hit adopts the snapshot's fixed-point lists into `dag` (skipping both
/// the fixed-point loop and its collecting kernel sweeps — see
/// [`leaf_steady_cache`]); a miss derives on `dag` exactly as
/// [`leaf_steady_plan`] does and stores the snapshot for later groups.
/// Adopted and derived states are bit-identical (collection is
/// support-determined), so the engine's fold is the same either way.
pub(super) fn leaf_steady_plan_cached<T, U>(
    template: &SeriesTemplate<U, T>,
    dag: &mut CommutatorDag<U>,
    zero_acc: &[U],
    first_slice: &[U],
    b: &[usize],
) -> Vec<usize>
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
        + AddAssign
        + Send
        + Sync
        + 'static,
{
    let table_id = template.0.feasible_decompositions_handle().table_id();
    if let Some(snap) = leaf_steady_cache::get(table_id, b)
        && dag.adopt_steady_lists(
            &snap.a_plan,
            b,
            Arc::clone(&snap.nonzeros),
            Arc::clone(&snap.dirty),
            Arc::clone(&snap.order),
        )
    {
        return snap.a_plan.clone();
    }
    let a_plan = leaf_steady_plan(template, dag, zero_acc, first_slice, b);
    let (nonzeros, dirty) = dag.steady_lists_snapshot();
    leaf_steady_cache::store(
        table_id,
        b,
        a_plan.clone(),
        nonzeros,
        dirty,
        template.0.class_order(),
    );
    a_plan
}

/// Fills `dag`'s steady-list state for `(a_plan, b)` from the leaf plan
/// cache when possible (adopt — no collecting sweep), falling back to the
/// regular collecting rebuild. Used by the scalar leaf fold, whose dag is
/// popped from the pool AFTER the group's plans were derived: the fold
/// must run on the derived fixed-point lists, and adoption provides them
/// without a rebuild regardless of which pooled dag the lane receives.
/// Returns `true` when the lists were adopted (the caller skips the
/// rebuild); `false` means the caller must call `ensure_lists_steady`.
pub(super) fn leaf_steady_adopt_into<U>(
    dag: &mut CommutatorDag<U>,
    table_id: u64,
    a_plan: &[usize],
    b: &[usize],
) -> bool
where
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
{
    match leaf_steady_cache::get(table_id, b) {
        Some(snap) if snap.a_plan == *a_plan => dag.adopt_steady_lists(
            &snap.a_plan,
            b,
            Arc::clone(&snap.nonzeros),
            Arc::clone(&snap.dirty),
            Arc::clone(&snap.order),
        ),
        _ => false,
    }
}

/// The tournament leaf round's group engine: folds `group_leaves` (2..=4)
/// adjacent leaf chunks, each from a zeroed accumulator, through the
/// STEADY-CHUNK leaf engine — one reachable-support plan per lane (see
/// [`leaf_steady_plan`]), then the lane's WHOLE chunk as ONE batch
/// dispatch — no per-fold warm-up walks anywhere.
///
/// When the group's displacement slices are support-uniform (every slice's
/// kernel-visible support ⊆ lane 0's first slice's support — the subset
/// test is the sound condition: unshared mask positions gather exact
/// zeros and contribute exact-zero adds) and `cohort` is on, the lanes
/// share ONE plan and ONE dispatch: the SoA cohort engine
/// ([`CommutatorDag::fold_batch_cohort`]) runs every lane's whole chunk
/// with the lanes' values as 4-lane vectors — up to 4 folds per plan walk
/// from the very first fold of the batch. Otherwise every lane folds
/// alone through the scalar batch engine ([`CommutatorDag::fold_batch`])
/// under its own plan (the group's lanes in parallel tasks).
///
/// # Bit-identity
///
/// The plan (decomposition lists, scatter indices, gating, BCH weights)
/// is a pure function of the basis and the atom supports, and the
/// iteration derives the same reachable support a* and the same
/// fixed-point node lists in every engine state. `fold_batch_cohort`
/// replicates `fold_batch`'s per-fold operation order exactly (the
/// oracle test pins the two engines bit-identical at shared masks), so
/// the cohort and scalar leaves of the SAME tree are bit-identical by
/// construction — which is what makes the tournament's cohort-vs-scalar
/// oracle same-tree sharp. The engine is NOT bit-identical to the
/// historical warm-up-then-batch leaf walk (the dense plan's gating
/// tickets order the shared units differently than the sparse per-fold
/// gating, moving a few ulps on adversarial data); the tournament's
/// reassociation caveat covers leaf-internal ordering, and the same
/// engine runs on both sides of the kill switch, so results stay
/// reproducible bit for bit at any pool size and either switch state.
pub(super) fn leaf_group_engine<T, U>(
    template: &SeriesTemplate<U, T>,
    source_dag: &CommutatorDag<U>,
    dag_pool: &Mutex<Vec<CommutatorDag<U>>>,
    basis_len: usize,
    rhss: &[&[U]],
    chunk: usize,
    first_leaf: usize,
    group_leaves: usize,
    cohort: bool,
) -> Vec<Vec<U>>
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
        + AddAssign
        + Send
        + Sync
        + 'static,
{
    use rayon::prelude::*;

    debug_assert!((1..=4).contains(&group_leaves), "leaf groups fold 1..=4 lanes");
    let first = first_leaf * chunk;
    let b0 = template.0.nonzero_coefficient_indices(rhss[first]);
    // Support-uniformity of the group's displacement slices — a property of
    // the DATA alone, deliberately independent of the cohort kill switch:
    // both switch states must derive the SAME plan for a uniform group, or
    // the shared-plan cohort dispatch and the per-lane scalar fallback
    // would order the shared units' contributions differently (the exact
    // bit-identity the kill-switch oracle tests pin).
    let uniform_support = group_leaves >= 2
        && leaf_group_supports_subset(template, rhss, chunk, first_leaf, group_leaves, &b0);
    let uniform = cohort && uniform_support;
    if uniform {
        // Shared plan, one SoA dispatch for every lane's whole chunk. The
        // probe dag derives the plan (and ends up holding the fixed-point
        // lists) and then serves as the cohort's lead dag.
        let zero = vec![U::default(); basis_len];
        let mut lead = take_pooled_dag(dag_pool, source_dag);
        let a_plan = leaf_steady_plan_cached(template, &mut lead, &zero, rhss[first], &b0);
        let mut cohort_lanes: Vec<CohortLane<U>> = (0..group_leaves)
            .map(|l| {
                let start = (first_leaf + l) * chunk;
                CohortLane {
                    acc: vec![U::default(); basis_len],
                    rhss: &rhss[start..(start + chunk).min(rhss.len())],
                }
            })
            .collect();
        lead.fold_batch_cohort(&template.0, &mut cohort_lanes, &a_plan, &b0);
        return_pooled_dag(dag_pool, lead);
        return cohort_lanes
            .into_iter()
            .map(|lane| {
                audit_coefficients_no_nan(&lane.acc);
                lane.acc
            })
            .collect();
    }
    // Per-lane scalar batch under each lane's own steady plan (the group's
    // lanes in parallel tasks; no shared plan to cohort with). Support-
    // uniform groups (cohort off) reuse the shared (b0, a*(b0)) plan so
    // the kill-switch runs are bit-identical to the cohort ones; divergent
    // groups derive each lane's plan from the lane's own displacement
    // supports.
    // Per-lane plans. Support-uniform groups (cohort off) reuse the shared
    // (b0, a*(b0)) plan so the kill-switch runs are bit-identical to the
    // cohort ones; divergent groups derive each lane's plan from the
    // lane's own displacement-support union.
    let plans: Vec<(Vec<usize>, Vec<usize>)> = if uniform_support {
        let zero = vec![U::default(); basis_len];
        let mut probe = take_pooled_dag(dag_pool, source_dag);
        let a = leaf_steady_plan_cached(template, &mut probe, &zero, rhss[first], &b0);
        return_pooled_dag(dag_pool, probe);
        vec![(a, b0.clone()); group_leaves]
    } else {
        let zero = vec![U::default(); basis_len];
        // Per-lane plans, derived in PARALLEL on per-lane dags. Each plan's
        // fixed-point loop runs two collecting rebuilds of the whole DAG
        // (serial node-by-node kernel calls — at 4x8 the dominating cost of
        // the leaf round); deriving the group's plans on ONE probe dag
        // serialized those chains, and at wide pools the serial planning —
        // not the parallel dispatch — set the round's wall time. The plans
        // are pure functions of (template, zero, first slice, union
        // support): the dag they are derived on is interchangeable, so the
        // parallel form is bit-identical.
        let lanes: Vec<CommutatorDag<U>> = (0..group_leaves)
            .map(|_| take_pooled_dag(dag_pool, source_dag))
            .collect();
        let plans: Vec<(Vec<usize>, Vec<usize>)> = lanes
            .into_par_iter()
            .enumerate()
            .map(|(l, mut dag)| {
                let start = (first_leaf + l) * chunk;
                let end = (start + chunk).min(rhss.len());
                // The chunk's union displacement support. Collected then
                // sorted+deduped: an incremental binary_search-insert into the
                // growing (unsorted) vector can miss contained values and keep
                // duplicates, and duplicate supports corrupt the batch engine's
                // displacement AccumRuns (each duplicate emits its own run — the
                // displacement is added once per duplicate, a real value bug the
                // exact-arithmetic oracle tests catch).
                let mut b: Vec<usize> = rhss[start..end]
                    .iter()
                    .flat_map(|r| template.0.nonzero_coefficient_indices(r))
                    .collect();
                b.sort_unstable();
                b.dedup();
                let a = leaf_steady_plan_cached(template, &mut dag, &zero, rhss[start], &b);
                // The plan dag returns to the pool already holding this
                // lane's atom supports and fixed-point lists, so the lane
                // fold's own ensure_lists_steady early-returns.
                return_pooled_dag(dag_pool, dag);
                (a, b)
            })
            .collect();
        plans
    };
    (0..group_leaves)
        .into_par_iter()
        .map(|l| {
            let start = (first_leaf + l) * chunk;
            let end = (start + chunk).min(rhss.len());
            let (a_plan, b) = (&plans[l].0, &plans[l].1);
            let mut series = template.0.clone();
            series.coefficients = vec![U::default(); basis_len];
            let mut dag = take_pooled_dag(dag_pool, source_dag);
            // The fold must run on the derived fixed-point lists. The plan
            // cache provides them by adoption (no collecting sweep)
            // regardless of which pooled dag the lane receives; adoption
            // installs the same class order and atom supports a fresh
            // derivation would, so the subsequent fold is bit-identical to
            // the post-rebuild fold.
            let table_id = template.0.feasible_decompositions_handle().table_id();
            if !leaf_steady_adopt_into(&mut dag, table_id, a_plan, b) {
                dag.ensure_lists_steady(&template.0, &series.coefficients, a_plan, rhss[start], b);
            }
            dag.fold_batch(&mut series, &rhss[start..end], a_plan, b);
            return_pooled_dag(dag_pool, dag);
            audit_coefficients_no_nan(&series.coefficients);
            series.coefficients
        })
        .collect()
}

/// The tournament merge round's COHORT group: folds 2-4 adjacent merge
/// pairs as one 4-lane cohort. Each lane folds ONE displacement (the right
/// chunk result) into its dense accumulator (the left chunk result) —
/// today's per-lane merge fold, with the four folds sharing one plan walk.
///
/// Eligibility: all lanes' (accumulator, displacement) supports must be
/// identical (`dense x dense` merges always are — the leaf results' support
/// is the same support-derived reachable set on every lane). The shared
/// DAG's node lists are brought to the merge supports with ONE collecting
/// pass ([`CommutatorDag::ensure_lists_steady`], value-independent — the
/// collected targets are determined by the supports' gating alone) and the
/// cohort engine runs the union plan; per-lane gating is vacuous for
/// identical supports. Divergent supports fall back to the scalar per-lane
/// merge fold.
pub(super) fn cohort_merge_group<T, U>(
    template: &SeriesTemplate<U, T>,
    source_dag: &CommutatorDag<U>,
    dag_pool: &Mutex<Vec<CommutatorDag<U>>>,
    level: &[Vec<U>],
    pairs: &[usize],
) -> Vec<Vec<U>>
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
        + AddAssign
        + Send
        + Sync
        + 'static,
{
    use rayon::prelude::*;

    // Per-lane supports (value-dependent — computed from the values).
    let supports: Vec<(Vec<usize>, Vec<usize>)> = pairs
        .iter()
        .map(|&p| {
            (
                template
                    .0
                    .nonzero_coefficient_indices(&level[2 * p]),
                template
                    .0
                    .nonzero_coefficient_indices(&level[2 * p + 1]),
            )
        })
        .collect();
    if !supports.iter().all(|s| s == &supports[0]) {
        // Divergent supports: the scalar merge fold per pair, nested-
        // parallel (each task builds its own series — the template safety
        // contract covers the off-thread drops).
        return pairs
            .par_iter()
            .map(|&p| scalar_merge_fold(template, source_dag, dag_pool, &level[2 * p], &level[2 * p + 1]))
            .collect();
    }
    let (a_nonzero, b_nonzero) = (&supports[0].0, &supports[0].1);
    let mut dag = take_pooled_dag(dag_pool, source_dag);
    dag.ensure_lists_steady(
        &template.0,
        &level[2 * pairs[0]],
        a_nonzero,
        &level[2 * pairs[0] + 1],
        b_nonzero,
    );
    let rhs_views: Vec<&[U]> = pairs.iter().map(|&p| level[2 * p + 1].as_slice()).collect();
    let mut cohort_lanes: Vec<CohortLane<U>> = pairs
        .iter()
        .zip(&rhs_views)
        .map(|(&p, rhs)| CohortLane {
            acc: level[2 * p].clone(),
            rhss: std::slice::from_ref(rhs),
        })
        .collect();
    dag.fold_batch_cohort(&template.0, &mut cohort_lanes, a_nonzero, b_nonzero);
    return_pooled_dag(dag_pool, dag);
    cohort_lanes
        .into_iter()
        .map(|lane| {
            audit_coefficients_no_nan(&lane.acc);
            lane.acc
        })
        .collect()
}

/// The tournament's scalar merge fold (today's per-pair merge task body):
/// one `fold_one_displacement` of `right` into a `left`-seeded accumulator.
pub(super) fn scalar_merge_fold<T, U>(
    template: &SeriesTemplate<U, T>,
    source_dag: &CommutatorDag<U>,
    dag_pool: &Mutex<Vec<CommutatorDag<U>>>,
    left: &[U],
    right: &[U],
) -> Vec<U>
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
        + AddAssign
        + Send
        + Sync
        + 'static,
{
    let mut series = template.0.clone();
    series.coefficients = left.to_vec();
    let mut dag = take_pooled_dag(dag_pool, source_dag);
    fold_one_displacement(&mut series, &mut dag, right);
    return_pooled_dag(dag_pool, dag);
    series.coefficients
}
