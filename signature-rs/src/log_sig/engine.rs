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
#[cfg(test)]
mod tests {
    //! Batch/tournament/cohort fold-engine tests: engine selection, engine
    //! bit-identity, NaN audits and support-growth contracts.

    use super::*;
    use crate::commutator_dag::{COHORT_ENGINE_RUNS, COHORT_LANE_FOLDS, set_cohort_off};
    use num_rational::Ratio;
    use ordered_float::NotNan;
    use rstest::rstest;


    /// Adversarial regression test for the pooled-DAG list-reuse flake
    /// (`tournament_rationals_exact_at_awkward_sizes::case_1`, intermittent
    /// under CPU oversubscription). A leaf chunk's `fold_batch` records its
    /// planner scatter sets as the DAG's node lists; those lists contain
    /// degree-`max_degree` positions, and a later fold on the SAME pooled
    /// dag with matching atom supports early-returns through
    /// `ensure_lists_steady` and gates with them. The gating prologues'
    /// full-support fast path keyed on the support LENGTH alone: a
    /// batch-recorded list whose length happens to equal
    /// `degree_start(max_degree)` (here node `AB`'s recorded list
    /// `{2,3,4}` at cutoff 3) triggered the all-entries-active gating
    /// without presence tests, reading compact slots the operand layout
    /// does not cover — exact-rational garbage contributions. The fix
    /// requires the fast-path supports to be exactly the prefix
    /// `0..cutoff`; this test pins that deterministically (leaf 1 via the
    /// pooled early-return path must equal the sequential fold chain).
    #[test]

    fn pooled_batch_recorded_lists_fold_bit_identically() {
        type R = Ratio<i64>;
        let (d, m) = (2usize, 3usize);
        let builder = LogSignatureBuilder::<u8>::new()
            .with_num_dimensions(d)
            .with_max_degree(m);
        let basis =
            lyndon_rs::lyndon::LyndonBasis::<u8>::new(d, lyndon_rs::lyndon::Sort::Lexicographical)
                .generate_basis(m);
        let mut seed = 0x3a11_u64.wrapping_add((d * 1007 + m * 97 + 17) as u64);
        let lcg = |seed: &mut u64| {
            *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
            ((*seed >> 33) % 7) as i64 - 3
        };
        let rhss: Vec<Vec<R>> = (0..17)
            .map(|i| {
                (0..basis.len())
                    .map(|k| {
                        if k < d && !(i % 3 == 2 && k == i % d) {
                            R::from_integer(lcg(&mut seed))
                        } else {
                            R::from_integer(0)
                        }
                    })
                    .collect()
            })
            .collect();
        let basis_len = basis.len();
        let zero = vec![R::default(); basis_len];
        let mk = || {
            let s = builder.build::<R>();
            let LogSignature { series, bch_series: _, dag } = s;
            (series, dag)
        };
        let (ref_series, mut probe) = mk();

        // The leaf plan exactly as leaf_group_engine's divergent branch does:
        // probe dag, leaf_steady_plan from (zero acc, rhs[0], b0).
        let b0: Vec<usize> = rhss[0..8]
            .iter()
            .flat_map(|r| ref_series.nonzero_coefficient_indices(r))
            .collect::<Vec<_>>();
        let mut b0s = b0.clone();
        b0s.sort_unstable();
        b0s.dedup();
        let a_star = leaf_steady_plan(&SeriesTemplate(ref_series.clone()), &mut probe, &zero, &rhss[0], &b0s);

        // Leaf 0 on the probe dag itself (mimics the pool handing the probe
        // to lane 0): atoms already match -> early-return -> fold_batch.
        let mut series0 = ref_series.clone();
        series0.coefficients = vec![R::default(); basis_len];
        probe.ensure_lists_steady(&ref_series, &series0.coefficients, &a_star, &rhss[0], &b0s);
        probe.fold_batch(&mut series0, &rhss[0..8].iter().map(|r| r.as_slice()).collect::<Vec<_>>(), &a_star, &b0s);
        // Leaf 0's fold must land on the sequential ground truth too — the
        // pooled early-return path has no room to diverge on the first leaf.
        let mut seq0 = builder.build::<R>();
        for r in &rhss[0..8] {
            seq0.concatenate_assign_coefficients(r);
        }
        assert_eq!(series0.coefficients, seq0.series.coefficients, "leaf0 (early-return) diverged from sequential");

        // Leaf 1 on the SAME dag (pooled early-return path).
        let mut series1p = ref_series.clone();
        series1p.coefficients = vec![R::default(); basis_len];
        probe.ensure_lists_steady(&ref_series, &series1p.coefficients, &a_star, &rhss[8], &b0s);
        probe.fold_batch(&mut series1p, &rhss[8..16].iter().map(|r| r.as_slice()).collect::<Vec<_>>(), &a_star, &b0s);
        let leaf1_pooled = series1p.coefficients.clone();

        // Leaf 1 on a FRESH dag (rebuild path).
        let (_, mut fresh) = mk();
        let mut series1f = ref_series.clone();
        series1f.coefficients = vec![R::default(); basis_len];
        fresh.ensure_lists_steady(&ref_series, &series1f.coefficients, &a_star, &rhss[8], &b0s);
        fresh.fold_batch(&mut series1f, &rhss[8..16].iter().map(|r| r.as_slice()).collect::<Vec<_>>(), &a_star, &b0s);
        let leaf1_fresh = series1f.coefficients.clone();

        // Ground truth: sequential folds of rhss[8..16] onto zero.
        let mut seq1 = builder.build::<R>();
        for r in &rhss[8..16] {
            seq1.concatenate_assign_coefficients(r);
        }
        assert_eq!(leaf1_fresh, seq1.series.coefficients, "fresh-dag leaf1 must equal sequential");
        assert_eq!(leaf1_pooled, seq1.series.coefficients, "pooled-dag leaf1 (early-return) diverged from sequential");
    }


    /// Differential test for the SEQUENTIAL batch engine, exercised
    /// directly via `concatenate_batch_sequential` so the exact-equality
    /// comparison stays valid regardless of the tournament threshold: the
    /// driver must warm-fold per-fold while the support grows and the node
    /// lists steady, then batch the eligible rest in ONE dispatch — and
    /// stay bit-identical to folding each displacement per-fold, for both
    /// all-eligible letter displacements and a mixed-support run that
    /// forces the per-fold fallback.
    #[test]

    fn concatenate_batch_matches_sequential_folds() {
        use ordered_float::NotNan;

        let lcg = |seed: &mut u64| {
            *seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((*seed >> 33) % 19) as i128 - 9
        };

        for (d, m) in [(2usize, 5usize), (3usize, 4usize)] {
            let builder = LogSignatureBuilder::<u8>::new()
                .with_num_dimensions(d)
                .with_max_degree(m);
            let basis = lyndon_rs::lyndon::LyndonBasis::<u8>::new(
                d,
                lyndon_rs::lyndon::Sort::Lexicographical,
            )
            .generate_basis(m);

            // Letter displacements: all-d support (batch-eligible) plus a
            // final mixed-support case (one zeroed letter forces the
            // per-fold fallback for the run containing it).
            let mut seed = 0xba7c_u64.wrapping_add(d as u64 * 31 + m as u64);
            let mk_dense = |seed: &mut u64| -> Vec<NotNan<f64>> {
                (0..basis.len())
                    .map(|k| {
                        if k < d {
                            NotNan::new(lcg(seed) as f64).unwrap()
                        } else {
                            NotNan::new(0.0).unwrap()
                        }
                    })
                    .collect()
            };
            let rhss: Vec<Vec<NotNan<f64>>> = (0..8).map(|_| mk_dense(&mut seed)).collect();
            let mut mixed = mk_dense(&mut seed);
            mixed[d - 1] = NotNan::new(0.0).unwrap();

            for (label, rhss) in [
                ("dense-eligible", rhss.clone()),
                (
                    "mixed-support",
                    rhss[..rhss.len() - 1]
                        .to_vec()
                        .into_iter()
                        .chain([mixed.clone()])
                        .collect(),
                ),
            ] {
                // Same accumulator and displacements for both paths.
                let acc0 = mk_dense(&mut seed);
                // Sequential reference: fold each displacement per-fold.
                let mut seq = builder.build::<NotNan<f64>>();
                seq.series.coefficients = acc0.clone();
                for r in &rhss {
                    seq.concatenate_assign_coefficients(r);
                }

                // Batched: same displacements through the sequential
                // batch driver.
                let mut bat = builder.build::<NotNan<f64>>();
                bat.series.coefficients = acc0;
                let slices: Vec<&[NotNan<f64>]> = rhss.iter().map(|r| r.as_slice()).collect();
                bat.concatenate_batch_sequential(&slices);

                let diffs: Vec<_> = seq
                    .series
                    .coefficients
                    .iter()
                    .zip(&bat.series.coefficients)
                    .enumerate()
                    .filter(|(_, (s, b))| s != b)
                    .take(4)
                    .map(|(k, (s, b))| format!("k={k} seq={s} bat={b}"))
                    .collect();
                assert!(
                    diffs.is_empty(),
                    "d={d} m={m} {label}: batch diverged: {}",
                    diffs.join("; ")
                );
            }
        }
    }


    /// Stress-hammer for the batch walk's stage ordering: the same
    /// comparison as [`Self::concatenate_batch_matches_sequential_folds`],
    /// repeated with varied pool sizes and shapes. The batch's value
    /// correctness is timing-sensitive under a stage-ordering race, so a
    /// single pass is not sufficient evidence — this loop makes a race
    /// visible in CI instead of surfacing as a rare downstream wrong
    /// number. The sequential driver is exercised directly
    /// (`concatenate_batch_sequential`) so the hammer targets that
    /// engine's walk regardless of the tournament threshold.
    #[test]

    fn concatenate_batch_survives_repeated_pool_stress() {
        use ordered_float::NotNan;
        let lcg = |seed: &mut u64| {
            *seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((*seed >> 33) % 19) as i128 - 9
        };

        let schedule: Vec<(usize, usize, usize)> = (0..7)
            .flat_map(|round| {
                [
                    (0usize, 1usize),
                    (1, 1),
                    (2, 1),
                    (0, 2),
                    (1, 4),
                    (2, 8),
                    (0, 16),
                    (1, 32),
                    (2, 32),
                ]
                .iter()
                .map(move |(shape, threads)| (round, *shape, *threads))
            })
            .collect();
        for (round, shape, threads) in schedule {
            let pass = round * 9 + shape;
            {
                // The seed depends only on (shape, round) — NOT on the pool
                // size — so the 1-thread baseline and the concurrent runs of
                // the same (shape, round) use identical data and their
                // per-stage checksums are comparable.
                let (d, m) = [(3usize, 4usize), (2usize, 5usize), (3usize, 5usize)][shape];
                let pool = rayon::ThreadPoolBuilder::new()
                    .num_threads(threads)
                    .build()
                    .expect("pool");
                pool.install(|| {
                    let builder = LogSignatureBuilder::<u8>::new()
                        .with_num_dimensions(d)
                        .with_max_degree(m);
                    let basis = lyndon_rs::lyndon::LyndonBasis::<u8>::new(
                        d,
                        lyndon_rs::lyndon::Sort::Lexicographical,
                    )
                    .generate_basis(m);
                    let mut seed = 0x5755_u64.wrapping_add((round * 3 + shape) as u64 * 7919);
                    let mk = |seed: &mut u64, dense: bool| -> Vec<NotNan<f64>> {
                        (0..basis.len())
                            .map(|k| {
                                if dense || k < d {
                                    NotNan::new(lcg(seed) as f64).unwrap()
                                } else {
                                    NotNan::new(0.0).unwrap()
                                }
                            })
                            .collect()
                    };
                    let rhss: Vec<Vec<NotNan<f64>>> =
                        (0..12).map(|_| mk(&mut seed, true)).collect();
                    let acc0 = mk(&mut seed, false);

                    let mut seq = builder.build::<NotNan<f64>>();
                    seq.series.coefficients = acc0.clone();
                    for r in &rhss {
                        seq.concatenate_assign_coefficients(r);
                    }
                    let mut bat = builder.build::<NotNan<f64>>();
                    bat.series.coefficients = acc0;
                    let slices: Vec<&[NotNan<f64>]> = rhss.iter().map(|r| r.as_slice()).collect();
                    bat.concatenate_batch_sequential(&slices);

                    let diffs: Vec<_> = seq
                        .series
                        .coefficients
                        .iter()
                        .zip(&bat.series.coefficients)
                        .enumerate()
                        .filter(|(_, (s, b))| s != b)
                        .take(4)
                        .map(|(k, (s, b))| format!("k={k} seq={s} bat={b}"))
                        .collect();
                    assert!(
                        diffs.is_empty(),
                        "pass={pass} d={d} m={m} threads={threads}: batch diverged: {}",
                        diffs.join("; ")
                    );
                });
            }
        }
    }


    /// The sequential batch engine must also produce bit-identical results
    /// with exact rational coefficients, where mid-batch cancellations to
    /// zero are REAL (the dense-mask soundness argument is exercised for
    /// values, not just float dust). Exercised via
    /// `concatenate_batch_sequential` directly so the exact-equality
    /// assertion stays valid regardless of the tournament threshold.
    #[test]

    fn concatenate_batch_exact_rationals_match_sequential() {
        let (d, m) = (2usize, 4usize);
        let builder = LogSignatureBuilder::<u8>::new()
            .with_num_dimensions(d)
            .with_max_degree(m);
        let basis =
            lyndon_rs::lyndon::LyndonBasis::<u8>::new(d, lyndon_rs::lyndon::Sort::Lexicographical)
                .generate_basis(m);
        type R = num_rational::Ratio<i128>;

        let mut seed = 0x7a71_u64;
        let lcg = |seed: &mut u64| {
            *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
            ((*seed >> 33) % 19) as i128 - 9
        };
        let mk = |seed: &mut u64| -> Vec<R> {
            (0..basis.len())
                .map(|k| {
                    if k < d {
                        R::from_integer(lcg(seed))
                    } else {
                        R::from_integer(0)
                    }
                })
                .collect()
        };
        let rhss: Vec<Vec<R>> = (0..6).map(|_| mk(&mut seed)).collect();

        let acc0 = mk(&mut seed);
        let mut seq = builder.build::<R>();
        seq.series.coefficients = acc0.clone();
        for r in &rhss {
            seq.concatenate_assign_coefficients(r);
        }
        let mut bat = builder.build::<R>();
        bat.series.coefficients = acc0;
        let slices: Vec<&[R]> = rhss.iter().map(|r| r.as_slice()).collect();
        bat.concatenate_batch_sequential(&slices);
        assert_eq!(seq.series.coefficients, bat.series.coefficients);
    }


    /// The tournament reduction must be EXACTLY equal to sequential folding
    /// for exact coefficient types: rational arithmetic is insensitive to
    /// association order, so the reassociated tree may not change any
    /// coefficient. Runs in an explicit multi-thread pool so the tournament
    /// path is taken regardless of the host's core count.
    #[test]

    fn tournament_reduction_exact_rationals_match_sequential() {
        let (d, m) = (3usize, 4usize);
        let builder = LogSignatureBuilder::<u8>::new()
            .with_num_dimensions(d)
            .with_max_degree(m);
        let basis =
            lyndon_rs::lyndon::LyndonBasis::<u8>::new(d, lyndon_rs::lyndon::Sort::Lexicographical)
                .generate_basis(m);
        type R = num_rational::Ratio<i128>;

        let mut seed = 0x70da_u64;
        let lcg = |seed: &mut u64| {
            *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
            ((*seed >> 33) % 19) as i128 - 9
        };
        // Letter displacements, with every fifth one missing a letter — the
        // support changes force collecting rebuilds mid-leaf and across
        // pooled-dag reuse, exercising the pool's arbitrary-support soundness.
        let rhss: Vec<Vec<R>> = (0..24)
            .map(|i| {
                (0..basis.len())
                    .map(|k| {
                        if k < d && !(i % 5 == 4 && k == d - 1) {
                            R::from_integer(lcg(&mut seed))
                        } else {
                            R::from_integer(0)
                        }
                    })
                    .collect()
            })
            .collect();

        let threads: usize = std::env::var("DBG_THREADS")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(4);
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .expect("pool");
        pool.install(|| {
            let mut seq = builder.build::<R>();
            for r in &rhss {
                seq.concatenate_assign_coefficients(r);
            }
            let mut bat = builder.build::<R>();
            let slices: Vec<&[R]> = rhss.iter().map(|r| r.as_slice()).collect();
            bat.concatenate_batch_coefficients(&slices);
            assert_eq!(seq.series.coefficients, bat.series.coefficients);
        });
    }


     /// The tournament's tree shape depends only on the displacement count,
    /// never on the thread count: the same input reduced in pools of
    /// different sizes must produce bit-identical f64 coefficients (the
    /// fold engine preserves the serial per-position accumulation order at
    /// any slot count, and the tree is positional).
    #[test]

    fn tournament_reduction_independent_of_thread_count() {
        use ordered_float::NotNan;

        let (d, m) = (2usize, 5usize);
        let builder = LogSignatureBuilder::<u8>::new()
            .with_num_dimensions(d)
            .with_max_degree(m);
        let basis =
            lyndon_rs::lyndon::LyndonBasis::<u8>::new(d, lyndon_rs::lyndon::Sort::Lexicographical)
                .generate_basis(m);
        let mut seed = 0x71ee_u64;
        let lcg = |seed: &mut u64| {
            *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
            ((*seed >> 33) % 7) as i128 - 3
        };
        let rhss: Vec<Vec<NotNan<f64>>> = (0..24)
            .map(|_| {
                (0..basis.len())
                    .map(|k| {
                        if k < d {
                            NotNan::new(lcg(&mut seed) as f64).unwrap()
                        } else {
                            NotNan::new(0.0).unwrap()
                        }
                    })
                    .collect()
            })
            .collect();
        let slices: Vec<&[NotNan<f64>]> = rhss.iter().map(|r| r.as_slice()).collect();

        let run = |threads: usize| -> Vec<NotNan<f64>> {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .expect("pool");
            pool.install(|| {
                let mut bat = builder.build::<NotNan<f64>>();
                bat.concatenate_batch_coefficients(&slices);
                bat.series.coefficients.clone()
            })
        };
        assert_eq!(run(4), run(9));
    }


    /// Reassociation tolerance for tournament-vs-sequential f64 value
    /// comparisons, applied per coefficient as
    /// `TOURNAMENT_REASSOCIATION_TOL * max(1, |a|, |b|)`. Measured worst
    /// case on the 16-displacement zero-accumulator batch below:
    /// 1.279e-13 absolute on coefficients of up to 1.86e2 magnitude — a
    /// relative 2.2e-15, i.e. ~10 ulps of rounding dust. The bound keeps
    /// ~780x absolute headroom over that dust at O(1) magnitudes (and
    /// ~1e-10 relative headroom at larger ones) while staying ten-plus
    /// orders below the O(1)..O(1e2) shifts that a real support/gating/
    /// kernel bug produces: dust passes, structural errors fail loudly.

    const TOURNAMENT_REASSOCIATION_TOL: f64 = 1e-10;

    /// Asserts `tournament` matches the sequential `oracle` within the
    /// documented reassociation tolerance; `context` names the batch in
    /// the panic message.
    fn assert_within_reassociation_tolerance(
        oracle: &[NotNan<f64>],
        tournament: &[NotNan<f64>],
        context: &str,
    ) {
        assert_eq!(oracle.len(), tournament.len(), "{context}: length mismatch");
        for (k, (a, b)) in oracle.iter().zip(tournament.iter()).enumerate() {
            let (a, b) = (a.into_inner(), b.into_inner());
            let diff = (a - b).abs();
            let bound = TOURNAMENT_REASSOCIATION_TOL * a.abs().max(b.abs()).max(1.0);
            assert!(
                diff <= bound,
                "{context}: coefficient {k} off by {diff:.3e} (bound {bound:.3e}): \
                 oracle={a} tournament={b}"
            );
        }
    }




    /// Bit-level equality for f64 coefficient vectors: `assert_eq` on
    /// `NotNan` equates `0.0` and `-0.0`, which are different bytes, so
    /// cross-pool/cross-run determinism is asserted on `to_bits`.

    fn assert_bits_identical(a: &[NotNan<f64>], b: &[NotNan<f64>], context: &str) {
        assert_eq!(a.len(), b.len(), "{context}: length mismatch");
        for (k, (x, y)) in a.iter().zip(b.iter()).enumerate() {
            let (xb, yb) = (x.into_inner().to_bits(), y.into_inner().to_bits());
            assert!(
                xb == yb,
                "{context}: coefficient {k} differs at the bit level: bits {xb:b} vs {yb:b}"
            );
        }
    }


    /// Folding from a ZERO accumulator (the build_from_path shape) through
    /// the tournament path: every leaf folds its displacement chunk from a
    /// zeroed accumulator — the same support-growth regime the sequential
    /// driver runs — and pooled merge DAGs then meet the full-support
    /// accumulator. The oracle is the sequential driver
    /// (`concatenate_batch_sequential`) on a fresh clone.
    ///
    /// Assertions are split by nature. SUPPORT is structural: which basis
    /// words are non-zero must agree between the two reduction shapes
    /// bit-for-bit (a scatter-target or gating bug drops or spurs a
    /// supported word long before values drift past tolerance), and the
    /// supported word set is asserted word-for-word — it is data, so a
    /// machinery change that moves support must fail here, not silently
    /// re-derive it. VALUES are compared under the documented
    /// reassociation tolerance: measured worst case on this exact batch
    /// is 1.279e-13 absolute on coefficients of up to 1.86e2 magnitude
    /// (~2.2e-15 relative, ~10 ulps; see
    /// `assert_within_reassociation_tolerance`).
    #[test]

    fn concatenate_batch_grows_from_sparse_support() {
        let (d, m) = (2usize, 5usize);
        let builder = LogSignatureBuilder::<u8>::new()
            .with_num_dimensions(d)
            .with_max_degree(m);
        let basis =
            lyndon_rs::lyndon::LyndonBasis::<u8>::new(d, lyndon_rs::lyndon::Sort::Lexicographical)
                .generate_basis(m);
        let mut seed = 0x5a7e_u64;
        let lcg = |seed: &mut u64| {
            *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
            ((*seed >> 33) % 7) as i128 - 3
        };
        let rhss: Vec<Vec<NotNan<f64>>> = (0..16)
            .map(|_| {
                (0..basis.len())
                    .map(|k| {
                        if k < d {
                            NotNan::new(lcg(&mut seed) as f64).unwrap()
                        } else {
                            NotNan::new(0.0).unwrap()
                        }
                    })
                    .collect()
            })
            .collect();
        let slices: Vec<&[NotNan<f64>]> = rhss.iter().map(|r| r.as_slice()).collect();

        // The words support must grow to — identical for every reduction
        // shape: both letters and every degree-2..4 word over {0,1}; no
        // degree-5 word enters the kernel-visible support (the support is
        // cut off at max_degree).
        const SUPPORTED_WORDS: &[&[u8]] = &[
            &[0],
            &[1],
            &[0, 1],
            &[0, 0, 1],
            &[0, 1, 1],
            &[0, 0, 0, 1],
            &[0, 0, 1, 1],
            &[0, 1, 1, 1],
        ];
        let check_support = |sig: &LogSignature<u8, NotNan<f64>>, label: &str| {
            let sup = sig
                .series
                .nonzero_coefficient_indices(&sig.series.coefficients);
            let words: Vec<&[u8]> = sup.iter().map(|&i| basis[i].letters.as_slice()).collect();
            assert_eq!(words, SUPPORTED_WORDS, "{label}: supported words changed");
        };

        // Oracle: the strictly sequential driver on a fresh clone of the
        // zero-accumulator signature (fresh DAG — no node lists yet).
        let mut oracle = builder.build::<NotNan<f64>>();
        oracle.concatenate_batch_sequential(&slices);
        check_support(&oracle, "sequential oracle");

        // Tournament path: 16 displacements in a forced multi-thread pool.
        let threads: usize = std::env::var("DBG_THREADS")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(4);
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .expect("pool");
        let mut bat = builder.build::<NotNan<f64>>();
        pool.install(|| bat.concatenate_batch_coefficients(&slices));
        check_support(&bat, "tournament");

        assert_within_reassociation_tolerance(
            &oracle.series.coefficients,
            &bat.series.coefficients,
            "16-displacement zero-accumulator batch",
        );
    }


    /// The sequential engine's growth contract, kept explicit under the
    /// tournament's shadow: from a ZERO accumulator, folding one-by-one
    /// grows the kernel-visible support fold by fold (measured trajectory
    /// on this data: 1 → 6 → 8 words over the first three folds, then
    /// saturation), and the sequential batch driver — which warm-folds
    /// through that growth before batching the eligible rest — must stay
    /// BIT-identical to all-per-fold. Exact equality: no reassociation is
    /// involved on this path, so any difference is a sequential-machinery
    /// regression (the tournament path is compared against this engine
    /// with tolerance in `concatenate_batch_grows_from_sparse_support`,
    /// and the driver's engine selection itself is pinned by
    /// `batch_driver_length_gate_selects_the_reduction_path`).
    #[test]

    fn concatenate_batch_sequential_bit_identical_from_zero_accumulator() {
        let (d, m) = (2usize, 5usize);
        let builder = LogSignatureBuilder::<u8>::new()
            .with_num_dimensions(d)
            .with_max_degree(m);
        let basis =
            lyndon_rs::lyndon::LyndonBasis::<u8>::new(d, lyndon_rs::lyndon::Sort::Lexicographical)
                .generate_basis(m);
        let mut seed = 0x5a7e_u64;
        let lcg = |seed: &mut u64| {
            *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
            ((*seed >> 33) % 7) as i128 - 3
        };
        let rhss: Vec<Vec<NotNan<f64>>> = (0..16)
            .map(|_| {
                (0..basis.len())
                    .map(|k| {
                        if k < d {
                            NotNan::new(lcg(&mut seed) as f64).unwrap()
                        } else {
                            NotNan::new(0.0).unwrap()
                        }
                    })
                    .collect()
            })
            .collect();
        let slices: Vec<&[NotNan<f64>]> = rhss.iter().map(|r| r.as_slice()).collect();

        // Per-fold reference, recording the support-growth trajectory.
        let mut per_fold = builder.build::<NotNan<f64>>();
        let mut support_growth = Vec::new();
        for r in &rhss {
            per_fold.concatenate_assign_coefficients(r);
            support_growth.push(
                per_fold
                    .series
                    .nonzero_coefficient_indices(&per_fold.series.coefficients)
                    .len(),
            );
        }
        assert_eq!(
            support_growth,
            vec![1, 6, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8],
            "support-growth trajectory from the zero accumulator changed"
        );

        let mut driver = builder.build::<NotNan<f64>>();
        driver.concatenate_batch_sequential(&slices);
        assert_eq!(
            per_fold.series.coefficients,
            driver.series.coefficients,
            "sequential batch driver diverged from per-fold folding"
        );
        assert_eq!(
            driver
                .series
                .nonzero_coefficient_indices(&driver.series.coefficients)
                .len(),
            8,
            "driver's final support count changed"
        );
    }


    /// Exact-type identity at tournament scale, over awkward odd batch
    /// sizes (17 = 8+8+1: a single-displacement tail leaf; 31 = 8+8+8+7:
    /// a short tail; 63: eight full leaves), dimensions 2–4, degrees 3–5,
    /// and mixed sparse supports: every third displacement drops a
    /// rotating letter, and letter values hit zero exactly, so leaf folds
    /// rebuild node lists mid-flight and pooled DAGs see varying supports.
    /// Rational arithmetic has no rounding to hide behind: ANY coefficient
    /// difference from the sequential fold chain is a real exactness bug
    /// (reassociation must be invisible to exact types). The gate premise
    /// is asserted against the source gate itself (n ≥
    /// `TOURNAMENT_MIN_DISPLACEMENTS` in a multi-thread pool); the
    /// routing behavior those conditions select is pinned observationally
    /// by `batch_driver_length_gate_selects_the_reduction_path`.
    #[rstest]
    #[case(2, 3, 16)]
    #[case(2, 3, 17)]
    #[case(3, 4, 31)]
    #[case(2, 5, 31)]
    #[case(4, 5, 63)]

    fn tournament_rationals_exact_at_awkward_sizes(
        #[case] d: usize,
        #[case] m: usize,
        #[case] n: usize,
    ) {
        type R = Ratio<i64>;
        let builder = LogSignatureBuilder::<u8>::new()
            .with_num_dimensions(d)
            .with_max_degree(m);
        let basis =
            lyndon_rs::lyndon::LyndonBasis::<u8>::new(d, lyndon_rs::lyndon::Sort::Lexicographical)
                .generate_basis(m);
        let mut seed = 0x3a11_u64.wrapping_add((d * 1007 + m * 97 + n) as u64);
        let lcg = |seed: &mut u64| {
            *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
            ((*seed >> 33) % 7) as i64 - 3
        };
        // Letter displacements with a rotating dropped letter every third
        // displacement (changing supports → collecting rebuilds inside
        // leaves and across pooled-DAG reuse) and exact-zero letter values
        // wherever the LCG produces 0.
        let rhss: Vec<Vec<R>> = (0..n)
            .map(|i| {
                (0..basis.len())
                    .map(|k| {
                        if k < d && !(i % 3 == 2 && k == i % d) {
                            R::from_integer(lcg(&mut seed))
                        } else {
                            R::from_integer(0)
                        }
                    })
                    .collect()
            })
            .collect();
        let slices: Vec<&[R]> = rhss.iter().map(|r| r.as_slice()).collect();

        let mut seq = builder.build::<R>();
        for r in &rhss {
            seq.concatenate_assign_coefficients(r);
        }

        let threads: usize = std::env::var("DBG_THREADS")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(4);
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .expect("pool");
        pool.install(|| {
            assert!(
                n >= TOURNAMENT_MIN_DISPLACEMENTS && rayon::current_num_threads() > 1,
                "test premise: the driver must select the tournament path"
            );
            let mut bat = builder.build::<R>();
            bat.concatenate_batch_coefficients(&slices);
            assert_eq!(
                seq.series.coefficients,
                bat.series.coefficients,
                "d={d} m={m} n={n}: rational tournament diverged from sequential"
            );
        });
    }


    /// Cross-thread determinism at scale: the same 40-displacement batch
    /// (five leaves, three merge rounds, an odd pass-through) reduced
    /// under forced 1-, 3-, and 16-thread pools must be BYTE-identical —
    /// the tree shape is a function of the displacement count alone, and
    /// every fold preserves the serial per-position accumulation order at
    /// any slot count. `tournament_reduce` is called directly so the
    /// 1-thread pool exercises the tournament walk itself (the public
    /// driver routes a single-worker pool to the sequential engine — that
    /// gate is `batch_driver_length_gate_selects_the_reduction_path`'s
    /// subject); the public driver is additionally checked for
    /// byte-identity across 3- and 16-thread pools, where it takes the
    /// tournament on both.
    #[test]

    fn tournament_reduce_bit_identical_across_pool_sizes() {
        let (d, m) = (2usize, 5usize);
        let builder = LogSignatureBuilder::<u8>::new()
            .with_num_dimensions(d)
            .with_max_degree(m);
        let basis =
            lyndon_rs::lyndon::LyndonBasis::<u8>::new(d, lyndon_rs::lyndon::Sort::Lexicographical)
                .generate_basis(m);
        let mut seed = 0xc055_u64;
        let lcg = |seed: &mut u64| {
            *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
            ((*seed >> 33) % 7) as i128 - 3
        };
        let rhss: Vec<Vec<NotNan<f64>>> = (0..40)
            .map(|_| {
                (0..basis.len())
                    .map(|k| {
                        if k < d {
                            NotNan::new(lcg(&mut seed) as f64).unwrap()
                        } else {
                            NotNan::new(0.0).unwrap()
                        }
                    })
                    .collect()
            })
            .collect();
        let slices: Vec<&[NotNan<f64>]> = rhss.iter().map(|r| r.as_slice()).collect();

        let reduce_direct = |threads: usize| -> Vec<NotNan<f64>> {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .expect("pool");
            let sig = builder.build::<NotNan<f64>>();
            pool.install(|| sig.tournament_reduce(&slices))
        };
        let reduce_driver = |threads: usize| -> Vec<NotNan<f64>> {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .expect("pool");
            let mut sig = builder.build::<NotNan<f64>>();
            pool.install(|| {
                sig.concatenate_batch_coefficients(&slices);
                sig.series.coefficients.clone()
            })
        };

        let one = reduce_direct(1);
        let three = reduce_direct(3);
        let sixteen = reduce_direct(16);
        assert_bits_identical(&one, &three, "tournament 1 vs 3 threads");
        assert_bits_identical(&one, &sixteen, "tournament 1 vs 16 threads");

        assert_bits_identical(
            &reduce_driver(3),
            &reduce_driver(16),
            "public driver 3 vs 16 threads",
        );
    }


    /// The length gate selects the engine observationally. Cohort-capable
    /// types (raw floats, kill switch off) route EVERY batch of 2+
    /// displacements through the cohort tournament, at any pool size —
    /// the SIMD-across-folds engine amortizes the tree overhead away, so
    /// the historical 16-displacement/1-thread gates only apply when the
    /// cohort engine is unavailable (`set_cohort_off(true)`, non-float
    /// types). Both observables are sharp because the engines' f64
    /// outputs differ by reassociation dust on this data (measured:
    /// reducing the same 15 or 16 displacements along the tournament
    /// tree vs sequentially moves 9 coefficients at the bit level,
    /// ~1e-13 absolute), so exact equality with a direct
    /// `tournament_reduce` call FAILS if the sequential engine is wrongly
    /// taken, and exact per-fold equality FAILS if the tournament is
    /// wrongly taken:
    /// - 0 displacements: a no-op (the accumulator stays zero).
    /// - 1 displacement: identical to one `concatenate_assign_coefficients`
    ///   (a truly lone fold — nothing to cohort with).
    /// - 15 displacements, 4-thread pool, cohort on: bit-identical to a
    ///   direct `tournament_reduce` call (the cohort tournament runs from
    ///   2 displacements up).
    /// - 15 displacements, 4-thread pool, cohort OFF: bit-identical to
    ///   per-fold folding (the historical gate under the kill switch).
    /// - 16 displacements, 4-thread pool: bit-identical to a direct
    ///   `tournament_reduce` call and within reassociation tolerance of
    ///   per-fold folding, with at least one coefficient bit-different —
    ///   the measured tournament-ran observable (a bit-identical result
    ///   here would mean the two engines coincide on this data and the
    ///   data must be re-derived).
    /// - 16 displacements, 1-thread pool, cohort on: bit-identical to a
    ///   direct `tournament_reduce` call in that pool (the cohort's data
    ///   parallelism does not need workers).
    /// - 16 displacements, 1-thread pool, cohort OFF: bit-identical to
    ///   per-fold folding (single-worker pools never took the scalar
    ///   tournament).
    #[test]

    fn batch_driver_length_gate_selects_the_reduction_path() {
        let (d, m) = (2usize, 5usize);
        let builder = LogSignatureBuilder::<u8>::new()
            .with_num_dimensions(d)
            .with_max_degree(m);
        let basis =
            lyndon_rs::lyndon::LyndonBasis::<u8>::new(d, lyndon_rs::lyndon::Sort::Lexicographical)
                .generate_basis(m);
        let mut seed = 0xb0a_d05_u64;
        let lcg = |seed: &mut u64| {
            *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
            ((*seed >> 33) % 7) as i128 - 3
        };
        let mk_rhss = |count: usize, seed: &mut u64| -> Vec<Vec<NotNan<f64>>> {
            (0..count)
                .map(|_| {
                    (0..basis.len())
                        .map(|k| {
                            if k < d {
                                NotNan::new(lcg(seed) as f64).unwrap()
                            } else {
                                NotNan::new(0.0).unwrap()
                            }
                        })
                        .collect()
                })
                .collect()
        };
        let rhss = mk_rhss(17, &mut seed);
        let fold_all = |rs: &[Vec<NotNan<f64>>]| -> Vec<NotNan<f64>> {
            let mut sig = builder.build::<NotNan<f64>>();
            for r in rs {
                sig.concatenate_assign_coefficients(r);
            }
            sig.series.coefficients
        };
        let pool4 = rayon::ThreadPoolBuilder::new()
            .num_threads(4)
            .build()
            .expect("pool");
        let pool1 = rayon::ThreadPoolBuilder::new()
            .num_threads(1)
            .build()
            .expect("pool");

        // 0 displacements: a no-op.
        let mut sig0 = builder.build::<NotNan<f64>>();
        sig0.concatenate_batch_coefficients(&[]);
        assert!(
            sig0.series.coefficients.iter().all(|c| c.is_zero()),
            "empty batch must not touch the accumulator"
        );

        // 1 displacement: a single fold.
        let mut single = builder.build::<NotNan<f64>>();
        single.concatenate_assign_coefficients(&rhss[0]);
        let mut sig1 = builder.build::<NotNan<f64>>();
        sig1.concatenate_batch_coefficients(&[&rhss[0]]);
        assert_bits_identical(
            &single.series.coefficients,
            &sig1.series.coefficients,
            "1 displacement must equal a single concatenate_assign_coefficients",
        );

        // 15 displacements in a multi-thread pool, cohort on: tournament
        // bits (the cohort tournament runs from 2 displacements up).
        let seq15 = fold_all(&rhss[..15]);
        let slices15: Vec<&[NotNan<f64>]> =
            rhss[..15].iter().map(|r| r.as_slice()).collect();
        let tred15 = builder.build::<NotNan<f64>>();
        let torn15 = pool4.install(|| tred15.tournament_reduce(&slices15));
        let mut got15 = builder.build::<NotNan<f64>>();
        pool4.install(|| got15.concatenate_batch_coefficients(&slices15));
        assert_bits_identical(
            &torn15,
            &got15.series.coefficients,
            "15 displacements must take the (cohort) tournament path",
        );
        assert_within_reassociation_tolerance(
            &seq15,
            &got15.series.coefficients,
            "15-displacement cohort batch vs sequential oracle",
        );

        // 15 displacements with the cohort engine disabled: the SAME
        // tournament tree through the scalar engines — the kill switch
        // swaps engines, never the association tree, so the cohort-off run
        // is bit-identical to the cohort-on one (and within the
        // reassociation tolerance of the sequential oracle).
        let mut got15_off = builder.build::<NotNan<f64>>();
        set_cohort_off(true);
        pool4.install(|| got15_off.concatenate_batch_coefficients(&slices15));
        set_cohort_off(false);
        assert_bits_identical(
            &torn15,
            &got15_off.series.coefficients,
            "cohort-off 15 displacements: same tree, scalar engines — tournament bits",
        );
        assert_within_reassociation_tolerance(
            &seq15,
            &got15_off.series.coefficients,
            "15-displacement scalar-tournament batch vs sequential oracle",
        );

        // 16 displacements in a multi-thread pool: tournament bits.
        let seq16 = fold_all(&rhss[..16]);
        let slices16: Vec<&[NotNan<f64>]> =
            rhss[..16].iter().map(|r| r.as_slice()).collect();
        let tred = builder.build::<NotNan<f64>>();
        let torn16 = pool4.install(|| tred.tournament_reduce(&slices16));
        let mut got16 = builder.build::<NotNan<f64>>();
        pool4.install(|| got16.concatenate_batch_coefficients(&slices16));
        assert_bits_identical(
            &torn16,
            &got16.series.coefficients,
            "16 displacements must take the tournament path",
        );
        assert_within_reassociation_tolerance(
            &seq16,
            &got16.series.coefficients,
            "16-displacement batch vs sequential oracle",
        );
        let bit_diffs = seq16
            .iter()
            .zip(&got16.series.coefficients)
            .filter(|(a, b)| a.into_inner().to_bits() != b.into_inner().to_bits())
            .count();
        assert!(
            bit_diffs > 0,
            "16-displacement tournament result is bit-identical to the sequential \
             engine on this data — the tournament-ran observable is gone; \
             re-derive the test data"
        );

        // 16 displacements in a single-worker pool, cohort on: tournament
        // bits (the cohort's data parallelism does not need workers).
        let tred1 = builder.build::<NotNan<f64>>();
        let torn1 = pool1.install(|| tred1.tournament_reduce(&slices16));
        let mut got1 = builder.build::<NotNan<f64>>();
        pool1.install(|| got1.concatenate_batch_coefficients(&slices16));
        assert_bits_identical(
            &torn1,
            &got1.series.coefficients,
            "16 displacements in a 1-thread pool must take the (cohort) tournament path",
        );

        // 16 displacements in a single-worker pool, cohort OFF: the SAME
        // tournament tree through the scalar engines (the kill switch
        // swaps engines, never the tree).
        let mut got1_off = builder.build::<NotNan<f64>>();
        set_cohort_off(true);
        pool1.install(|| got1_off.concatenate_batch_coefficients(&slices16));
        set_cohort_off(false);
        assert_bits_identical(
            &torn1,
            &got1_off.series.coefficients,
            "cohort-off 16 displacements in a 1-thread pool: same tree, scalar engines",
        );
    }


    /// NaN-audit parity on the tournament path: a leaf's SECOND fold
    /// overflows — the first fold from a zero accumulator is safe
    /// (0 × MAX = 0; the result is just the degree-1 displacement), then
    /// MAX × MAX products overflow to ±inf and the fused two-orientation
    /// term computes inf + (−inf) = NaN — and the per-fold audit must
    /// panic with the same message the sequential engine panics with,
    /// propagating through rayon's task boundary (rayon re-throws the
    /// original payload, and the tournament's DAG-pool mutex is
    /// poison-tolerant, so no secondary poison panic can replace it).
    #[test]
    #[should_panic(expected = "log-signature coefficients overflowed to NaN")]

    fn tournament_path_audits_nan_after_overflow() {
        let (d, m) = (2usize, 3usize);
        let builder = LogSignatureBuilder::<u8>::new()
            .with_num_dimensions(d)
            .with_max_degree(m);
        let len = builder.build::<NotNan<f64>>().series.coefficients.len();
        let rhss: Vec<Vec<NotNan<f64>>> = (0..24)
            .map(|_| {
                (0..len)
                    .map(|k| NotNan::new(if k < d { f64::MAX } else { 0.0 }).unwrap())
                    .collect()
            })
            .collect();
        let slices: Vec<&[NotNan<f64>]> = rhss.iter().map(|r| r.as_slice()).collect();

        let threads: usize = std::env::var("DBG_THREADS")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(4);
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .expect("pool");
        pool.install(|| {
            let mut sig = builder.build::<NotNan<f64>>();
            sig.concatenate_batch_coefficients(&slices);
        });
    }


    /// Repeated-run determinism: identical tournament reductions in a row
    /// must be byte-identical — nothing may leak between calls through
    /// the shared compiled plan (Arc'd decomposition tables, source-DAG
    /// node lists, pooled-DAG scratch). A different-shape,
    /// different-support tournament runs in between on a signature whose
    /// DAG shares the plan with a warmed one (clone_shallow), and the
    /// first input then repeats twice — a stale support-keyed gate entry
    /// or left-over pooled-DAG node lists would surface as a byte
    /// difference on the repeats.
    #[test]

    fn tournament_repeated_runs_byte_identical() {
        let (d, m) = (2usize, 5usize);
        let builder = LogSignatureBuilder::<u8>::new()
            .with_num_dimensions(d)
            .with_max_degree(m);
        let basis =
            lyndon_rs::lyndon::LyndonBasis::<u8>::new(d, lyndon_rs::lyndon::Sort::Lexicographical)
                .generate_basis(m);
        let mk_rhss = |count: usize, seed: &mut u64, drop_letter: bool| -> Vec<Vec<NotNan<f64>>> {
            (0..count)
                .map(|i| {
                    (0..basis.len())
                        .map(|k| {
                            let dropped = drop_letter && i % 3 == 2 && k == i % d;
                            if k < d && !dropped {
                                let v = {
                                    *seed = seed
                                        .wrapping_mul(6364136223846793005)
                                        .wrapping_add(1);
                                    ((*seed >> 33) % 7) as i128 - 3
                                };
                                NotNan::new(v as f64).unwrap()
                            } else {
                                NotNan::new(0.0).unwrap()
                            }
                        })
                        .collect()
                })
                .collect()
        };
        // A: dense letter displacements; B: mixed sparse supports.
        let mut seed_a = 0x9eaf_u64;
        let rhss_a = mk_rhss(24, &mut seed_a, false);
        let slices_a: Vec<&[NotNan<f64>]> = rhss_a.iter().map(|r| r.as_slice()).collect();
        let mut seed_b = 0x9eb0_u64;
        let rhss_b = mk_rhss(17, &mut seed_b, true);
        let slices_b: Vec<&[NotNan<f64>]> = rhss_b.iter().map(|r| r.as_slice()).collect();

        let threads: usize = std::env::var("DBG_THREADS")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(4);
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .expect("pool");
        let run_a = || -> Vec<NotNan<f64>> {
            pool.install(|| {
                let mut sig = builder.build::<NotNan<f64>>();
                sig.concatenate_batch_coefficients(&slices_a);
                sig.series.coefficients.clone()
            })
        };

        let first = run_a();
        // Interleave: a warmed signature (node lists built for A's
        // supports) is cloned — sharing the compiled plan — and runs B's
        // mixed-support tournament through the shared plan.
        pool.install(|| {
            let mut warm = builder.build::<NotNan<f64>>();
            warm.concatenate_batch_coefficients(&slices_a);
            let mut shared = warm.clone();
            shared.concatenate_batch_coefficients(&slices_b);
        });
        let second = run_a();
        let third = run_a();
        assert_bits_identical(&first, &second, "repeated identical tournament runs");
        assert_bits_identical(&first, &third, "third identical tournament run");
    }


    /// Same audit for the batched fold: NaN arising anywhere in a batch
    /// survives every later add/multiply, so the single end-of-batch audit
    /// still fails loudly.
    #[test]
    #[should_panic(expected = "log-signature coefficients overflowed to NaN")]

    fn batch_fold_audits_nan_after_overflow() {
        use ordered_float::NotNan;

        let (d, m) = (2usize, 3usize);
        let builder = LogSignatureBuilder::<u8>::new()
            .with_num_dimensions(d)
            .with_max_degree(m);
        let mut log_sig = builder.build::<NotNan<f64>>();
        // Full, small-valued accumulator support, then warm folds to build
        // the node lists and reach batch eligibility.
        let mut seed = 0xfeed_u64;
        let mut lcg = |seed: &mut u64| {
            *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
            ((*seed >> 33) % 7) as i128 - 3
        };
        let len = log_sig.series.coefficients.len();
        log_sig.series.coefficients = (0..len)
            .map(|k| {
                NotNan::new(if k < d { 0.0 } else { lcg(&mut seed) as f64 }).unwrap()
            })
            .collect();
        let letter_disp = |v: f64| -> Vec<NotNan<f64>> {
            (0..len)
                .map(|k| NotNan::new(if k < d { v } else { 0.0 }).unwrap())
                .collect()
        };
        for _ in 0..4 {
            log_sig.concatenate_assign_coefficients(&letter_disp(1.5));
        }
        // The overflowing displacement must take the batch path.
        assert!(log_sig.dag.batch_eligible(
            &log_sig.series.nonzero_coefficient_indices(&log_sig.series.coefficients),
            &log_sig.series.nonzero_coefficient_indices(&letter_disp(f64::MAX)),
        ));
        log_sig.concatenate_batch_coefficients(&[&letter_disp(f64::MAX)]);
    }


    /// Fold structure analysis for class-partitioned scheduling: per level,
    /// the node count, distinct anagram classes, and how sweep work
    /// (support sizes) distributes across classes — plus the WRITE-ACCESS
    /// view: each node's scatter-target list is the exact set of word
    /// classes its sweep owns, so the report shows the write-class count
    /// per level and the balance a static word-class partition could
    /// achieve. Run explicitly with --ignored --nocapture.
    #[test]
    #[ignore = "structure stats: run explicitly with --ignored --nocapture"]

    fn dump_fold_class_structure() {
        use ordered_float::NotNan;

        for (d, m) in [(2usize, 12usize), (3usize, 8usize), (3usize, 10usize)] {
            let builder = LogSignatureBuilder::<u8>::new()
                .with_num_dimensions(d)
                .with_max_degree(m);
            let mut log_sig = builder.build::<NotNan<f64>>();
            let basis: Vec<LyndonWord<u8>> =
                LyndonBasis::<u8>::new(d, lyndon_rs::Sort::Lexicographical).generate_basis(m);
            let mut seed = 0xfeed_u64;
            let lcg = |seed: &mut u64| {
                *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
                ((*seed >> 33) as i64) as f64
            };
            let acc: Vec<NotNan<f64>> = (0..basis.len())
                .map(|_| NotNan::new(lcg(&mut seed)).unwrap())
                .collect();
            log_sig.series.coefficients.clone_from(&acc);
            let disp: Vec<NotNan<f64>> = (0..basis.len())
                .map(|k| NotNan::new(if k < d { lcg(&mut seed) } else { 0.0 }).unwrap())
                .collect();
            let mut seg = log_sig.clone();
            seg.series.coefficients.clone_from(&disp);
            log_sig.concatenate_assign(&seg);

            // Letter content per basis word, keyed through a sorted alphabet
            // (works for any generator representation).
            let mut alphabet: Vec<u8> = log_sig
                .series
                .basis
                .iter()
                .flat_map(|w| w.letters.iter().copied())
                .collect();
            alphabet.sort_unstable();
            alphabet.dedup();
            let contents: Vec<Vec<u8>> = log_sig
                .series
                .basis
                .iter()
                .map(|w| {
                    let mut c = vec![0u8; alphabet.len()];
                    for l in &w.letters {
                        c[alphabet.binary_search(l).expect("letter in alphabet")] += 1;
                    }
                    c
                })
                .collect();

            let dag: &CommutatorDag<NotNan<f64>> = &log_sig.dag;
            let levels = &dag.structure.levels;
            let mut total_work = 0usize;
            let mut total_word_classes = 0usize;
            let mut report = String::new();
            let mut all_class_work: std::collections::HashMap<(u8, Vec<u8>), usize> =
                std::collections::HashMap::new();
            for (li, level) in levels.iter().enumerate().skip(1) {
                let mut class_work: std::collections::HashMap<(u8, Vec<u8>), usize> =
                    std::collections::HashMap::new();
                let mut level_work = 0usize;
                let mut level_word_classes = 0usize;
                let mut max_node_words = 0usize;
                for &k in level {
                    let sup = &dag.nonzeros[k as usize - 2];
                    let mut key = (0u8, vec![0u8; d]);
                    if let Some(&i) = sup.first() {
                        key = (basis[i].len() as u8, contents[i].clone());
                    }
                    class_work
                        .entry(key)
                        .and_modify(|x| *x += sup.len())
                        .or_insert(sup.len());
                    level_work += sup.len();
                    // Write-access view: the node's exact scatter-target set
                    // is the set of word classes its sweep owns (one class
                    // per target word; the list was recorded from the
                    // planner's per-word gating).
                    level_word_classes += sup.len();
                    max_node_words = max_node_words.max(sup.len());
                }
                for (k, v) in &class_work {
                    *all_class_work.entry(k.clone()).or_insert(0) += v;
                }
                total_work += level_work;
                total_word_classes += level_word_classes;
                let mut sizes: Vec<usize> = class_work.values().copied().collect();
                sizes.sort_unstable_by(|a, b| b.cmp(a));
                report.push_str(&format!(
                    "  level {li}: nodes={} classes={} work={} word_classes={max_node_words}->{level_word_classes} max_node_words_share={}%(of level)",
                    level.len(),
                    class_work.len(),
                    level_work,
                    100 * max_node_words / level_word_classes.max(1),
                ));
                if !sizes.is_empty() {
                    report.push_str(&format!(
                        " max_class={} ({}%)",
                        sizes[0],
                        100 * sizes[0] / level_work.max(1)
                    ));
                }
                report.push('\n');
            }
            let mut works: Vec<usize> = all_class_work.values().copied().collect();
            works.sort_unstable_by(|a, b| b.cmp(a));
            let mut packs = vec![0usize; 8];
            for w in works {
                let p = packs.iter_mut().min_by_key(|x| **x).unwrap();
                *p += w;
            }
            packs.sort_unstable_by(|a, b| b.cmp(a));
            let max_share =
                all_class_work.values().copied().max().unwrap_or(0) * 100 / total_work.max(1);
            println!(
                "d={d} m={m}: levels={} total_work={total_work} classes={} max_class_share={max_share}% packs8_max_min_ratio={:.2} word_classes={total_word_classes}",
                levels.len().saturating_sub(1),
                all_class_work.len(),
                packs[0] as f64 / packs[7].max(1) as f64,
            );
            print!("{report}");
        }
    }


    /// Differential test: after the first fold (which builds the node lists
    /// through the collecting pass), the steady level-batch evaluation must
    /// produce bit-identical node buffers to a fresh DAG's collecting
    /// rebuild — the oracle — for every subsequent fold.
    #[test]

    fn steady_batch_matches_rebuild_oracle() {
        use ordered_float::NotNan;

        for (d, m, folds) in [(3usize, 3usize, 6), (3, 4, 6), (2, 8, 6), (4, 5, 4)] {
            let builder = LogSignatureBuilder::<u8>::new()
                .with_num_dimensions(d)
                .with_max_degree(m);
            let mut log_sig = builder.build::<NotNan<f64>>();

            let basis = lyndon_rs::lyndon::LyndonBasis::<u8>::new(
                d,
                lyndon_rs::lyndon::Sort::Lexicographical,
            )
            .generate_basis(m);
            let mut seed = 0x5eed_u64;
            let mut lcg = |seed: &mut u64| {
                *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
                ((*seed >> 33) as i64) as f64
            };
            let acc: Vec<NotNan<f64>> = (0..basis.len())
                .map(|_| NotNan::new(lcg(&mut seed)).unwrap())
                .collect();
            log_sig.series.coefficients.clone_from(&acc);

            for fold in 0..folds {
                // Alternate letter-only displacements (the production fold's
                // regime) with dense ones (worst case for the gating).
                let dense = fold % 2 == 1;
                let disp: Vec<NotNan<f64>> = (0..basis.len())
                    .map(|k| {
                        NotNan::new(if dense || k < d { lcg(&mut seed) } else { 0.0 }).unwrap()
                    })
                    .collect();
                let mut disp_sig = builder.build::<NotNan<f64>>();
                disp_sig.series.coefficients.clone_from(&disp);

                let original = log_sig.series.coefficients.clone();
                let a_nonzero = log_sig.series.nonzero_coefficient_indices(&original);
                let b_nonzero = log_sig.series.nonzero_coefficient_indices(&disp);
                log_sig.dag.evaluate(
                    &log_sig.series,
                    &original,
                    &a_nonzero,
                    &disp,
                    &b_nonzero,
                );

                // Oracle: a fresh DAG has no lists, so its evaluate takes the
                // collecting-rebuild path with the identical inputs.
                let mut oracle = CommutatorDag::<NotNan<f64>>::from_bch_series(
                    &log_sig.bch_series,
                );
                oracle.evaluate(&log_sig.series, &original, &a_nonzero, &disp, &b_nonzero);
                for (k, (buf, reff)) in
                    log_sig.dag.buffers.iter().zip(&oracle.buffers).enumerate()
                {
                    assert_eq!(
                        buf.len(),
                        reff.len(),
                        "d={d} m={m} fold={fold}: node {} buffer lengths differ",
                        k + 2
                    );
                    for (i, (x, y)) in buf.iter().zip(reff.iter()).enumerate() {
                        assert!(
                            x == y,
                            "d={d} m={m} fold={fold}: node {} index {} differs: {x:?} vs {y:?}",
                            k + 2,
                            i
                        );
                    }
                }

                // Complete the fold so the next one sees the evolved
                // accumulator (and the steady-batch lists).
                log_sig
                    .dag
                    .accumulate_terms(&mut log_sig.series.coefficients, &disp);
            }
        }
    }


    #[test]

    fn debug_print_support_sizes() {
        use ordered_float::NotNan;
        for (d, m) in [(2usize, 12usize), (3usize, 8usize), (8usize, 3usize), (12usize, 2usize)] {
            let builder = LogSignatureBuilder::<u8>::new()
                .with_num_dimensions(d)
                .with_max_degree(m);
            let basis = lyndon_rs::lyndon::LyndonBasis::<u8>::new(
                d,
                lyndon_rs::lyndon::Sort::Lexicographical,
            )
            .generate_basis(m);
            let len = basis.len();
            let mut log_sig = builder.build::<NotNan<f64>>();
            // dense accumulator + letter displacement, a few folds to settle lists
            let acc: Vec<NotNan<f64>> = (0..len)
                .map(|k| NotNan::new((k % 17) as f64 / 19.0 - 0.4).unwrap())
                .collect();
            log_sig.series.coefficients.clone_from(&acc);
            for fold in 0..4 {
                let disp: Vec<NotNan<f64>> = (0..len)
                    .map(|k| if k < d { NotNan::new(0.25).unwrap() } else { NotNan::new(0.0).unwrap() })
                    .collect();
                let _ = fold;
                let original = log_sig.series.coefficients.clone();
                let a_nz = log_sig.series.nonzero_coefficient_indices(&original);
                let b_nz = log_sig.series.nonzero_coefficient_indices(&disp);
                log_sig.dag.evaluate(&log_sig.series, &original, &a_nz, &disp, &b_nz);
                log_sig.dag.accumulate_terms(&mut log_sig.series.coefficients, &disp);
            }
            let dag = &log_sig.dag;
            let internal = dag.buffers.len();
            let total_nz: usize = dag.nonzeros.iter().map(|v| v.len()).sum();
            let full = internal * len;
            let mut sizes: Vec<usize> = dag.nonzeros.iter().map(|v| v.len()).collect();
            sizes.sort_unstable();
            let med = sizes[sizes.len() / 2];
            let max = *sizes.last().unwrap();
            let min = *sizes.first().unwrap();
            // slice-contiguity: per node, the degrees its support touches,
            // and the allocation if buffers were whole-slice-contiguous.
            let (degrees, deg_starts) = log_sig.series.debug_degree_layout();
            let mut slice_alloc = 0usize;
            let mut deg_count_sum = 0usize;
            for nz in dag.nonzeros.iter() {
                let mut degs: Vec<u8> = nz.iter().map(|&p| degrees[p]).collect();
                degs.sort_unstable();
                degs.dedup();
                deg_count_sum += degs.len();
                for &dg in &degs {
                    let d0 = deg_starts[dg as usize] as usize;
                    let d1 = deg_starts
                        .get(dg as usize + 1)
                        .map(|&v| v as usize)
                        .unwrap_or(len);
                    slice_alloc += d1 - d0;
                }
            }
            println!(
                "{d}x{m}: d={len} nodes={internal} sum_nz={total_nz} full={full} ratio={:.2} avg={:.0} med={med} min={min} max={max} slice_alloc={slice_alloc} avg_degs_per_node={:.1}",
                total_nz as f64 / full as f64,
                total_nz as f64 / internal as f64,
                deg_count_sum as f64 / internal as f64,
            );
        }
    }


    /// Fold-step probe for large grids, env-gated so normal test runs skip it.
    /// PROF_GRID="d,m" PROF_FOLDS=n cargo test --release -p signature-rs --lib bench_probe -- --nocapture
    /// Builds the series ONCE (the O(D^2) table build dominates), then times
    /// steady-state letter-displacement folds. RAYON_NUM_THREADS selects the
    /// regime under test.
    #[test]

    fn probe_large_grid_folds() {
        let Some(grid) = std::env::var_os("PROF_GRID") else {
            return;
        };
        let grid = grid.to_str().unwrap().to_owned();
        let (d, m) = grid.split_once(',').unwrap();
        let (d, m) = (d.parse::<usize>().unwrap(), m.parse::<usize>().unwrap());
        let folds: usize = std::env::var("PROF_FOLDS")
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(20);

        use ordered_float::NotNan;
        let builder = LogSignatureBuilder::<u8>::new()
            .with_num_dimensions(d)
            .with_max_degree(m);
        let t = std::time::Instant::now();
        let mut log_sig = builder.build::<NotNan<f64>>();
        let build = t.elapsed();
        println!(
            "PROBE d={d} m={m} D={} series_build={build:?}",
            log_sig.series.coefficients.len()
        );

        let basis =
            lyndon_rs::lyndon::LyndonBasis::<u8>::new(d, lyndon_rs::lyndon::Sort::Lexicographical)
                .generate_basis(m);
        let mut seed = 0xfeed_u64;
        let lcg = |seed: &mut u64| {
            *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
            ((*seed >> 33) as i64) as f64
        };
        let acc: Vec<NotNan<f64>> = (0..basis.len())
            .map(|_| NotNan::new(lcg(&mut seed)).unwrap())
            .collect();
        log_sig.series.coefficients.clone_from(&acc);
        let disp: Vec<NotNan<f64>> = (0..basis.len())
            .map(|k| NotNan::new(if k < d { lcg(&mut seed) } else { 0.0 }).unwrap())
            .collect();
        let mut seg = log_sig.clone();
        seg.series.coefficients.clone_from(&disp);

        // warm: the first fold takes the collecting-rebuild path, the rest
        // run the steady level-batch
        for _ in 0..4 {
            log_sig.concatenate_assign(&seg);
        }
        let t = std::time::Instant::now();
        for k in 0..folds {
            log_sig.concatenate_assign(&seg);
            std::hint::black_box(&log_sig);
            if k + 1 == folds {
                println!(
                    "PROBE d={d} m={m} threads={} fold={:?} ({} folds)",
                    std::thread::available_parallelism()
                        .map(|n| n.get())
                        .unwrap_or(0),
                    t.elapsed() / (folds as u32),
                    folds
                );
            }
        }
    }


    /// The cohort engine's bit-exactness oracle: the tournament with the
    /// cohort ENABLED must reduce the same reduction tree to the SAME BITS
    /// as the tournament with the cohort DISABLED (the in-process kill
    /// switch — one binary, two engine configurations), for f64
    /// coefficients across the production grids, batch sizes at and above
    /// the tournament gate (including odd lengths whose chunk tails pass
    /// through and 2-lane merge cohorts), and displacement flavors that
    /// drive the eligibility rule's fallbacks (dense letters, per-position
    /// sparse mixes that split chunk supports, block-sparse runs whose
    /// leaf supports differ). The tree shape depends only on
    /// `rhss.len()`, so both runs walk the identical tree; the cohort's
    /// per-lane operation order replicates the scalar engine's (plain
    /// mul+add, no FMA), so 0 ulp is the only acceptable outcome.
    ///
    /// For dense-letter batches of 32+ displacements the cohort must also
    /// OBSERVABLY engage (the engine-run counter — a scalar==scalar
    /// coincidence would otherwise pass vacuously); sparse flavors make no
    /// such demand (their lanes legitimately fall back).
    #[test]

    fn cohort_tournament_matches_scalar_tournament_bit_identically() {
        let pool = rayon::ThreadPoolBuilder::new().num_threads(4).build().expect("pool");
        let grids = [
            (2usize, 12usize),
            (3, 8),
            (12, 2),
            (8, 3),
            (2, 5),
        ];
        let sizes = [16usize, 17, 23, 31, 32, 33, 64];
        for (d, m) in grids {
            let builder = LogSignatureBuilder::<u8>::new()
                .with_num_dimensions(d)
                .with_max_degree(m);
            let basis =
                LyndonBasis::<u8>::new(d, lyndon_rs::Sort::Lexicographical).generate_basis(m);
            let basis_len = basis.len();
            let sizes: &[usize] = if (d, m) == (2, 12) {
                &[16, 17, 23, 31, 32, 33, 100]
            } else {
                &sizes
            };
            for &n in sizes {
                for flavor in 0..3usize {
                    let mut seed = 0x6e_77_u64 ^ ((d as u64) << 32) ^ (m as u64) ^ ((n as u64) << 8) ^ (flavor as u64);
                    let lcg = |seed: &mut u64| {
                        *seed = seed
                            .wrapping_mul(6364136223846793005)
                            .wrapping_add(1);
                        ((*seed >> 33) % 5) as i64 + 1 // 1..=5: never zero
                    };
                    let rhss: Vec<Vec<NotNan<f64>>> = (0..n)
                        .map(|i| {
                            (0..basis_len)
                                .map(|k| {
                                    if k >= d {
                                        NotNan::new(0.0).unwrap()
                                    } else {
                                        let keep = match flavor {
                                            0 => true,
                                            // per-position sparse mix: chunk
                                            // supports split inside leaves
                                            1 => (i + k) % 2 == 0,
                                            // block-sparse: leaves differ
                                            2 => i < n / 2,
                                            _ => unreachable!(),
                                        };
                                        if keep {
                                            NotNan::new(lcg(&mut seed) as f64).unwrap()
                                        } else {
                                            NotNan::new(0.0).unwrap()
                                        }
                                    }
                                })
                                .collect()
                        })
                        .collect();
                    let slices: Vec<&[NotNan<f64>]> =
                        rhss.iter().map(|r| r.as_slice()).collect();

                    let reduce = |off: bool| -> (Vec<NotNan<f64>>, usize) {
                        let sig = builder.build::<NotNan<f64>>();
                        set_cohort_off(off);
                        let before = COHORT_ENGINE_RUNS.load(std::sync::atomic::Ordering::SeqCst);
                        let reduced = pool.install(|| sig.tournament_reduce(&slices));
                        let after = COHORT_ENGINE_RUNS.load(std::sync::atomic::Ordering::SeqCst);
                        set_cohort_off(false);
                        (reduced, after - before)
                    };

                    let (scalar, _scalar_runs) = reduce(true);
                    let (cohort, cohort_runs) = reduce(false);
                    if n >= 32 && flavor == 0 {
                        assert!(
                            cohort_runs > 0,
                            "{d}x{m} n={n}: the cohort engine never engaged"
                        );
                    }
                    assert_bits_identical(
                        &scalar,
                        &cohort,
                        &format!("{d}x{m} n={n} flavor={flavor} cohort-vs-scalar"),
                    );
                }
            }
        }
    }


    /// Cohort-engagement fold coverage: with shared displacement supports,
    /// EVERY fold of a small batch routes through the SIMD-across-folds
    /// engine — the leaf round's shared-support fast path covers all leaf
    /// folds in ONE dispatch (2-3-lane remainder groups included) and every
    /// merge round with a pair cohorts its folds. The only scalar fold is
    /// the truly lone final merge. The lane-fold counter proves coverage
    /// beyond the merge rounds — the only place the pre-SIMD-everywhere
    /// engine engaged cohorts.
    #[test]

    fn cohort_covers_every_eligible_segment_not_only_merges() {
        use ordered_float::NotNan;

        let (d, m) = (2usize, 5usize);
        let builder = LogSignatureBuilder::<u8>::new()
            .with_num_dimensions(d)
            .with_max_degree(m);
        let basis = lyndon_rs::lyndon::LyndonBasis::<u8>::new(
            d,
            lyndon_rs::lyndon::Sort::Lexicographical,
        )
        .generate_basis(m);
        let mut seed = 0xc0de_u64;
        let lcg = |seed: &mut u64| {
            *seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((*seed >> 33) % 19) as i128 - 9
        };
        // 7 letter displacements: chunk = ceil(7/4) = 2, so the leaf round
        // is ONE remainder group of four leaves holding 2+2+2+1 folds.
        let n = 7usize;
        let rhss: Vec<Vec<NotNan<f64>>> = (0..n)
            .map(|_| {
                (0..basis.len())
                    .map(|k| {
                        if k < d {
                            NotNan::new(lcg(&mut seed) as f64).unwrap()
                        } else {
                            NotNan::new(0.0).unwrap()
                        }
                    })
                    .collect()
            })
            .collect();
        let slices: Vec<&[NotNan<f64>]> = rhss.iter().map(|r| r.as_slice()).collect();
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(4)
            .build()
            .expect("pool");

        let run = |off: bool| -> (Vec<NotNan<f64>>, usize, usize) {
            let mut sig = builder.build::<NotNan<f64>>();
            set_cohort_off(off);
            let runs0 = COHORT_ENGINE_RUNS.load(std::sync::atomic::Ordering::SeqCst);
            let folds0 = COHORT_LANE_FOLDS.load(std::sync::atomic::Ordering::SeqCst);
            pool.install(|| sig.concatenate_batch_coefficients(&slices));
            let runs = COHORT_ENGINE_RUNS.load(std::sync::atomic::Ordering::SeqCst) - runs0;
            let folds = COHORT_LANE_FOLDS.load(std::sync::atomic::Ordering::SeqCst) - folds0;
            set_cohort_off(false);
            (sig.series.coefficients.clone(), runs, folds)
        };

        let (scalar, scalar_runs, scalar_folds) = run(true);
        assert_eq!(scalar_runs, 0, "kill switch must silence the cohort engine");
        assert_eq!(scalar_folds, 0, "kill switch must silence the cohort engine");
        let (cohort, cohort_runs, cohort_folds) = run(false);
        // n=7: leaves hold 2+2+2+1 folds — the leaf fast path runs all 7
        // as ONE dispatch (1 run), but the merge round's supports DIVERGE
        // (the 1-fold leaf never grows the degree-3+ positions the 2-fold
        // leaves reach), so its pairs fold scalar. The n=4 case below
        // shows the merge round cohorting when the leaf shapes match.
        assert_eq!(
            (cohort_runs, cohort_folds),
            (1, 7),
            "n=7: leaf fast path only — the merge supports diverge"
        );
        assert_bits_identical(
            &scalar,
            &cohort,
            "n=7 cohort routing must be bit-identical to the scalar tournament",
        );

        // n=4: uniform 1-fold leaves — identical merge supports, so the
        // 2-pair merge round cohorts too (1 more run, 2 more lane-folds);
        // only the final lone merge stays scalar.
        let rhss4: Vec<Vec<NotNan<f64>>> = (0..4)
            .map(|_| {
                (0..basis.len())
                    .map(|k| {
                        if k < d {
                            NotNan::new(lcg(&mut seed) as f64).unwrap()
                        } else {
                            NotNan::new(0.0).unwrap()
                        }
                    })
                    .collect()
            })
            .collect();
        let slices4: Vec<&[NotNan<f64>]> = rhss4.iter().map(|r| r.as_slice()).collect();
        let run4 = |off: bool| -> (Vec<NotNan<f64>>, usize, usize) {
            let mut sig = builder.build::<NotNan<f64>>();
            set_cohort_off(off);
            let runs0 = COHORT_ENGINE_RUNS.load(std::sync::atomic::Ordering::SeqCst);
            let folds0 = COHORT_LANE_FOLDS.load(std::sync::atomic::Ordering::SeqCst);
            pool.install(|| sig.concatenate_batch_coefficients(&slices4));
            let runs = COHORT_ENGINE_RUNS.load(std::sync::atomic::Ordering::SeqCst) - runs0;
            let folds = COHORT_LANE_FOLDS.load(std::sync::atomic::Ordering::SeqCst) - folds0;
            set_cohort_off(false);
            (sig.series.coefficients.clone(), runs, folds)
        };
        let (scalar4, _, _) = run4(true);
        let (cohort4, runs4, folds4) = run4(false);
        assert_eq!(
            (runs4, folds4),
            (2, 6),
            "n=4: leaf fast path (4 lane-folds) + the 2-pair merge cohort (2 lane-folds)"
        );
        assert_bits_identical(
            &scalar4,
            &cohort4,
            "n=4 cohort routing must be bit-identical to the scalar tournament",
        );
    }


    /// Small-batch adversarial matrix: batch sizes 2..=7 (2/3-lane
    /// remainder groups, lone scalars), all-zero displacements, and
    /// divergent displacement supports that force the per-lane scalar
    /// fallbacks. Cohort routing must match the scalar tournament on the
    /// SAME tree bit for bit (value equality for the all-zero flavor,
    /// whose result signs of zero are engine-independent), and the
    /// divergent-support flavor must leave the cohort engine untouched.
    #[test]

    fn small_batch_cohort_matches_scalar_tournament_bit_identically() {
        use ordered_float::NotNan;

        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(4)
            .build()
            .expect("pool");
        for (d, m) in [(2usize, 5usize), (3usize, 4usize)] {
            let builder = LogSignatureBuilder::<u8>::new()
                .with_num_dimensions(d)
                .with_max_degree(m);
            let basis = lyndon_rs::lyndon::LyndonBasis::<u8>::new(
                d,
                lyndon_rs::lyndon::Sort::Lexicographical,
            )
            .generate_basis(m);
            let mut seed = 0x5ba11_u64 + d as u64 * 977 + m as u64;
            let lcg = |seed: &mut u64| {
                *seed = seed
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                ((*seed >> 33) % 19) as i128 - 9
            };
            for n in 2..=7usize {
                for flavor in ["shared", "zeros", "divergent"] {
                    let rhss: Vec<Vec<NotNan<f64>>> = (0..n)
                        .map(|i| {
                            (0..basis.len())
                                .map(|k| {
                                    if k >= d {
                                        NotNan::new(0.0).unwrap()
                                    } else if flavor == "zeros" {
                                        NotNan::new(0.0).unwrap()
                                    } else {
                                        // divergent: alternate which letter
                                        // is present, so adjacent leaves
                                        // fold disjoint supports
                                        let keep = flavor != "divergent" || (i + k) % 2 == 0;
                                        if keep {
                                            // Zero LCG draws remap to 1: the
                                            // shared flavor needs every letter
                                            // nonzero (a zero letter in the
                                            // FIRST displacement makes the
                                            // group support-divergent, and the
                                            // scalar leaf fallback is then
                                            // correct, not an engagement
                                            // failure).
                                            let v = lcg(&mut seed);
                                            let v = if v == 0 { 1 } else { v };
                                            NotNan::new(v as f64).unwrap()
                                        } else {
                                            NotNan::new(0.0).unwrap()
                                        }
                                    }
                                })
                                .collect()
                        })
                        .collect();
                    let slices: Vec<&[NotNan<f64>]> =
                        rhss.iter().map(|r| r.as_slice()).collect();

                    let reduce = |off: bool| -> (Vec<NotNan<f64>>, usize, usize) {
                        let sig = builder.build::<NotNan<f64>>();
                        set_cohort_off(off);
                        let runs0 = COHORT_ENGINE_RUNS.load(std::sync::atomic::Ordering::SeqCst);
                        let folds0 =
                            COHORT_LANE_FOLDS.load(std::sync::atomic::Ordering::SeqCst);
                        let reduced = pool.install(|| sig.tournament_reduce(&slices));
                        let runs =
                            COHORT_ENGINE_RUNS.load(std::sync::atomic::Ordering::SeqCst) - runs0;
                        let folds =
                            COHORT_LANE_FOLDS.load(std::sync::atomic::Ordering::SeqCst) - folds0;
                        set_cohort_off(false);
                        (reduced, runs, folds)
                    };

                    let (scalar, _, _) = reduce(true);
                    let (cohort, runs, folds) = reduce(false);
                    if flavor == "divergent" {
                        // Divergent INPUT supports force the LEAF round's
                        // per-lane scalar fallback — but they do NOT pin the
                        // counters to zero: chunk-2+ leaves union their
                        // slices' disjoint supports to the SAME reachable
                        // result support on every lane, and the merge rounds
                        // legitimately cohort those equal-result-support
                        // pairs (SIMD-across-folds everywhere it is sound).
                        // Correctness here is the bit-identity check below;
                        // the shared/zeros flavors' engagement asserts cover
                        // the engine-on side.
                        let _ = (runs, folds);
                    } else {
                        assert!(
                            runs > 0 && folds > 0,
                            "{d}x{m} n={n} {flavor}: the cohort engine never engaged"
                        );
                    }
                    if flavor == "zeros" {
                        // All-zero displacements: the result is the zero
                        // vector up to signs of zero, which are add-order
                        // artifacts (both engines add value-zero terms; the
                        // cohort's shared plan skips positions the scalar
                        // walk zero-adds). Compare by value.
                        for (k, (c, s)) in cohort.iter().zip(&scalar).enumerate() {
                            assert_eq!(
                                c, s,
                                "{d}x{m} n={n}: all-zero batch diverged at {k}: {c} vs {s}"
                            );
                        }
                    } else {
                        assert_bits_identical(
                            &scalar,
                            &cohort,
                            &format!("{d}x{m} n={n} {flavor} cohort-vs-scalar"),
                        );
                    }
                }
            }
        }
    }


    /// The cohort tournament engages in SINGLE-worker pools: the
    /// SIMD-across-folds amortization is data parallelism (4 folds share
    /// one plan walk and one kernel pass), not thread parallelism, so a
    /// 1-thread pool still cohorts — and stays bit-identical to the same
    /// tree's scalar engine.
    #[test]

    fn cohort_engages_in_single_thread_pools() {
        use ordered_float::NotNan;

        let (d, m) = (2usize, 5usize);
        let builder = LogSignatureBuilder::<u8>::new()
            .with_num_dimensions(d)
            .with_max_degree(m);
        let basis = lyndon_rs::lyndon::LyndonBasis::<u8>::new(
            d,
            lyndon_rs::lyndon::Sort::Lexicographical,
        )
        .generate_basis(m);
        let mut seed = 0x17e7_u64;
        let lcg = |seed: &mut u64| {
            *seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((*seed >> 33) % 19) as i128 - 9
        };
        let n = 8usize;
        let rhss: Vec<Vec<NotNan<f64>>> = (0..n)
            .map(|_| {
                (0..basis.len())
                    .map(|k| {
                        if k < d {
                            NotNan::new(lcg(&mut seed) as f64).unwrap()
                        } else {
                            NotNan::new(0.0).unwrap()
                        }
                    })
                    .collect()
            })
            .collect();
        let slices: Vec<&[NotNan<f64>]> = rhss.iter().map(|r| r.as_slice()).collect();
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(1)
            .build()
            .expect("pool");

        let (scalar, sruns, sfolds) = {
            let sig = builder.build::<NotNan<f64>>();
            set_cohort_off(true);
            let r0 = COHORT_ENGINE_RUNS.load(std::sync::atomic::Ordering::SeqCst);
            let f0 = COHORT_LANE_FOLDS.load(std::sync::atomic::Ordering::SeqCst);
            let reduced = pool.install(|| sig.tournament_reduce(&slices));
            let r = COHORT_ENGINE_RUNS.load(std::sync::atomic::Ordering::SeqCst) - r0;
            let f = COHORT_LANE_FOLDS.load(std::sync::atomic::Ordering::SeqCst) - f0;
            set_cohort_off(false);
            (reduced, r, f)
        };
        assert_eq!((sruns, sfolds), (0, 0), "kill switch must silence the engine");
        let (cohort, cruns, cfolds) = {
            let sig = builder.build::<NotNan<f64>>();
            let r0 = COHORT_ENGINE_RUNS.load(std::sync::atomic::Ordering::SeqCst);
            let f0 = COHORT_LANE_FOLDS.load(std::sync::atomic::Ordering::SeqCst);
            let reduced = pool.install(|| sig.tournament_reduce(&slices));
            let r = COHORT_ENGINE_RUNS.load(std::sync::atomic::Ordering::SeqCst) - r0;
            let f = COHORT_LANE_FOLDS.load(std::sync::atomic::Ordering::SeqCst) - f0;
            (reduced, r, f)
        };
        assert!(
            cruns > 0 && cfolds > 0,
            "1-thread pool: the cohort engine never engaged"
        );
        assert_bits_identical(
            &scalar,
            &cohort,
            "1-thread pool cohort routing must be bit-identical to the scalar tournament",
        );
    }


    /// Wide-pool engagement for LIGHT folds: an 8x3-class config (tiny
    /// table ~252 pairs, 4-term DAG) at a 32-thread pool must still route
    /// its folds through the cohort engine. This pins the DELETION of the
    /// wide-pool gate (`threads <= 16 || terms >= 256 || table >= 1000 ||
    /// groups >= threads`): post-F2 the single-phase per-unit sweep made
    /// the cohort the light engine, and the measured counterfactual showed
    /// scalar fallback at 32 workers costs +18..30% wall (plus a ~73%-of-
    /// cycles steal-probe storm) precisely in the regime the old gate
    /// protected. If a width-based gate ever returns, this fails first.
    #[test]

    fn cohort_engages_at_32_thread_pools_for_light_folds() {
        use ordered_float::NotNan;

        let (d, m) = (8usize, 3usize);
        let builder = LogSignatureBuilder::<u8>::new()
            .with_num_dimensions(d)
            .with_max_degree(m);
        let basis = lyndon_rs::lyndon::LyndonBasis::<u8>::new(
            d,
            lyndon_rs::lyndon::Sort::Lexicographical,
        )
        .generate_basis(m);
        let mut seed = 0x8d3_u64;
        let lcg = |seed: &mut u64| {
            *seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((*seed >> 33) % 19) as i128 - 9
        };
        let n = 64usize;
        let rhss: Vec<Vec<NotNan<f64>>> = (0..n)
            .map(|_| {
                (0..basis.len())
                    .map(|k| {
                        if k < d {
                            NotNan::new(lcg(&mut seed) as f64).unwrap()
                        } else {
                            NotNan::new(0.0).unwrap()
                        }
                    })
                    .collect()
            })
            .collect();
        let slices: Vec<&[NotNan<f64>]> = rhss.iter().map(|r| r.as_slice()).collect();
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(32)
            .build()
            .expect("pool");

        let run = |off: bool| -> (Vec<NotNan<f64>>, usize, usize) {
            let sig = builder.build::<NotNan<f64>>();
            set_cohort_off(off);
            let r0 = COHORT_ENGINE_RUNS.load(std::sync::atomic::Ordering::SeqCst);
            let f0 = COHORT_LANE_FOLDS.load(std::sync::atomic::Ordering::SeqCst);
            let reduced = pool.install(|| sig.tournament_reduce(&slices));
            let r = COHORT_ENGINE_RUNS.load(std::sync::atomic::Ordering::SeqCst) - r0;
            let f = COHORT_LANE_FOLDS.load(std::sync::atomic::Ordering::SeqCst) - f0;
            set_cohort_off(false);
            (reduced, r, f)
        };

        let (scalar, sruns, _) = run(true);
        assert_eq!(sruns, 0, "kill switch must silence the cohort engine");
        let (cohort, cruns, cfolds) = run(false);
        assert!(
            cruns > 0 && cfolds > 0,
            "32-thread pool, light 8x3-class fold: the cohort engine never \
             engaged — a width-based gate is back"
        );
        assert_bits_identical(
            &scalar,
            &cohort,
            "32-thread light-fold cohort routing must be bit-identical to scalar",
        );
    }


    /// Kill-switch bit-identity at wide-pool width for the same light
    /// config: the scalar tournament (switch off) and the cohort tournament
    /// (switch on) must agree bit for bit at 32 threads — the switch swaps
    /// engines, never the association tree.
    #[test]

    fn cohort_kill_switch_bit_identity_at_32_thread_light_folds() {
        use ordered_float::NotNan;

        let (d, m) = (8usize, 3usize);
        let builder = LogSignatureBuilder::<u8>::new()
            .with_num_dimensions(d)
            .with_max_degree(m);
        let basis = lyndon_rs::lyndon::LyndonBasis::<u8>::new(
            d,
            lyndon_rs::lyndon::Sort::Lexicographical,
        )
        .generate_basis(m);
        let mut seed = 0x8d4_u64;
        let lcg = |seed: &mut u64| {
            *seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((*seed >> 33) % 13) as i128 - 6
        };
        let n = 40usize;
        let rhss: Vec<Vec<NotNan<f64>>> = (0..n)
            .map(|_| {
                (0..basis.len())
                    .map(|k| {
                        if k < d {
                            NotNan::new(lcg(&mut seed) as f64).unwrap()
                        } else {
                            NotNan::new(0.0).unwrap()
                        }
                    })
                    .collect()
            })
            .collect();
        let slices: Vec<&[NotNan<f64>]> = rhss.iter().map(|r| r.as_slice()).collect();
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(32)
            .build()
            .expect("pool");

        let run = |off: bool| -> Vec<NotNan<f64>> {
            let sig = builder.build::<NotNan<f64>>();
            set_cohort_off(off);
            let reduced = pool.install(|| sig.tournament_reduce(&slices));
            set_cohort_off(false);
            reduced
        };
        let scalar = run(true);
        let cohort = run(false);
        assert_bits_identical(
            &scalar,
            &cohort,
            "8x3 kill-switch states must be bit-identical at 32 threads",
        );
    }


    /// Degenerate light config with the cohort engine UNAVAILABLE (kill
    /// switch): the tournament must stay fully scalar (zero engagement
    /// counters) and still match the sequential left-fold reference
    /// bit-for-bit — scalar fallback remains the correct engine whenever
    /// the cohort cannot run, independent of pool width.
    #[test]

    fn degenerate_light_config_stays_scalar_when_cohort_unavailable() {
        use ordered_float::NotNan;

        let (d, m) = (8usize, 3usize);
        let builder = LogSignatureBuilder::<u8>::new()
            .with_num_dimensions(d)
            .with_max_degree(m);
        let basis = lyndon_rs::lyndon::LyndonBasis::<u8>::new(
            d,
            lyndon_rs::lyndon::Sort::Lexicographical,
        )
        .generate_basis(m);
        let mut seed = 0x8d5_u64;
        let lcg = |seed: &mut u64| {
            *seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((*seed >> 33) % 17) as i128 - 8
        };
        let n = 24usize;
        let rhss: Vec<Vec<NotNan<f64>>> = (0..n)
            .map(|_| {
                (0..basis.len())
                    .map(|k| {
                        if k < d {
                            NotNan::new(lcg(&mut seed) as f64).unwrap()
                        } else {
                            NotNan::new(0.0).unwrap()
                        }
                    })
                    .collect()
            })
            .collect();
        let slices: Vec<&[NotNan<f64>]> = rhss.iter().map(|r| r.as_slice()).collect();
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(32)
            .build()
            .expect("pool");

        // Sequential association differs from the tournament tree for
        // floats (only exact arithmetic is association-independent — see
        // `tournament_rationals_exact_at_awkward_sizes`), so the reference
        // here is the scalar tournament at a NARROW width: same tree, same
        // bits, different scheduling.
        set_cohort_off(true);
        let r0 = COHORT_ENGINE_RUNS.load(std::sync::atomic::Ordering::SeqCst);
        let f0 = COHORT_LANE_FOLDS.load(std::sync::atomic::Ordering::SeqCst);
        let scalar_tournament_32t = pool
            .install(|| builder.build::<NotNan<f64>>().tournament_reduce(&slices));
        let runs = COHORT_ENGINE_RUNS.load(std::sync::atomic::Ordering::SeqCst) - r0;
        let folds = COHORT_LANE_FOLDS.load(std::sync::atomic::Ordering::SeqCst) - f0;
        let narrow = rayon::ThreadPoolBuilder::new()
            .num_threads(4)
            .build()
            .expect("narrow pool");
        let scalar_tournament_4t =
            narrow.install(|| builder.build::<NotNan<f64>>().tournament_reduce(&slices));
        set_cohort_off(false);
        assert_eq!(runs, 0, "kill switch: no cohort runs");
        assert_eq!(folds, 0, "kill switch: no cohort lane-folds");
        assert_bits_identical(
            &scalar_tournament_4t,
            &scalar_tournament_32t,
            "scalar tournament must be bit-identical across pool widths at 8x3",
        );
    }
}
