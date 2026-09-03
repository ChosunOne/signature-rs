use lie_rs::LieSeriesGenerator;
use lie_rs::{LieSeries, bch_series_generator::BchSeriesGenerator};
use lyndon_rs::lyndon::LyndonWord;
use lyndon_rs::{
    generators::Generator,
    lyndon::{LyndonBasis, Sort},
};
use ndarray::{ArrayView, Axis, Dimension, RemoveAxis};

use num_traits::{FromPrimitive, One, Zero};
use std::ops::{Mul, Sub};
use std::{
    fmt::Debug,
    hash::Hash,
    ops::{AddAssign, Div, Index, IndexMut, MulAssign, Neg, SubAssign},
    sync::Mutex,
};

use crate::commutator_dag::{CohortLane, CommutatorDag, cohort_capable};

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
/// tournament reduction over the sequential fold chain. Below this the
/// per-fold/batch sequential path is cheaper than the tournament's tree
/// overhead (extra DAG warm-ups, pool coordination, one final merge fold).
const TOURNAMENT_MIN_DISPLACEMENTS: usize = 16;

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
const TOURNAMENT_LEAF_CHUNK: usize = 8;

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
const TOURNAMENT_MAX_LEAVES: usize = 32;

/// Cohort engagement cap for narrow folds: above this many rayon workers,
/// a fold with few commutator terms does not form cohort groups. The
/// cohort's top-level parallelism is its group count (each group runs one
/// 4-lane engine whose stages ramp up internally), while the scalar
/// tournament keeps one task per leaf/pair. Measured crossover on the
/// e2e matrix at n=1000 (32 leaves → 8 groups): at 32 workers the
/// 533-term 2x12 fold still wins (−12.5%) because its engine stages fill
/// the pool, but the 52-term 3x8 fold regresses (+7.5%) — the 8 groups
/// cannot feed 32 workers through ~4-wide per-step stage batches. 256
/// terms sits 5× below and 5× above the two measured data points.
const COHORT_WIDE_POOL_TERMS: usize = 256;

/// Below [`COHORT_WIDE_POOL_TERMS`] terms, cohort groups form only up to
/// this many workers (every measured combo at ≤16 workers passes the
/// no-regression gate; the big cohort wins are at 4–8 workers).
const COHORT_MAX_THREADS: usize = 16;

/// NaN audit for an arbitrary coefficient slice — the
/// `LogSignature::audit_no_nan` check factored off `&self` so the
/// tournament's locally owned accumulators (plain `LieSeries` values, not
/// `LogSignature`s) can run the same per-fold audit.
#[inline]
fn audit_coefficients_no_nan<U: PartialEq>(coefficients: &[U]) {
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
/// it through [`fold_batch_sequential`]'s warm-up phase.
fn fold_one_displacement<T, U>(
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
fn fold_batch_sequential<T, U>(
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
struct SeriesTemplate<U, T>(LieSeries<T, U>);

// SAFETY: see the struct-level safety contract.
unsafe impl<U: Send + Sync, T> Send for SeriesTemplate<U, T> {}
// SAFETY: see the struct-level safety contract.
unsafe impl<U: Send + Sync, T> Sync for SeriesTemplate<U, T> {}

/// Pops a pooled DAG, falling back to a fresh shallow clone of `source`
/// (the caller's compiled plan) when the pool is empty.
fn take_pooled_dag<U>(
    pool: &Mutex<Vec<CommutatorDag<U>>>,
    source: &CommutatorDag<U>,
) -> CommutatorDag<U> {
    pool.lock()
        .unwrap_or_else(|poisoned| poisoned.into_inner())
        .pop()
        .unwrap_or_else(|| source.clone_shallow())
}

/// Returns a finished DAG to the pool for the next fold.
fn return_pooled_dag<U>(pool: &Mutex<Vec<CommutatorDag<U>>>, dag: CommutatorDag<U>) {
    pool.lock()
        .unwrap_or_else(|poisoned| poisoned.into_inner())
        .push(dag);
}

/// One cohort lane's evolving leaf-round state. `U`-only (no `T`), so the
/// vector of states crosses into nested rayon tasks for the warm-up and
/// scalar-step phases — the `LieSeries` instances those tasks build live
/// and die inside their own task (non-last `Arc` drops, the tournament's
/// [`SeriesTemplate`] safety contract).
struct CohortLeafState<U> {
    /// The lane's accumulator (public basis order).
    acc: Vec<U>,
    /// The lane's private DAG (pooled; `Some` until returned).
    dag: Option<CommutatorDag<U>>,
    /// Chunk range of this lane's leaf (indices into the batch's `rhss`).
    start: usize,
    len: usize,
    /// Folds already executed.
    done: usize,
}

impl<U> CohortLeafState<U> {
    fn rest<'a>(&self, rhss: &'a [&'a [U]]) -> &'a [&'a [U]] {
        &rhss[self.start + self.done..self.start + self.len]
    }

    fn is_done(&self) -> bool {
        self.done >= self.len
    }
}

/// The tournament leaf round's COHORT group: folds four adjacent leaf
/// chunks as one 4-lane cohort — each lane warms up its chunk through the
/// exact scalar engine (`fold_one_displacement` steps under the
/// [`batch_rest_eligible`] check, the sequential driver's own decision
/// sequence), and once every lane is batch-steady with IDENTICAL atom
/// supports, the lanes' remaining folds run as ONE SoA cohort dispatch
/// ([`CommutatorDag::fold_batch_cohort`]) — the same reduction tree, the
/// same per-lane fold sequence, 4 folds per plan walk.
///
/// # Soundness of the shared plan (the eligibility rule)
///
/// The fold plan (decomposition lists, scatter indices, gating, BCH
/// weights) is a pure function of the basis and the atom supports — data
/// independent. A cohort therefore shares one plan exactly when its lanes'
/// atom supports are identical and each lane's DAG state is the built
/// fixed point for them ([`batch_rest_eligible`]'s `batch_eligible`, which
/// also pins the accumulator support to the reachable set — the dense
/// steady state). Under that condition the per-lane gating is vacuous
/// (every lane's tickets are the same precomputed list), and the cohort
/// executes the UNION plan every lane would execute alone, bit for bit.
/// Lanes whose supports diverge mid-warm-up (exact-cancellation corners,
/// per-leaf displacement-support differences) fall back to the scalar
/// per-lane engine — never to a shared plan they don't match.
fn cohort_leaf_group<T, U>(
    template: &SeriesTemplate<U, T>,
    source_dag: &CommutatorDag<U>,
    dag_pool: &Mutex<Vec<CommutatorDag<U>>>,
    basis_len: usize,
    rhss: &[&[U]],
    chunk: usize,
    first_leaf: usize,
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

    // Phase 1: warm each lane's chunk up to the batch-steady state through
    // the sequential engine's exact decision sequence (check, fold one,
    // re-check). The lane's series lives only inside its nested task.
    let mut lanes: Vec<CohortLeafState<U>> = (0..4)
        .into_par_iter()
        .map(|l| {
            let start = (first_leaf + l) * chunk;
            let len = chunk.min(rhss.len() - start);
            let mut series = template.0.clone();
            series.coefficients = vec![U::default(); basis_len];
            let mut dag = take_pooled_dag(dag_pool, source_dag);
            let mut done = 0usize;
            while done < len {
                if batch_rest_eligible(&template.0, &series.coefficients, &dag, &rhss[start + done..start + len])
                    .is_some()
                {
                    break;
                }
                fold_one_displacement(&mut series, &mut dag, rhss[start + done]);
                done += 1;
            }
            CohortLeafState {
                acc: series.coefficients,
                dag: Some(dag),
                start,
                len,
                done,
            }
        })
        .collect();

    // Phase 2: drive the lanes to completion. All-steady lanes with
    // identical supports run the cohort tail; anything else scalar-steps
    // one fold per round (bit-identical per lane to the sequential engine —
    // per-fold arithmetic is value- and bit-equal to the batch's per-fold
    // step; only the plan reuse differs).
    loop {
        let active: Vec<usize> = (0..lanes.len()).filter(|&l| !lanes[l].is_done()).collect();
        if active.is_empty() {
            break;
        }
        if active.len() == 1 {
            // A lone lane gets the plain sequential engine (the tournament
            // would gain nothing from a 1-lane cohort).
            let l = active[0];
            let mut series = template.0.clone();
            series.coefficients = std::mem::take(&mut lanes[l].acc);
            let mut dag = lanes[l].dag.take().expect("lane dag present");
            fold_batch_sequential(&mut series, &mut dag, lanes[l].rest(rhss));
            return_pooled_dag(dag_pool, dag);
            lanes[l].acc = series.coefficients;
            lanes[l].done = lanes[l].len;
            break;
        }
        let checks: Vec<Option<(Vec<usize>, Vec<usize>)>> = active
            .iter()
            .map(|&l| {
                let dag = lanes[l].dag.as_ref().expect("lane dag present");
                batch_rest_eligible(&template.0, &lanes[l].acc, dag, lanes[l].rest(rhss))
            })
            .collect();
        let shared_plan = checks[0].is_some()
            && checks.iter().all(|c| c == &checks[0]);
        if shared_plan {
            // Cohort tail: one shared plan, every lane's remaining folds.
            let (a_nonzero, b_nonzero) = checks[0].clone().expect("checked above");
            let lead = active[0];
            let mut lead_dag = lanes[lead].dag.take().expect("lane dag present");
            let mut cohort_lanes: Vec<CohortLane<U>> = active
                .iter()
                .map(|&l| CohortLane {
                    acc: std::mem::take(&mut lanes[l].acc),
                    rhss: lanes[l].rest(rhss),
                })
                .collect();
            lead_dag.fold_batch_cohort(&template.0, &mut cohort_lanes, &a_nonzero, &b_nonzero);
            for (i, &l) in active.iter().enumerate() {
                lanes[l].acc = std::mem::take(&mut cohort_lanes[i].acc);
                lanes[l].done = lanes[l].len;
            }
            return_pooled_dag(dag_pool, lead_dag);
            for &l in active.iter().skip(1) {
                return_pooled_dag(dag_pool, lanes[l].dag.take().expect("lane dag present"));
            }
            continue;
        }
        // Plans diverge (or some lane is still warming up): one scalar fold
        // per active lane, nested-parallel. A lane that was already steady
        // stays steady (value-only changes); the group converges or falls
        // back to scalar for good.
        lanes
            .par_iter_mut()
            .enumerate()
            .filter(|(l, lane)| active.contains(l) && !lane.is_done())
            .for_each(|(_l, lane)| {
                let mut series = template.0.clone();
                series.coefficients = std::mem::take(&mut lane.acc);
                let dag = lane.dag.as_mut().expect("lane dag present");
                fold_one_displacement(&mut series, dag, rhss[lane.start + lane.done]);
                lane.acc = series.coefficients;
                lane.done += 1;
            });
    }

    // Return every remaining dag and collect the lanes' results in leaf
    // order.
    for lane in &mut lanes {
        if let Some(dag) = lane.dag.take() {
            return_pooled_dag(dag_pool, dag);
        }
    }
    lanes.into_iter().map(|lane| lane.acc).collect()
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
fn cohort_merge_group<T, U>(
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
fn scalar_merge_fold<T, U>(
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

/// Builder for constructing log signatures from path data.
///
/// The log signature is a mathematical transform that captures the geometry
/// of a path by expressing it as a series in the free Lie algebra. This builder
/// allows configuring the parameters and constructing log signatures from various inputs.
pub struct LogSignatureBuilder<T> {
    /// The maximum degree of terms to include in the log signature computation.
    pub max_degree: usize,
    /// The Lyndon basis configuration for the underlying algebra.
    pub lyndon_basis: LyndonBasis<T>,
}

impl<T: Debug> Debug for LogSignatureBuilder<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("LogSignatureBuilder")
            .field("max_degree", &self.max_degree)
            .field("lyndon_basis", &self.lyndon_basis)
            .finish()
    }
}

impl<T> Copy for LogSignatureBuilder<T> {}

impl<T> Clone for LogSignatureBuilder<T> {
    fn clone(&self) -> Self {
        *self
    }
}

impl<T> Default for LogSignatureBuilder<T> {
    fn default() -> Self {
        Self {
            max_degree: usize::default(),
            lyndon_basis: LyndonBasis::default(),
        }
    }
}

impl<T> LogSignatureBuilder<T> {
    /// Creates a new log signature builder with default settings.
    #[must_use]
    pub fn new() -> Self {
        Self {
            ..Default::default()
        }
    }

    /// Sets the maximum degree of terms to include in the log signature.
    ///
    /// Higher degrees capture more complex geometric features but increase
    /// computational complexity exponentially.
    #[must_use]
    pub fn with_max_degree(mut self, max_degree: usize) -> Self {
        self.max_degree = max_degree;
        self
    }

    /// Sets the number of dimensions for the path data.
    ///
    /// This determines the size of the generator alphabet and should match
    /// the dimensionality of the input path data.
    #[must_use]
    pub fn with_num_dimensions(mut self, num_dimensions: usize) -> Self {
        self.lyndon_basis.alphabet_size = num_dimensions;
        self
    }

    /// Returns the maximum degree setting.
    #[must_use]
    pub fn max_degree(&self) -> usize {
        self.max_degree
    }

    /// Returns the number of dimensions setting.
    #[must_use]
    pub fn num_dimensions(&self) -> usize {
        self.lyndon_basis.alphabet_size
    }
}

impl<T: Debug + Clone + Eq + Hash + Ord + Generator + Send + Sync> LogSignatureBuilder<T> {
    /// Builds an empty log signature with the configured parameters.
    ///
    /// The resulting log signature has the proper basis structure but with
    /// all coefficients set to zero.
    #[must_use]
    pub fn build<
        U: Clone
            + Default
            + AddAssign
            + Div<Output = U>
            + FromPrimitive
            + Neg<Output = U>
            + One
            + Ord
            + Hash
            + Send
            + Sync
            + Zero
            + SubAssign
            + MulAssign
            + Send
            + Sync,
    >(
        &self,
    ) -> LogSignature<T, U> {
        let bch_basis = LyndonBasis::<u8>::new(2, Sort::Lexicographical);
        let bch_series = BchSeriesGenerator::new(bch_basis, self.max_degree).generate_lie_series();
        let dag = CommutatorDag::from_bch_series(&bch_series);
        let basis = self.lyndon_basis.generate_basis(self.max_degree);
        let coefficients = vec![U::default(); basis.len()];
        let series = LieSeries::<T, U>::new(basis, coefficients);
        LogSignature::<T, U> {
            series,
            bch_series,
            dag,
        }
    }

    /// Computes the log signature of a path from multidimensional array data.
    ///
    /// The path should be provided as a 2D array where each row represents a point
    /// and each column represents a coordinate dimension. The log signature is
    /// computed incrementally over consecutive path segments.
    ///
    /// For paths long enough to qualify for the batch driver's tournament
    /// reduction (see [`LogSignature::concatenate_batch_coefficients`]), the
    /// folds are reassociated along a balanced tree — a shape that depends
    /// only on the number of displacements, never on the thread count — so
    /// f32/f64 rounding may differ by a few ulps from strictly sequential
    /// folding, while exact coefficient types remain bit-identical.
    #[must_use]
    pub fn build_from_path<
        D: Dimension + RemoveAxis,
        U: Clone
            + Default
            + AddAssign
            + Div<Output = U>
            + FromPrimitive
            + Neg<Output = U>
            + One
            + Ord
            + Hash
            + Send
            + Sync
            + Zero
            + SubAssign
            + MulAssign
            + Sub<Output = U>
            + 'static,
    >(
        &self,
        path: &ArrayView<U, D>,
    ) -> LogSignature<T, U> {
        let mut log_sig = self.build();

        // One flat displacement buffer (each row is the elementwise
        // difference of two consecutive path points); the batch's
        // displacement slices borrow it. The per-window `Vec`s this
        // replaces (an owned difference array plus a `Vec` per segment)
        // were a measurable slice of the small-grid e2e's allocator
        // traffic.
        let rows = path.len_of(Axis(0));
        let n_displacements = rows.saturating_sub(1);
        let row_len = path.len().checked_div(rows).unwrap_or(0);
        let mut displacements = vec![U::default(); n_displacements * row_len];
        for (wi, window) in path.axis_windows(Axis(0), 2).into_iter().enumerate() {
            let start = window.index_axis(Axis(0), 0);
            let end = window.index_axis(Axis(0), 1);
            let dst = &mut displacements[wi * row_len..(wi + 1) * row_len];
            for (c, (e, s)) in dst.iter_mut().zip(end.iter().zip(start.iter())) {
                *c = e.clone() - s.clone();
            }
        }
        let slices: Vec<&[U]> = displacements
            .chunks(row_len.max(1))
            .take(n_displacements)
            .collect();
        log_sig.concatenate_batch_coefficients(&slices);

        log_sig
    }
}

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
    dag: CommutatorDag<U>,
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
    /// Two strategies, selected per call:
    ///
    /// - **Sequential** (fewer than `TOURNAMENT_MIN_DISPLACEMENTS`
    ///   displacements, or a single-worker pool): the folds run one-by-one
    ///   until the accumulator's support is the full basis and the node
    ///   lists are steady for it; the remaining folds — which must share
    ///   one displacement support — then run as ONE continuous batch
    ///   dispatch (see `CommutatorDag::fold_batch`), which keeps the
    ///   accumulator class-ordered and reuses one plan and one job table
    ///   for the whole batch, turning the per-fold gather/accumulate
    ///   phases into parallel in-graph stages. This path is bit-identical
    ///   to folding each displacement with
    ///   [`Self::concatenate_assign_coefficients`].
    /// - **Tournament** (at least `TOURNAMENT_MIN_DISPLACEMENTS`
    ///   displacements and `rayon::current_num_threads() > 1`): BCH
    ///   concatenation is associative, so the folds over `rhss` form a
    ///   reduction whose tree shape is free to choose. A balanced binary
    ///   tournament over contiguous leaf chunks (`Self::tournament_reduce`)
    ///   makes each round's folds independent — leaf chunks fold disjoint
    ///   displacement ranges, merge rounds fold adjacent chunk results — so
    ///   whole rounds run in parallel, and only the merge rounds' folds
    ///   (at most `TOURNAMENT_MAX_LEAVES - 1` of them) fold a fully dense
    ///   accumulator instead of all `n` folds in one chain. The tree shape
    ///   depends ONLY on `rhss.len()` — never on the thread count — so a
    ///   given input always reduces along the same tree and results are
    ///   reproducible at any pool size.
    ///
    /// Rounding caveat: f32/f64 accumulation is not associative, so the
    /// tournament's reassociated tree is NOT bit-identical to the
    /// sequential path — observed worst-case shifts are a few ulps
    /// (~1e-13 absolute on coefficients of magnitude ~1e2 in adversarial
    /// test batches). Within one strategy results are bit-stable: the
    /// sequential path is bit-identical to folding each displacement with
    /// `Self::concatenate_assign_coefficients`, and the tournament's
    /// tree shape depends only on `rhss.len()`, so a given input reduces
    /// along the same tree at any pool size. Exact coefficient types
    /// (e.g. `Ratio<i128>`) are insensitive to association order and
    /// produce identical results on every path.
    ///
    /// No nested-recursion guard is needed: the tournament's leaves call
    /// the free [`fold_batch_sequential`] engine directly — the engine
    /// contains no driver-selection logic, so it cannot re-enter this
    /// method — and neither it nor `concatenate_assign_coefficients`
    /// re-enters this driver.
    pub fn concatenate_batch_coefficients(&mut self, rhss: &[&[U]])
    where
        U: Send + Sync,
    {
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
    fn concatenate_batch_sequential(&mut self, rhss: &[&[U]])
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
    fn tournament_reduce(&self, rhss: &[&[U]]) -> Vec<U>
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

        // Adaptive leaf chunk: large enough that each leaf's sequential
        // warm-up amortizes over its chunk, capped at ~32 leaves to keep
        // in-flight oversubscription modest — see `TOURNAMENT_MAX_LEAVES`.
        let chunk =
            TOURNAMENT_LEAF_CHUNK.max(rhss.len().div_ceil(TOURNAMENT_MAX_LEAVES));

        // The cohort capability: raw repr-transparent float coefficients +
        // AVX2 + kill switch (see `cohort_capable`). False keeps every
        // round on the scalar paths below — the pre-cohort tournament,
        // bit for bit.
        let cohort_on = cohort_capable::<U>();
        // Wide-pool engagement rule (see `COHORT_WIDE_POOL_TERMS`): with a
        // wide pool, only heavy folds form cohort groups outright; light
        // folds need enough groups in the round to cover the pool on their
        // own. `groups` is per round — the leaf round's group count is
        // `leaf_count/4`, a merge round's is `ceil(n_pairs/4)` — so the
        // round gates below re-check with their own counts.
        let threads = rayon::current_num_threads().max(1);
        let terms = source_dag.structure.terms.len();
        let cohort_round_ok = |groups: usize| {
            threads <= COHORT_MAX_THREADS
                || terms >= COHORT_WIDE_POOL_TERMS
                || groups >= threads
        };

        // A leaf folded through the sequential batch engine (the scalar
        // leaf task body — also the partial group's path).
        let scalar_leaf = |leaf: usize| -> Vec<U> {
            let start = leaf * chunk;
            let end = (start + chunk).min(rhss.len());
            let mut series = template.0.clone();
            series.coefficients = vec![U::default(); basis_len];
            let mut dag = take_pooled_dag(&dag_pool, source_dag);
            fold_batch_sequential(&mut series, &mut dag, &rhss[start..end]);
            return_pooled_dag(&dag_pool, dag);
            series.coefficients
        };

        // Leaf round: chunks fold their disjoint displacement ranges
        // independently (the tree's bottom level). Full groups of four run
        // as cohort groups (one task per group: per-lane scalar warm-ups +
        // one shared-plan SoA tail); the tail group of one-to-three leaves
        // and every lane of a non-capable configuration run the scalar
        // path.
        let leaf_count = rhss.len().div_ceil(chunk);
        let mut level: Vec<Vec<U>> =
            if cohort_on && leaf_count >= 4 && cohort_round_ok(leaf_count / 4) {
            (0..leaf_count.div_ceil(4))
                .into_par_iter()
                .flat_map_iter(|g| {
                    let first = g * 4;
                    let last = (first + 4).min(leaf_count);
                    if last - first == 4 {
                        cohort_leaf_group(
                            template,
                            source_dag,
                            &dag_pool,
                            basis_len,
                            rhss,
                            chunk,
                            first,
                        )
                    } else {
                        (first..last)
                            .into_par_iter()
                            .map(scalar_leaf)
                            .collect::<Vec<_>>()
                    }
                })
                .collect()
        } else {
            (0..leaf_count).into_par_iter().map(scalar_leaf).collect()
        };

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
            level = if cohort_on && cohort_round_ok(n_pairs.div_ceil(4)) {
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

#[cfg(test)]
mod test {
    use ndarray::{Array2, array};
    use num_rational::Ratio;
    use num_traits::ToPrimitive;
    use ordered_float::NotNan;
    use rstest::rstest;

    use super::*;
    use crate::commutator_dag::{COHORT_ENGINE_RUNS, set_cohort_off};

    #[test]
    fn dag_node_lists_cover_buffer_support() {
        use ordered_float::NotNan;

        let (d, m) = (3usize, 5usize);
        let builder: LogSignatureBuilder<u8> = LogSignatureBuilder::<u8>::new()
            .with_num_dimensions(d)
            .with_max_degree(m);
        let mut log_sig: LogSignature<u8, NotNan<f64>> = builder.build::<NotNan<f64>>();

        // Dense accumulator, then several letter-displacement folds: the
        // first fold goes through the collecting rebuild, later folds reuse
        // the fixed-point lists.
        let basis: Vec<LyndonWord<u8>> =
            lyndon_rs::lyndon::LyndonBasis::<u8>::new(d, lyndon_rs::lyndon::Sort::Lexicographical)
                .generate_basis(m);
        let mut seed: u64 = 0xabcd_u64;
        let lcg = |seed: &mut u64| {
            *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
            ((*seed >> 33) as i64) as f64
        };
        let acc: Vec<NotNan<f64>> = (0..basis.len())
            .map(|_| NotNan::new(lcg(&mut seed)).unwrap())
            .collect();
        log_sig.series.coefficients.clone_from(&acc);

        for _ in 0..4 {
            let disp: Vec<NotNan<f64>> = (0..basis.len())
                .map(|k| NotNan::new(if k < d { lcg(&mut seed) } else { 0.0 }).unwrap())
                .collect();
            let mut disp_sig: LogSignature<u8, NotNan<f64>> = builder.build::<NotNan<f64>>();
            disp_sig.series.coefficients.clone_from(&disp);
            log_sig.concatenate_assign(&disp_sig);

            let dag: &CommutatorDag<NotNan<f64>> = &log_sig.dag;
            let cutoff: usize = log_sig
                .series
                .commutator_basis
                .iter()
                .take_while(|w| w.degree() != m)
                .count();
            for k in 2..dag.structure.nodes.len() {
                let buffer: &Vec<NotNan<f64>> = &dag.buffers[k - 2];
                let list: &Vec<usize> = &dag.nonzeros[k - 2];
                let mut sorted: Vec<usize> = list.clone();
                sorted.sort_unstable();
                sorted.dedup();
                assert_eq!(
                    sorted.len(),
                    list.len(),
                    "node {k}: list contains duplicate indices"
                );
                for (i, c) in buffer.iter().enumerate().take(cutoff) {
                    assert!(
                        c.is_zero() || list.contains(&i),
                        "node {k}: index {i} is non-zero but not listed"
                    );
                }
            }
        }
    }

    #[rstest]
    #[case(
        3,
        3,
        array![
            [0.0, 0., 0.],
            [1.0, 2.0, 3.0],
        ],
        vec![
            1.000000,
            2.000000,
            3.000000,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
        ])]
    #[case(
        3,
        3,
        array![
            [0.0, 0., 0.],
            [1.0, 2.0, 3.0],
            [6.0, 5.0, 4.0],
        ],
        vec![
            6.,
            5.,
            4.,
            -3.5,
            -7.,
            -3.5,
            2.333333,
            4.666667,
            -0.583333,
            3.5,
            0.,
            2.333333,
            0.583333,
            1.166667,
        ])]
    #[case(
        3,
        3,
        array![
            [0.0, 0., 0.],
            [1.0, 2.0, 3.0],
            [6.0, 5.0, 4.0],
            [7.0, 8.0, 9.0],
            [12.0, 11.0, 10.0]
        ],
        vec![
            12.000000,
            11.000000,
            10.000000,
            -6.500000,
            -13.000000,
            -6.500000,
            -1.166667,
            -2.333333,
            4.416667,
            6.500000,
            16.500000,
            15.333333,
            -4.416667,
            7.666667,
        ])]
    #[case(
        3,
        3,
        array![
            [0.0, 0., 0.],
            [-0.077, 0.042, -0.067],
            [-0.154, 0.675, 0.006],
            [0.916, 1.177, -0.139],
            [1.095, 0.823, -0.261]
        ],
        vec![
           1.095000,
           0.823000,
          -0.261000,
          -0.690006,
          -0.040871,
          -0.124105,
           0.098690,
          -0.004304,
           0.146613,
           0.024960,
           0.044713,
          -0.000215,
          -0.038903,
           0.001568,
        ])]
    #[case(
        3,
        4,
        array![
            [0., 0., 0.],
            [1.0, 2.0, 3.0],
            [6.0, 5.0, 4.0],
        ],
        vec![
            6.000000,
            5.000000,
            4.000000,
            -3.500000,
            -7.000000,
            -3.500000,
            2.333333,
            4.666667,
            -0.583333,
            3.500000,
            0.000000,
            2.333333,
            0.583333,
            1.166667,
            1.458333,
            2.916667,
            -3.791667,
            -3.208333,
            -12.250000,
            -9.333333,
            -1.458333,
            1.750000,
            3.500000,
            2.916667,
            -3.791667,
            6.708333,
            2.625000,
            7.291667,
            1.750000,
            1.750000,
            -3.208333,
            0.875000,
        ]
    )]
    #[case(
        3,
        4,
        array![
            [   0.000,    0.000,    0.000],
            [   1.000,    2.000,    3.000],
            [   6.000,    5.000,    4.000],
            [   7.000,    8.000,    9.000],
            [  12.000,   11.000,   10.000]
        ],
        vec![
            12.000000,
            11.000000,
            10.000000,
            -6.500000,
            -13.000000,
            -6.500000,
            -1.166667,
            -2.333333,
            4.416667,
            6.500000,
            16.500000,
            15.333333,
            -4.416667,
            7.666667,
            -0.041667,
            -0.083333,
            1.208333,
            2.291667,
            4.750000,
            4.666667,
            0.041667,
            -2.250000,
            -4.500000,
            -11.083333,
            -4.291667,
            -12.291667,
            -19.875000,
            -22.208333,
            -13.250000,
            -2.250000,
            7.791667,
            -6.625000,
        ]
    )]
    #[case(
        4,
        4,
        array![
            [   0.000,    0.000,    0.000,    0.000],
            [   1.000,    2.000,    3.000,    4.000],
            [   6.000,    5.000,    4.000,    3.000],
            [   7.000,    8.000,    9.000,    8.000],
            [  12.000,   11.000,   10.000,    9.000],
        ],
        vec![
            12.000000,
            11.000000,
            10.000000,
            9.000000,
            -6.500000,
            -13.000000,
            -13.500000,
            -6.500000,
            -7.000000,
            -0.500000,
            -1.166667,
            -2.333333,
            10.500000,
            4.416667,
            6.500000,
            18.583333,
            16.500000,
            15.333333,
            26.666667,
            5.166667,
            20.833333,
            7.750000,
            -4.416667,
            6.166667,
            7.666667,
            17.500000,
            6.250000,
            0.833333,
            8.333333,
            -6.083333,
            -0.041667,
            -0.083333,
            -14.625000,
            1.208333,
            2.291667,
            -19.125000,
            4.750000,
            4.666667,
            -23.625000,
            21.083333,
            12.916667,
            2.375000,
            0.041667,
            8.083333,
            -2.250000,
            -4.500000,
            -18.750000,
            -11.083333,
            -4.291667,
            -24.500000,
            4.583333,
            8.416667,
            5.750000,
            1.541667,
            -12.291667,
            -19.875000,
            -9.958333,
            -22.208333,
            -13.250000,
            -25.375000,
            1.083333,
            -21.458333,
            9.125000,
            -14.833333,
            -16.666667,
            -9.500000,
            -38.250000,
            -31.791667,
            -19.000000,
            -12.416667,
            -22.458333,
            -4.500000,
            -2.250000,
            -11.500000,
            7.791667,
            -8.166667,
            22.666667,
            10.166667,
            4.250000,
            -6.625000,
            -15.250000,
            -8.166667,
            7.916667,
            -18.458333,
            -9.500000,
            -14.583333,
            -3.333333,
            -5.125000,
            6.708333,
            -2.166667,
        ]
    )]
    fn test_log_sig_builder_from_path(
        #[case] num_dimensions: usize,
        #[case] max_degree: usize,
        #[case] path: Array2<f64>,
        #[case] expected_coefficients: Vec<f64>,
    ) {
        use lyndon_rs::generators::ENotation;
        use ndarray::s;
        use ordered_float::NotNan;

        let builder = LogSignatureBuilder::<ENotation>::new()
            .with_max_degree(max_degree)
            .with_num_dimensions(num_dimensions);
        let path = path.mapv(|v| NotNan::new(v).expect("value to be a number"));
        let log_sig = builder.build_from_path(&path.slice(s![.., ..]));
        for (i, c) in log_sig.series.commutator_basis.iter().enumerate() {
            println!("{i}: {} \t {c}", log_sig.series.basis[i]);
        }
        for (i, c) in log_sig.series.coefficients.iter().enumerate() {
            println!("[{i}]: {c}");
        }
        assert_eq!(
            log_sig.series.coefficients.len(),
            expected_coefficients.len()
        );
        for (i, &c) in expected_coefficients.iter().enumerate() {
            assert!(
                (c - log_sig.series.coefficients[i].to_f64().unwrap()).abs() < 0.0001,
                "{i}: {} != {c}",
                log_sig.series.coefficients[i].to_f64().unwrap()
            );
        }
    }

    #[test]
    fn test_log_sig_concat() {
        let builder = LogSignatureBuilder::<u8>::new()
            .with_num_dimensions(2)
            .with_max_degree(3);
        let mut a = builder.build();
        let mut b = builder.build();
        a.series.coefficients = [1, 2, 3, 4, 5].map(Ratio::from_integer).to_vec();
        b.series.coefficients = [6, 7, 8, 9, 10].map(Ratio::from_integer).to_vec();
        let c = a.concatenate(&b);
        let expected_coefficients = [
            Ratio::new(7, 1),
            Ratio::new(9, 1),
            Ratio::new(17, 2),
            Ratio::new(121, 12),
            Ratio::new(185, 12),
        ];
        assert_eq!(c.series.coefficients, expected_coefficients);
    }

    #[rstest]
    #[case(5, 5, array![
        [1., 2., 3., 4., 5.],
        [6., 7., 8., 9., 10.],
        [11., 12., 13., 14., 15.],
        [16., 17., 18., 19., 20.],
    ])]
    fn test_log_sig_concat_from_path(
        #[case] num_dimensions: usize,
        #[case] max_degree: usize,
        #[case] path: Array2<f64>,
    ) {
        use ndarray::s;
        use ordered_float::NotNan;

        let builder = LogSignatureBuilder::<u8>::new()
            .with_num_dimensions(num_dimensions)
            .with_max_degree(max_degree);
        let path = path.mapv(|v| NotNan::new(v).expect("value to be a number"));

        dbg!(&path.slice(s![0..2, ..]));

        let mut concatenated_log_sig = builder.build_from_path(&path.slice(s![0..=1, ..]));
        let log_sig_2 = builder.build_from_path(&path.slice(s![1.., ..]));

        concatenated_log_sig.concatenate_assign(&log_sig_2);

        let log_sig = builder.build_from_path(&path.slice(s![.., ..]));

        for (concat_c, full_path_c) in concatenated_log_sig
            .series
            .coefficients
            .iter()
            .zip(log_sig.series.coefficients.iter())
        {
            assert_eq!(concat_c, full_path_c);
        }
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

    /// The length gate selects the engine observationally. Below
    /// `TOURNAMENT_MIN_DISPLACEMENTS` — and in any single-worker pool —
    /// the driver must return the SEQUENTIAL engine's exact bits; at or
    /// above the threshold in a multi-thread pool it must return the
    /// tournament tree's exact bits. Both observables are sharp because
    /// the engines' f64 outputs differ by reassociation dust on this data
    /// (measured: reducing the same 15 or 16 displacements along the
    /// tournament tree vs sequentially moves 9 coefficients at the bit
    /// level, ~1e-13 absolute), so exact equality below the gate FAILS if
    /// the tournament is wrongly taken, and exact equality with a direct
    /// `tournament_reduce` call at the gate FAILS if the sequential
    /// engine is wrongly taken:
    /// - 0 displacements: a no-op (the accumulator stays zero).
    /// - 1 displacement: identical to one `concatenate_assign_coefficients`.
    /// - 15 displacements, 4-thread pool: bit-identical to per-fold
    ///   folding (the tournament's dust would fail this).
    /// - 16 displacements, 4-thread pool: bit-identical to a direct
    ///   `tournament_reduce` call (the sequential engine's dust would
    ///   fail this) and within reassociation tolerance of per-fold
    ///   folding, with at least one coefficient bit-different — the
    ///   measured tournament-ran observable (a bit-identical result here
    ///   would mean the two engines coincide on this data and the data
    ///   must be re-derived).
    /// - 16 displacements, 1-thread pool: bit-identical to per-fold
    ///   folding (single-worker pools never take the tournament).
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

        // 15 displacements in a multi-thread pool: sequential bits.
        let seq15 = fold_all(&rhss[..15]);
        let slices15: Vec<&[NotNan<f64>]> =
            rhss[..15].iter().map(|r| r.as_slice()).collect();
        let mut got15 = builder.build::<NotNan<f64>>();
        pool4.install(|| got15.concatenate_batch_coefficients(&slices15));
        assert_bits_identical(
            &seq15,
            &got15.series.coefficients,
            "15 displacements must take the sequential path",
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

        // 16 displacements in a single-worker pool: sequential bits.
        let mut got1 = builder.build::<NotNan<f64>>();
        pool1.install(|| got1.concatenate_batch_coefficients(&slices16));
        assert_bits_identical(
            &seq16,
            &got1.series.coefficients,
            "16 displacements in a 1-thread pool must take the sequential path",
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

    /// The raw-float fast path does not panic on per-operation overflow (see
    /// `lie_rs::raw_mul`'s NaN policy); the fold audits the accumulator once
    /// per step and fails loudly instead of persisting a NaN through the
    /// `NotNan` invariant. All-`MAX` accumulator and letter displacement:
    /// every product overflows to `inf`, and the fused two-orientation term
    /// computes `inf + (-inf) = NaN` deterministically.
    #[test]
    #[should_panic(expected = "log-signature coefficients overflowed to NaN")]
    fn fold_audits_nan_after_overflow() {
        use ordered_float::NotNan;

        let (d, m) = (2usize, 3usize);
        let builder = LogSignatureBuilder::<u8>::new()
            .with_num_dimensions(d)
            .with_max_degree(m);
        let mut log_sig = builder.build::<NotNan<f64>>();
        log_sig.series.coefficients =
            vec![NotNan::new(f64::MAX).unwrap(); log_sig.series.coefficients.len()];
        let mut seg = builder.build::<NotNan<f64>>();
        seg.series.coefficients = (0..seg.series.coefficients.len())
            .map(|k| {
                NotNan::new(if k < d { f64::MAX } else { 0.0 }).unwrap()
            })
            .collect();
        log_sig.concatenate_assign(&seg);
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
            for nz in &dag.nonzeros {
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
}
