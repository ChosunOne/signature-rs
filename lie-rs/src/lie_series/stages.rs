use super::*;
use super::gating::KernelGating;

/// One stage of a class-ordered fold schedule: a kernel sweep level (packs
/// of `(job, unit)` tasks over one level's jobs) or a caller-defined block
/// stage (one atomic claim per block index, dispatched through a shared
/// closure — the gather/accumulate phases of a batched fold). A batch of
/// folds is a linear chain of stages; every slot waits at every stage
/// boundary, so a stage's reads observe exactly the writes of all earlier
/// stages (each completion is published with `Release`, observed with
/// `Acquire`).
pub struct ClassBatchStage<'a, U> {
    /// Sweep stages: `(job index within the stage, unit index in that
    /// job's gating)` tasks in sweep order. Block stages: unused.
    pub(super) tasks: Arc<[(u32, u32)]>,
    /// Balanced pack ranges into `tasks` (sweep stages) or one block per
    /// pack (block stages).
    pub(super) packs: Arc<[(usize, usize)]>,
    /// Sweep stages: per-job gating.
    pub(super) gateways: Vec<KernelGating>,
    /// Sweep stages: the stage's jobs, indexed by the tasks' job ids.
    pub(super) jobs: &'a [KernelJob<'a, U>],
    /// Block stages: `block(block_index, fold_index)` per claimed task.
    /// Blocks write disjoint ranges (the caller's contract); the stage
    /// counters order all cross-stage dependencies.
    pub(super) block: Option<&'a (dyn Fn(usize, usize) + Send + Sync)>,
    /// Block stages: which fold of the batch this stage serves (the
    /// block task reads the fold's inputs through it). Sweep stages: 0.
    pub(super) fold: usize,
    /// Sweep stages: run the 4-lane SoA cohort sweep instead of the scalar
    /// one. A cohort stage's jobs carry SoA-interleaved operand/result
    /// pointers (the four lanes' values of class position `p` live at
    /// `[p*4 .. p*4+4]`, lane `l` at `p*4 + l`), and every other plan
    /// input (tasks, packs, gateways, support lists) is lane-invariant —
    /// see `CommutatorDag::fold_batch_cohort` for how such a plan is
    /// built.
    pub(super) cohort: bool,
}

/// The per-block task/pack shape shared by every block stage of one
/// batch: `blocks` disjoint-range tasks, one pack per task. Built once
/// per batch, then each block stage clones the `Arc`s (refcount bump —
/// no per-stage allocation).
pub struct BlockShape {
    pub(super) tasks: Arc<[(u32, u32)]>,
    pub(super) packs: Arc<[(usize, usize)]>,
}

impl BlockShape {
    /// The block shape over `blocks` disjoint-range tasks.
    pub fn new(blocks: usize) -> Self {
        Self {
            tasks: (0..blocks as u32).map(|i| (i, 0)).collect(),
            packs: (0..blocks).map(|i| (i, i + 1)).collect(),
        }
    }
}

impl<'a, U> ClassBatchStage<'a, U> {
    /// A block stage over `shape`'s disjoint-range tasks: pack `p` runs
    /// `task(p, fold)` exactly once.
    pub fn block_stage(
        shape: &BlockShape,
        task: &'a (dyn Fn(usize, usize) + Send + Sync),
        fold: usize,
    ) -> Self {
        Self {
            tasks: Arc::clone(&shape.tasks),
            packs: Arc::clone(&shape.packs),
            gateways: Vec::new(),
            jobs: &[],
            block: Some(task),
            fold,
            cohort: false,
        }
    }

    /// The cohort (SoA) sweep variant of a planned sweep stage: identical
    /// tasks, packs, gateways and job views (all lane-invariant — see the
    /// `cohort` field), dispatched to the 4-lane SoA kernel instead of the
    /// scalar sweep. The stage's jobs' operand/result pointers must be
    /// SoA-interleaved with lane stride 4 and the plan must have been
    /// built with `cohort = true` — see `plan_class_sweep_stages`.
    pub fn cohort_variant(stage: &Self) -> Self {
        assert!(
            stage.cohort,
            "cohort_variant requires a plan built with cohort = true              (SoA-interleaved buffers)"
        );
        Self {
            tasks: Arc::clone(&stage.tasks),
            packs: Arc::clone(&stage.packs),
            gateways: stage.gateways.clone(),
            jobs: stage.jobs,
            block: None,
            fold: 0,
            cohort: true,
        }
    }
}


/// Plans the sweep stages for the DAG's dependency levels: per-job gating
/// from the (fixed) support lists, then ONE single-phase per-unit sweep
/// stage per level (tasks = units, balanced pack cuts weighted by unit
/// pair counts). `cohort` builds the 4-lane SoA variant. The stages
/// borrow `levels`; a batched fold plans once and reuses the same stages
/// for every fold in the batch.
///
/// `want_scatter_sets` gates the exact per-job batch scatter sets: the
/// batch path sizes its compact buffers from them, but the per-fold path
/// discards them — computing them there is pure overhead (a `seen` bitvec
/// plus a full entry pass per job, O(jobs × (space + entries)) per level) —
/// so it passes `false` and receives an empty vec.
pub fn plan_class_sweep_stages<'a, T, U>(
    a_series: &LieSeries<T, U>,
    order: &ClassOrder,
    levels: &'a [Vec<KernelJob<'a, U>>],
    cache: &mut GatingCache,
    want_scatter_sets: bool,
    cohort: bool,
) -> (Vec<ClassBatchStage<'a, U>>, Vec<Vec<Vec<u32>>>)
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
    let threads = rayon::current_num_threads().max(1);
    let mut stages = Vec::with_capacity(levels.len());
    // The jobs' TRUE batch scatter sets: the class positions each node's
    // sweep writes under its (fixed) operand supports. The recorded node
    // lists are only a union-level fixed point (the batch eligibility
    // checks the union, not each list) — a level-0 job's list, recorded
    // under an earlier, smaller accumulator support, can miss positions
    // the batch's first fold scatters. Compact buffers must therefore be
    // sized from this exact set. Under the per-unit division the set is
    // exact and free: the gating's `unit_words` list IS the write set
    // (ascending — units in degree order, word sets ascending per unit).
    let mut scatter_sets: Vec<Vec<Vec<u32>>> =
        Vec::with_capacity(if want_scatter_sets { levels.len() } else { 0 });
    // Pass 1: gate every job up front (the support lists are fixed for the
    // whole fold/batch, so no stage's gating depends on kernel results)
    // and total the planned sweep work — the fold-level slot budget that
    // every stage's pack cut below shares. A stage pack exists to give a
    // walk slot work during that stage, so packs must never outnumber the
    // slots the policy will actually run: with the budget derived from the
    // whole fold's work, a tiny fold (12x2: ~66 entries) plans ONE pack per
    // stage and the serial walk drains a single claim instead of nine.
    let level_gateways: Vec<Vec<KernelGating>> = levels
        .iter()
        .map(|level| {
            level
                .iter()
                .map(|j| {
                    LieSeries::<T, U>::kernel_prologue_cached_class(
                        a_series,
                        j.a_nonzero,
                        j.b_nonzero,
                        order,
                        cache,
                    )
                })
                .collect()
        })
        .collect();
    // Planned sweep work per fold: one term per active entry PLUS one
    // multiply-add per contribution pair (the single-phase unit sweep's
    // two work units — the same magnitudes the two-phase planner funded).
    let sweep_entries: usize = level_gateways
        .iter()
        .map(|gateways| {
            gateways
                .iter()
                .map(|g| g.total_pairs + g.tickets.len())
                .sum::<usize>()
        })
        .sum();
    let slots = work_adaptive_slots(threads, usize::MAX, sweep_entries);
    for (level, gateways) in levels.iter().zip(level_gateways) {
        if want_scatter_sets {
            // Exact scatter sets for this level's jobs: the active word
            // positions themselves (each is written exactly once by the
            // sweep, in ascending order).
            let level_sets: Vec<Vec<u32>> = gateways
                .iter()
                .map(|gateway| gateway.unit_words.to_vec())
                .collect();
            scatter_sets.push(level_sets);
        }
        // Flatten (job, unit) tasks: the per-unit division's atomic
        // scheduling units. Units are emitted per job in table (degree)
        // order and packs cut only between units — race-free (unit word
        // sets partition the write slots) and bit-exact (a word's adds
        // stay inside its one unit, in table-entry order).
        let total: usize = gateways.iter().map(|g| g.total_pairs).sum();
        let unit_count: usize = gateways.iter().map(|g| g.units.len()).sum();
        let mut tasks: Vec<(u32, u32)> = Vec::with_capacity(unit_count);
        for (ji, gating) in gateways.iter().enumerate() {
            for ui in 0..gating.units.len() {
                tasks.push((ji as u32, ui as u32));
            }
        }
        // Balanced prefix cuts at unit boundaries: pack `p` of `p_count`
        // takes tasks while `cur < total * (p + 1) / p_count`. Packs are
        // capped by the task count, not the job count: a narrow upper
        // level's few big jobs must still spread across every slot (unit
        // word sets are conflict-free write regions, and each word's
        // accumulation sequence stays whole inside whichever pack owns its
        // unit).
        //
        // The balance weight is the unit's pair count (its WORK proxy: one
        // term computation per active entry plus one multiply-add per
        // pair).
        // Smallest work worth splitting a stage's sweep across an extra
        // pack: below ~8 contributions per pack the per-pack claim/gate
        // cost dominates the pack's compute.
        const MIN_PACK_WORK: usize = 8;
        let p_count = slots
            .min(unit_count.max(1))
            // Work-scaled cap: a stage too small to give every slot a
            // pack of at least `MIN_PACK_WORK` contributions must not
            // spawn micro-packs — the claim/gate churn outweighs the
            // parallelism on narrow grids (e.g. 12x2, whose whole fold is
            // ~66 entries per stage). Balanced bounds below keep the
            // resulting packs even. The leading `slots` cap is the
            // fold-level work budget (see pass 1): packs beyond the
            // policy's slot count would only be drained serially.
            .min(total.div_ceil(MIN_PACK_WORK))
            .max(1);
        let mut packs: Vec<(usize, usize)> = Vec::with_capacity(p_count);
        if p_count > 1 && !tasks.is_empty() {
            let mut start = 0usize;
            let mut cur = 0usize;
            for (ti, &(ji, ui)) in tasks.iter().enumerate() {
                let unit = &gateways[ji as usize].units[ui as usize];
                cur += unit.pairs as usize;
                let next_bound = total * (packs.len() + 1) / p_count;
                if cur >= next_bound && ti + 1 < tasks.len() {
                    packs.push((start, ti + 1));
                    start = ti + 1;
                }
            }
            packs.push((start, tasks.len()));
        } else if !tasks.is_empty() {
            packs.push((0, tasks.len()));
        }
        stages.push(ClassBatchStage {
            tasks: tasks.into(),
            packs: packs.into(),
            gateways,
            jobs: level.as_slice(),
            block: None,
            fold: 0,
            cohort,
        });
    }
    (stages, scatter_sets)
}
