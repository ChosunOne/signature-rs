use super::*;

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
