use super::*;

/// A park/wake gate for one stage boundary. Waiters poll the stage
/// counter briefly (stage work is microsecond-scale; polling costs the
/// finisher nothing), and only a waiter whose gate stays closed parks —
/// registering itself so the finisher knows a kernel wake is warranted.
/// A finisher whose gate has no parked waiters pays a plain store.
/// (`done` remains the counter — it is what waiters check — the word only
/// carries the wake signal; `futex_wait` re-validates the word's value
/// atomically under the kernel lock, and the bounded park timeout
/// self-heals the benign race where a waiter registers just after the
/// finisher looked.)
struct FutexGate {
    word: AtomicU32,
    parked: AtomicUsize,
}

impl FutexGate {
    const fn new() -> Self {
        FutexGate {
            word: AtomicU32::new(0),
            parked: AtomicUsize::new(0),
        }
    }

    /// Sleeps while this gate's word is still `0`, bounded — a signal that
    /// raced past the wake check self-heals on timeout.
    fn park(&self) {
        self.parked.fetch_add(1, Ordering::AcqRel);
        loop {
            if self.word.load(Ordering::Acquire) != 0 {
                break;
            }
            #[cfg(target_os = "linux")]
            unsafe {
                let timeout = libc::timespec {
                    tv_sec: 0,
                    tv_nsec: 200_000, // 200 µs
                };
                let _ = libc::syscall(
                    libc::SYS_futex,
                    self.word.as_ptr(),
                    libc::FUTEX_WAIT | libc::FUTEX_PRIVATE_FLAG,
                    0u32,
                    &timeout,
                    0u32,
                    0u32,
                );
            }
            #[cfg(not(target_os = "linux"))]
            std::thread::sleep(std::time::Duration::from_micros(200));
        }
        self.parked.fetch_sub(1, Ordering::AcqRel);
    }

    /// Signals completion; wakes parked waiters if there are any.
    fn finish(&self) {
        self.word.store(1, Ordering::Release);
        if self.parked.load(Ordering::Acquire) > 0 {
            #[cfg(target_os = "linux")]
            unsafe {
                let _ = libc::syscall(
                    libc::SYS_futex,
                    self.word.as_ptr(),
                    libc::FUTEX_WAKE | libc::FUTEX_PRIVATE_FLAG,
                    i32::MAX,
                );
            }
        }
    }
}

/// Publishes one completed pack of `target`; the finisher that completes
/// the stage signals its gate. Returns true iff this call was the last
/// publisher.
#[inline]
fn finish_pack(done: &AtomicUsize, gate: &FutexGate, target: usize) -> bool {
    // AcqRel: the read half acquires the earlier finishers' releases (so
    // the last finisher's completion observes the completed stage), and
    // the write half releases this pack's writes to the next stage's readers.
    done.fetch_add(1, Ordering::AcqRel) + 1 == target && {
        gate.finish();
        true
    }
}

/// The batch walk's read-only view: the stage chain plus the sweep's
/// class-ordered tables. Operand slices ride inside the stages' jobs.
struct FoldWalk<'a, U> {
    stages: &'a [&'a ClassBatchStage<'a, U>],
    entries: &'a [Entry],
    /// Class position → Lyndon degree (the compact-buffer shifts index by
    /// degree; a contribution's entry positions are class positions).
    degrees: &'a [u8],
    coeffs: &'a [U],
    /// The class-space decomposition rel slice (rows are addressed by
    /// `(decomp_start, decomp_start+1)` spans into this slice; the
    /// coefficients share its indexing).
    decomp_rels: &'a [u32],
}

impl<'a, U: Clone + Neg<Output = U> + Mul<Output = U> + AddAssign + std::hash::Hash + 'static>
    FoldWalk<'a, U>
{
    /// Claims the next unclaimed pack of stage `s` off the stage's atomic
    /// cursor and runs it (a sweep pack or one block); `false` means the
    /// stage is fully claimed and its remaining work is in flight on other
    /// slots. Claims go to whoever asks, so a pack whose original runner is
    /// queued behind pool contention is run by any running slot — a waiter
    /// never waits on a runnable pack. (`Relaxed` suffices: claims only
    /// hand out disjoint pack indices; the produced data itself is
    /// published by the caller's `done` increment with `Release`.)
    #[inline]
    fn claim_and_run_pack(&self, s: usize, claim: &AtomicUsize) -> bool {
        let stage = self.stages[s];
        // Read-only exhaustion peek first: the gate waits drain this cursor
        // in a loop, and a *failed* claim must not RMW the line — 32 slots
        // spin-claiming an exhausted cursor would otherwise keep the cache
        // line exclusive-bouncing for the whole gate wait. Only a claim
        // that looks available touches the cursor.
        if claim.load(Ordering::Relaxed) >= stage.packs.len() {
            return false;
        }
        let pack = claim.fetch_add(1, Ordering::Relaxed);
        if pack >= stage.packs.len() {
            return false;
        }
        match stage.block {
            // Block stage: pack `p` runs block `p` (one task per pack).
            Some(task) => task(stage.tasks[pack].0 as usize, stage.fold),
            None if stage.cohort => self.cohort_sweep_pack_range(stage, pack),
            None => self.sweep_pack_range(stage, pack),
        }
        true
    }

    /// Sweeps pack `pack` of a sweep stage: single-phase per-unit
    /// direct-add. Packs are single-writer regions — a unit's word set
    /// partitions the write slots (ownership invariant) and packs never
    /// split a unit — so the per-word accumulation order is exactly the
    /// serial sweep's regardless of which slot claims which pack.
    #[inline]
    fn sweep_pack_range(&self, stage: &ClassBatchStage<'a, U>, pack: usize) {
        let (start, end) = stage.packs[pack];
        let jobs = stage.jobs;
        let mut t = start;
        while t < end {
            let ji = stage.tasks[t].0 as usize;
            let mut e = t + 1;
            while e < end && stage.tasks[e].0 == stage.tasks[t].0 {
                e += 1;
            }
            let job = &jobs[ji];
            // SAFETY: within a stage the jobs' result buffers are pairwise
            // disjoint (caller contract), and within one job each unit's
            // word set is a single-writer region (two units never write the
            // same word, packs never split a unit), so concurrent packs of
            // the same job read-modify-write disjoint words through raw
            // pointers — no `&mut` is created, because two packs of one job
            // would otherwise hold aliasing `&mut` slices (UB under the
            // aliasing model even with disjoint ranges). Cross-stage
            // accesses are ordered by the stage counters. The buffer
            // outlives the walk (stable allocations, no resize during the
            // dispatch).
            let result: *mut U = job.result;
            let gating = &stage.gateways[ji];
            let a_shift = job.a_shift;
            let b_shift = job.b_shift;
            let r_shift = job.r_shift;
            for &(ji_task, ui) in &stage.tasks[t..e] {
                let _ = ji_task;
                let unit = &gating.units[ui as usize];
                let td = unit.td as usize;
                // The unit's degree-`td` slice starts at `rs` in the
                // working layout; the job's compact buffer stores the
                // target slice at `rs - r_shift[td]`.
                // SAFETY: the shift tables have >= max_degree+1 entries
                // and `td` is a degree of the table (see the plan).
                let base = unit.rs as usize - unsafe { *r_shift.add(td) } as usize;
                for tp in unit.ticket_start as usize..unit.ticket_end as usize {
                    let ticket = gating.tickets[tp];
                    let e_idx = (ticket & TICKET_INDEX_MASK) as usize;
                    let entry = self.entries[e_idx];
                    let p_active = ticket & TICKET_P_ACTIVE != 0;
                    let q_active = ticket & TICKET_Q_ACTIVE != 0;
                    let (i, j) = (entry.i as usize, entry.j as usize);
                    // Compact-layout operand shifts (see `KernelJob`): the
                    // contribution's operand degrees come from the
                    // class-position degree table.
                    let p = self.degrees[i] as usize;
                    let q = self.degrees[j] as usize;
                    // SAFETY: the shift tables have >= max_degree+1 entries
                    // and p/q are degrees of the table; i and j are class
                    // positions in the operands' supports (the tickets
                    // resolved presence). `raw_mul` skips the float
                    // wrappers' per-op NaN checks (raw-float fast path).
                    let term = unsafe {
                        if p_active {
                            let a_sh_p = *a_shift.add(p) as usize;
                            let b_sh_q = *b_shift.add(q) as usize;
                            let mut t =
                                raw_mul(&*job.a.add(i - a_sh_p), &*job.b.add(j - b_sh_q));
                            if q_active {
                                let a_sh_q = *a_shift.add(q) as usize;
                                let b_sh_p = *b_shift.add(p) as usize;
                                raw_add_assign(
                                    &mut t,
                                    &-raw_mul(&*job.a.add(j - a_sh_q), &*job.b.add(i - b_sh_p)),
                                );
                            }
                            t
                        } else {
                            let a_sh_q = *a_shift.add(q) as usize;
                            let b_sh_p = *b_shift.add(p) as usize;
                            -raw_mul(&*job.a.add(j - a_sh_q), &*job.b.add(i - b_sh_p))
                        }
                    };
                    let from = entry.decomp_start as usize;
                    let to = self.entries[e_idx + 1].decomp_start as usize;
                    for (k, &rel) in self.decomp_rels[from..to].iter().enumerate() {
                        // SAFETY: single-writer RMW (unit word-set
                        // ownership, packs never split a unit), in bounds
                        // (the gating addresses the job's compact result
                        // space; `rs + rel - r_shift[td]` is the compact
                        // slot of the row's target word).
                        unsafe {
                            raw_add_assign(
                                &mut *result.add(base + rel as usize),
                                &raw_mul(&self.coeffs[from + k], &term),
                            );
                        }
                    }
                }
            }
            t = e;
        }
    }

    /// The cohort (4-lane SoA) sweep for one pack — dispatched to the
    /// lane-typed kernel in [`cohort`] (the stage's `cohort` flag routes
    /// here; the engine that built the stage verified the cohort
    /// capability for `U`). The stage's jobs carry SoA-interleaved
    /// operand/result/term pointers; every other plan input (tasks, packs,
    /// gateways, support lists) is lane-invariant — see the `cohort`
    /// module's docs for the layout and the bit-exactness argument.
    #[inline]
    fn cohort_sweep_pack_range(&self, stage: &ClassBatchStage<'a, U>, pack: usize) {
        cohort::sweep_pack_range(
            self.entries,
            self.degrees,
            self.coeffs,
            self.decomp_rels,
            stage,
            pack,
        );
    }
}

/// The slot policy's work quantum, in planned active-entry tickets per
/// fold unit per slot: one more coordinated slot must be funded by this
/// much real sweep work. Calibration: a swept entry costs ~40 cycles
/// (measured marginal on the class kernel) and a coordinated slot pays
/// the per-stage claim/gate protocol (~2-5 µs across a stage boundary),
/// so 3750 entries × 40 cy ≈ 150 Kcy ≈ 50 µs of sweep work per slot at
/// ~3 GHz — the break-even where an extra slot stops paying for the
/// barriers it joins. Tuned between the two candidates (1500 ≈ 20 µs and
/// 3750 ≈ 50 µs equivalents): both fix the tiny regimes, 3750 also wins
/// the one regime they disagree on (3x8 at a 32t pool: e2e 161 ms vs
/// 171 ms). Measured regimes: 12x2 (66 entries/fold) and 8x3 (~700/fold)
/// walk serially at every pool size (their e2e used to run 6.5x slower at
/// 16t than 1t — barrier machinery dwarfing the work); 3x8 (~74K/fold)
/// lands near its measured 8-16t sweet spot (19 slots at a 32t pool);
/// 2x12 (~247K/fold) still fills the pool.
const SLOT_WORK_QUANTUM: usize = 3750;

/// SLOT_POLICY_DEBUG=1: one stderr line per walk dispatch with the
/// policy's inputs and decision (tuning/ops diagnostics; off by default).
fn slot_policy_debug() -> bool {
    static DEBUG: OnceLock<bool> = OnceLock::new();
    *DEBUG.get_or_init(|| std::env::var_os("SLOT_POLICY_DEBUG").is_some())
}

/// The stage chain's planned SWEEP work in swept-entry units: sweep
/// stages contribute their gateways' precomputed active-entry tickets (the
/// exact per-entry work unit the sweeps iterate — presence-resolved, so
/// canceled entries don't inflate the estimate); block stages contribute
/// nothing (their element count depends on the caller's block cut, which
/// is itself chosen from this estimate). Pure plan-time data: no kernel
/// results are needed.
pub fn planned_sweep_entries<U>(stages: &[&ClassBatchStage<'_, U>]) -> usize {
    stages
        .iter()
        .map(|s| match s.block {
            Some(_) => 0,
            None => s
                .gateways
                .iter()
                .map(|g| g.total_pairs + g.tickets.len())
                .sum::<usize>(),
        })
        .sum()
}

/// The stage chain's planned work in swept-entry units: sweep stages
/// contribute their gateways' precomputed active-entry tickets, block
/// stages their element count (one block = one disjoint class-position
/// range per pack). Pure plan-time data: no kernel results are needed.
fn planned_stage_entries<U>(stages: &[&ClassBatchStage<'_, U>]) -> usize {
    stages
        .iter()
        .map(|s| match s.block {
            Some(_) => s.packs.len(),
            None => s
                .gateways
                .iter()
                .map(|g| g.total_pairs + g.tickets.len())
                .sum::<usize>(),
        })
        .sum()
}

/// Work-adaptive slot count for the stage-chain walk: the walk's parallel
/// width is chosen from the PLANNED work, not just the pool size. Every
/// slot joins every stage boundary's counter/gate protocol, so a slot is
/// worth spawning only when its share of the per-fold-unit work covers the
/// coordination cost: `per_unit_work / QUANTUM` slots, capped by the pool
/// and by the widest stage's pack count (a slot beyond every stage's packs
/// only adds barrier arrivals, never runs work). Callers that cut their
/// own stage structure (e.g. a batch's block count) should derive it from
/// the same policy so packs never outnumber useful slots.
pub fn work_adaptive_slots(threads: usize, max_packs: usize, per_unit_work: usize) -> usize {
    threads
        .min(max_packs)
        .min((per_unit_work / SLOT_WORK_QUANTUM).max(1))
        .max(1)
}

/// Runs a linear chain of fold stages as ONE parallel dispatch (a serial
/// loop on the calling thread whenever the planned work cannot fund more
/// than one coordinated slot). Every slot walks the stages in order:
/// before waiting on a stage's counter it drains that stage's claim cursor
/// — a queued stage's packs are run by whoever is working, never waited
/// for — so the walk stays live under any pool contention. The per-stage
/// counters (Release on completion, Acquire at the boundary) are the only
/// synchronization: a stage's reads observe exactly the writes of all
/// earlier stages, and the per-word accumulation order is exactly the
/// serial schedule's (packs never split a unit; blocks never split a
/// range).
///
/// `fold_units` is how many per-fold repetitions the stage chain contains
/// (1 for the per-fold path; a batch of `rhss.len()` folds repeats its
/// gather/sweep/accumulate sub-chain once per displacement) — the slot
/// policy normalizes the planned work by it, because the barrier cost is
/// paid per stage, i.e. per fold unit, not once per dispatch.
pub fn run_class_batch<T, U>(
    a_series: &LieSeries<T, U>,
    order: &ClassOrder,
    stages: &[&ClassBatchStage<'_, U>],
    fold_units: usize,
) where
    T: Clone + Ord + Generator + Hash,
    U: Clone
        + Neg<Output = U>
        + Mul<Output = U>
        + AddAssign
        + std::hash::Hash
        + Send
        + Sync
        + 'static,
{
    run_class_batch_with_work(
        a_series,
        order,
        stages,
        fold_units,
        planned_stage_entries(stages),
    );
}

/// [`run_class_batch`] with the caller's planned-work estimate: a cohort
/// chain's sweep stages carry the 4-lane work that the plain stage-based
/// estimate (entry tickets, lane-agnostic plan data) cannot see, so the
/// cohort engine supplies `planned_work = scalar_tickets × lanes` and the
/// slot policy funds the walk from the true per-unit work.
pub fn run_class_batch_with_work<T, U>(
    a_series: &LieSeries<T, U>,
    order: &ClassOrder,
    stages: &[&ClassBatchStage<'_, U>],
    fold_units: usize,
    planned_work: usize,
) where
    T: Clone + Ord + Generator + Hash,
    U: Clone
        + Neg<Output = U>
        + Mul<Output = U>
        + AddAssign
        + std::hash::Hash
        + Send
        + Sync
        + 'static,
{
    // The walk is fully internal for sweep stages: relabeled entries gate
    // the presence tests and index the class-ordered operands. Results are
    // written per job through the job's own buffer (see `sweep_pack_range`).
    let walk = FoldWalk {
        stages,
        entries: order.entries_cls(),
        degrees: order.degree_cls(),
        coeffs: a_series.feasible_decompositions.decomp_coeffs(),
        decomp_rels: order.decomp_cls(),
    };

    // Pack-slot walk with dynamic claims: `slots` tasks walk every stage,
    // and each stage's packs are handed out by an atomic cursor at run
    // time. A slot therefore never waits on a *queued* runner: before it
    // blocks on a stage counter it drains that stage's claim cursor itself,
    // so every unclaimed pack is either claimed-and-in-flight (a running
    // task publishes it via the counter) or claimable by the waiting slot.
    // That keeps the walk live even when the pool is shared with outer
    // parallelism and slot tasks queue behind other work — the fixed
    // slot→pack assignment this replaces starved exactly there.
    let stage_pack_counts: Vec<usize> = stages.iter().map(|s| s.packs.len()).collect();
    // One padded cell per counter: adjacent stages' counters must not share
    // a cache line — 32 slots RMW-ing two hot counters per gate otherwise
    // drag every *other* stage's counters along with them.
    #[repr(align(128))]
    struct Pad(AtomicUsize);
    let done: Vec<Pad> = stage_pack_counts
        .iter()
        .map(|_| Pad(AtomicUsize::new(0)))
        .collect();
    let claims: Vec<Pad> = stage_pack_counts
        .iter()
        .map(|_| Pad(AtomicUsize::new(0)))
        .collect();
    let gates: Vec<FutexGate> = stage_pack_counts.iter().map(|_| FutexGate::new()).collect();
    let threads = rayon::current_num_threads().max(1);
    let max_packs = stage_pack_counts.iter().copied().max().unwrap_or(1);
    let planned_entries = planned_work;
    let per_unit_work = planned_entries / fold_units.max(1);
    let slots = work_adaptive_slots(threads, max_packs, per_unit_work);
    if slot_policy_debug() {
        eprintln!(
            "slot_policy: fold_units={fold_units} planned_entries={planned_entries} \
             per_unit={per_unit_work} threads={threads} max_packs={max_packs} slots={slots}"
        );
    }
    let walk_for_slot = |_slot: usize| {
        // Hoisted per-stage sync references: the drain/spin loops below run
        // per claim attempt / per spin iteration, and re-indexing the
        // counter Vecs there cost ~7% of the wide-pool profile in bounds
        // checks + deref chains alone (flamegraph, 4x8/32t). The Vecs are
        // never resized during the walk, so shared references hoisted once
        // per slot are equivalent.
        let stage_done: Vec<&AtomicUsize> = done.iter().map(|p| &p.0).collect();
        let stage_claims: Vec<&AtomicUsize> = claims.iter().map(|p| &p.0).collect();
        let stage_gates: Vec<&FutexGate> = gates.iter().collect();
        for s in 0..stages.len() {
            if s > 0 {
                // EVERY slot waits at every stage boundary — including
                // slots that claimed nothing at the previous stage — so no
                // task of stage `s` can read a buffer written by stage
                // `s - 1` before all of stage `s - 1`'s packs have released
                // (Release in the run callers, Acquire here). While
                // waiting, keep draining the claim cursor: if the slots
                // that claimed the stage's remaining packs are queued
                // behind pool contention, this running slot runs those
                // packs itself instead of blocking on them. Spin briefly,
                // then yield to co-scheduled rayon work, then sleep — at
                // this point every remaining pack is in flight, so the
                // wait is bounded by real work.
                let need = stage_pack_counts[s - 1];
                let (prev_done, prev_gate, prev_claim) =
                    (stage_done[s - 1], stage_gates[s - 1], stage_claims[s - 1]);
                // Drain once on arrival: after this the cursor is exhausted
                // — every pack is claimed-and-in-flight — so parking is
                // safe (no pack waits on this slot). The cursor never
                // grows back, so no re-drain is ever needed.
                while walk.claim_and_run_pack(s - 1, prev_claim) {
                    finish_pack(prev_done, prev_gate, need);
                }
                // Poll briefly (most gates open within microseconds —
                // polling costs the finisher nothing), then YIELD instead
                // of parking: a parked slot BLOCKS its rayon worker, and
                // the tournament runs many concurrent dispatches on one
                // pool (e2e: 32 leaf-group tasks × slot tasks on 32
                // workers) — queued sibling slot tasks sit behind the
                // parked holder, which measured as 40% slot utilization
                // at 32 threads (task-clock vs wall, 4x8 concat). A
                // yielding slot re-queues its task so the worker runs the
                // other dispatches' pack work; the drain-on-arrival rule
                // guarantees every pack is claimed-and-in-flight, so the
                // gate always opens from real work. Parking stays as the
                // long-wait fallback (frees the core when a stage is
                // genuinely slow; the finisher then pays the wake).
                let mut spins = 0usize;
                let mut yields = 0usize;
                while prev_done.load(Ordering::Acquire) < need {
                    if spins < 1 << 12 {
                        std::hint::spin_loop();
                        spins += 1;
                    } else if yields < 1 << 10 {
                        yields += 1;
                        rayon::yield_now();
                    } else {
                        prev_gate.park();
                    }
                }
            }
            // Run this stage's packs off the claim cursor and publish each
            // completion: `done[s]` reaches exactly `packs[s]` when — and
            // only when — every pack has been run by whichever slot
            // claimed it. (A slot that finds the cursor drained reports
            // nothing, or the stage would appear complete before its work
            // is done.)
            let (this_done, this_gate, this_claim) =
                (stage_done[s], stage_gates[s], stage_claims[s]);
            let this_need = stage_pack_counts[s];
            while walk.claim_and_run_pack(s, this_claim) {
                finish_pack(this_done, this_gate, this_need);
            }
        }
    };

    // The serial walk is entered by the POLICY's slot count, not the pool
    // size: `slots == 1` runs the identical walk on the calling thread —
    // the claims/counters/gates are self-contained per slot, so the
    // single-slot semantics are exactly the single-worker pool's long-
    // standing serial path, minus the par_iter dispatch cost.
    if slots <= 1 {
        walk_for_slot(0);
    } else {
        use rayon::prelude::*;
        (0..slots).into_par_iter().for_each(walk_for_slot);
    }
}

/// Multi-level batch kernel in the class-contiguous working mode: the whole
/// fold runs as ONE parallel dispatch.
///
/// `levels[l]` holds the jobs of dependency level `l`; every job's operands
/// are class-ordered slices of buffers produced by strictly earlier levels
/// (or the fold's inputs), its support lists are class-indexed, and its
/// result buffer receives the class-ordered sum directly. Jobs of one level
/// must have pairwise disjoint result buffers.
///
/// The call plans every level up front (gating from the fixed support
/// lists, then balanced class-contiguous packs) and runs one parallel walk
/// of sweep stages. Each stage's packs are claimed dynamically off an
/// atomic cursor, every completion is published with `Release`, and every
/// slot waits (`Acquire`) at every stage boundary — so cross-level operand
/// reads are ordered by the stage counters, the per-word accumulation
/// order is exactly the serial sweep's (units are never split across
/// packs), and a slot never blocks on runnable fold work: before waiting
/// it drains the level's claim cursor itself.
pub fn commutator_coefficients_class_fold_with_cache<T, U>(
    a_series: &LieSeries<T, U>,
    order: &ClassOrder,
    levels: &mut [Vec<KernelJob<'_, U>>],
    cache: &mut GatingCache,
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
        + Mul<Output = U>
        + AddAssign
        + Send
        + Sync
        + 'static,
{
    debug_assert_eq!(
        order.inv().len(),
        a_series.basis.len(),
        "class ordering does not describe this series' basis"
    );

    // The per-fold path discards the scatter sets (only the batch path
    // sizes compact buffers from them): skip their exact-set computation —
    // the gating and pack cuts the sweep actually reads are unchanged.
    let (stages, _scatter_sets) =
        plan_class_sweep_stages(a_series, order, levels, cache, false, false);
    let stage_refs: Vec<&ClassBatchStage<U>> = stages.iter().collect();
    // One fold unit: the stage chain is a single fold's sweep stages, so
    // the slot policy sees the fold's own planned work.
    run_class_batch(a_series, order, &stage_refs, 1);
}

#[cfg(test)]
mod tests {
    /// The work-adaptive slot policy must reproduce the measured regimes:
    /// tiny folds walk serially at any pool size, mid folds land on their
    /// measured sweet spot, wide folds fill the pool, and the pool/pack
    /// caps still bind.
    #[test]
    fn work_adaptive_slots_matches_measured_regimes() {
        use super::work_adaptive_slots;
        // 12x2 per-fold (66 entries/fold): serial at every pool size.
        assert_eq!(work_adaptive_slots(32, 32, 66), 1);
        assert_eq!(work_adaptive_slots(2, 32, 66), 1);
        // 8x3 (~700 entries/fold, + block elements in the batch path):
        // 1-2 slots.
        assert_eq!(work_adaptive_slots(32, 32, 700), 1);
        assert!(work_adaptive_slots(32, 32, 2000) <= 2, "8x3 regime: 2000 entries must fund at most 2 slots");
        // 3x8 (~74K planned entries/fold): the QUANTUM lands it near the
        // measured 8-16t sweet spot — 19 slots at a 32t pool (74_000/3750).
        assert_eq!(work_adaptive_slots(32, 32, 74_000), 19);
        // ... and stays in that neighborhood across the regime's ticket
        // variance (the real 3x8 fold's planned count wobbles around 74K).
        let s = work_adaptive_slots(32, 32, 70_000);
        assert!((16..=21).contains(&s), "3x8 regime: 70K entries must land in 16..=21 slots, got {s}");
        // 2x12 (~257K entries/fold): the full pool.
        assert_eq!(work_adaptive_slots(32, 32, 256_858), 32);
        // Pool and pack caps bind before the work term.
        assert_eq!(work_adaptive_slots(8, 3, 256_858), 3);
        assert_eq!(work_adaptive_slots(1, 32, 256_858), 1);
        assert_eq!(work_adaptive_slots(32, 1, 256_858), 1);
        // Degenerate: no planned work → serial.
        assert_eq!(work_adaptive_slots(32, 32, 0), 1);
    }
}
