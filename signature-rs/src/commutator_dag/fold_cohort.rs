//! The 4-lane SIMD-across-folds cohort engine: its lane type, the kill
//! switch, the observability counters and the SoA-interleaved batch fold.

use lie_rs::{
    BlockShape, ClassBatchStage, ClassOrderedCommutation, KernelJob, LieSeries,
    plan_class_sweep_stages, planned_sweep_entries, run_class_batch_with_work,
    work_adaptive_slots,
};
use lyndon_rs::generators::Generator;
use num_traits::{One, Zero};
use std::cell::UnsafeCell;
use std::hash::Hash;
use std::ops::{AddAssign, MulAssign, Neg};
use std::sync::atomic::{AtomicBool, AtomicUsize};
use std::sync::{Arc, OnceLock};

use super::block_work::{
    AccumRun, BlockWork, RhsView, SendBlockWork, SendRaw, SendRhsViews, ZeroRun,
    for_each_position_run,
};
use super::{CommutatorDag, DagNode, SharedNodeLists, TermSource};

/// One lane of a fold cohort: the lane's accumulator (public basis order,
/// updated in place by [`CommutatorDag::fold_batch_cohort`]) and the
/// displacement slices it still has to fold. Lanes of one cohort share the
/// (support-derived, value-independent) fold plan; they differ only in
/// coefficient values — see the `fold_batch_cohort` doc.
pub(crate) struct CohortLane<'a, U> {
    /// The lane's accumulator values, public basis order, full basis
    /// length. In-out: holds the pre-fold values on entry, the folded
    /// result on return.
    pub acc: Vec<U>,
    /// The lane's remaining displacement slices (public basis order).
    pub rhss: &'a [&'a [U]],
}

/// The cohort engine's kill switch. Initialized once from the environment
/// (`SIG_NO_COHORT=1` — the A/B benchmark's off switch and the scalar
/// oracle's handle), then settable in-process: the bit-exactness oracle
/// test flips it to run the scalar tournament against the cohort one
/// (nextest runs each test in its own process, so the process-global flag
/// cannot leak across tests).
static COHORT_OFF: AtomicBool = AtomicBool::new(false);
static COHORT_ENV_INIT: OnceLock<()> = OnceLock::new();

/// Cohort engine executions ([`CommutatorDag::fold_batch_cohort`] calls —
/// leaf-group warm-up/tail cohorts plus merge groups). The oracle test's
/// observability that the cohort path actually engaged for a configuration.
pub(crate) static COHORT_ENGINE_RUNS: AtomicUsize = AtomicUsize::new(0);

/// Individual lane-folds executed through the cohort engine: for every
/// dispatch, the number of (step, active-lane) pairs — the folds that
/// shared one plan walk instead of walking the scalar engine one fold at
/// a time. The engagement test's fold-coverage observable: cohort
/// lane-folds / total folds proves the SIMD-across-folds engine ran
/// every eligible concatenated segment, not only the merge rounds.
pub(crate) static COHORT_LANE_FOLDS: AtomicUsize = AtomicUsize::new(0);

/// True when the COEFFICIENT TYPE and CPU support the cohort engine,
/// irrespective of the kill switch: a repr-transparent raw float type
/// (mirroring the kernel's raw-float fast path — every other type, e.g.
/// `Ratio<i128>`, keeps the scalar kernel bit-for-bit) and AVX2 (the
/// conservative floor for the 4×f64 lane vectors; the kernel itself is
/// plain autovectored array code — no FMA, so results stay bit-identical
/// at any vector width).
///
/// This is the driver's tree-selection predicate: the
/// tournament-vs-sequential choice must NOT depend on the kill switch —
/// the switch swaps the ENGINES inside the tournament (cohort vs scalar),
/// never the association tree, so the oracle tests can pin the two
/// bit-identical.
///
/// All checks are per-call cheap: two `TypeId` compares (constant-folded
/// per monomorphization) and one CPU-feature cache probe.
pub(crate) fn cohort_type_capable<U: 'static>() -> bool {
    if !lie_rs::cohort_supported::<U>() {
        return false;
    }
    #[cfg(target_arch = "x86_64")]
    {
        std::arch::is_x86_feature_detected!("avx2")
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        false
    }
}

/// True when the cohort (4-lane SIMD-across-folds) engine may run for
/// coefficient type `U`: [`cohort_type_capable`] plus the kill switch
/// must not disable the path.
pub(crate) fn cohort_capable<U: 'static>() -> bool {
    COHORT_ENV_INIT.get_or_init(|| {
        if std::env::var_os("SIG_NO_COHORT").is_some_and(|v| v == *"1") {
            COHORT_OFF.store(true, std::sync::atomic::Ordering::Relaxed);
        }
    });
    if COHORT_OFF.load(std::sync::atomic::Ordering::Relaxed) {
        return false;
    }
    cohort_type_capable::<U>()
}

/// Flips the cohort kill switch in-process (the oracle test's scalar
/// handle; the environment init runs first so an explicit flip is never
/// overwritten).
#[cfg(test)]
pub(crate) fn set_cohort_off(off: bool) {
    COHORT_ENV_INIT.get_or_init(|| {});
    COHORT_OFF.store(off, std::sync::atomic::Ordering::Relaxed);
}

impl<U> CommutatorDag<U>
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
    /// Cohort (4-lane SIMD-across-folds) variant of [`Self::fold_batch`]:
    /// runs every lane's remaining folds in lockstep through ONE shared
    /// plan, carrying the lanes' coefficient values as SoA-interleaved
    /// 4-lane vectors (lane `l`'s value of class position `p` lives at
    /// element `p*4 + l` of the working buffers).
    ///
    /// The plan (decomposition lists, scatter indices, node order, BCH
    /// weights, gating) is a pure function of the basis and the atom
    /// supports — data-independent — so one walk of it serves every lane:
    /// the lanes differ only in coefficient VALUES, and every plan input
    /// (tasks, packs, gateways, support lists, shift tables) is
    /// lane-invariant. The lane arithmetic replicates `fold_batch`'s per-
    /// fold operation order exactly (plain mul+add, no FMA — Rust does not
    /// contract), so each lane's result is BIT-IDENTICAL to running that
    /// lane's folds through `fold_batch` on its own.
    ///
    /// # Contract (checked by the callers)
    ///
    /// - `lanes.len()` is in `2..=4` and every lane's `acc` is full basis
    ///   length;
    /// - every lane's accumulator support equals `a_nonzero` and every
    ///   lane's remaining displacements share the support `b_nonzero`
    ///   (for multi-fold lanes the accumulator's support must also equal
    ///   the reachable set — the full `batch_eligible` condition — so no
    ///   lane's gating can be outgrown mid-tail; for one-fold lanes the
    ///   density is irrelevant);
    /// - the node lists are the built fixed point for `(a_nonzero,
    ///   b_nonzero)` — identical across lanes by the plan's data
    ///   independence, and exactly what `lists_steady_for`/
    ///   `ensure_lists_steady` establish on the shared DAG.
    ///
    /// Lanes with fewer remaining folds than the cohort's step count are
    /// masked off per fold step (their gather/accumulate stages skip
    /// them; the shared sweeps' lane values for masked lanes are dead —
    /// never consumed). Masked folds accumulate per-lane scalar (the
    /// vector path requires all four lanes active: a zero-multiplier
    /// mask could turn an overflowing stale lane value into a NaN add,
    /// diverging from the scalar engine's skip).
    pub(crate) fn fold_batch_cohort<T>(
        &mut self,
        series: &LieSeries<T, U>,
        lanes: &mut [CohortLane<'_, U>],
        a_nonzero: &[usize],
        b_nonzero: &[usize],
    ) where
        T: Clone + Ord + Generator + Hash + Eq,
        U: Send + Sync,
    {
        let lane_n = lanes.len();
        debug_assert!((2..=4).contains(&lane_n), "cohort lanes must be in 2..=4");
        debug_assert!({
            let d = series.coefficients.len();
            lanes.iter().all(|l| l.acc.len() == d)
        });
        // The cohort's fold-step count: the longest lane's remaining folds.
        let steps = lanes.iter().map(|l| l.rhss.len()).max().unwrap_or(0);
        if steps == 0 {
            return;
        }
        COHORT_ENGINE_RUNS.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        // Per-step activity: lane l folds step f iff it still has one.
        let active_masks: Vec<u8> = (0..steps)
            .map(|f| {
                (0..lane_n)
                    .filter(|&l| f < lanes[l].rhss.len())
                    .fold(0u8, |m, l| m | (1 << l))
            })
            .collect();
        COHORT_LANE_FOLDS.fetch_add(
            active_masks.iter().map(|m| m.count_ones() as usize).sum(),
            std::sync::atomic::Ordering::Relaxed,
        );

        let internal = self.structure.nodes.len() - 2;
        self.ensure_buffers(series.coefficients.len(), internal);
        let order = self.class_order.get_or_init(|| series.class_order());
        let inv = order.inv();
        let perm = order.perm();
        let d = series.coefficients.len();
        const LANES: usize = lie_rs::COHORT_LANES;

        // The cohort's class-space working state, SoA-interleaved: the
        // accumulators stay class-ordered across steps (step f's sweep reads
        // them as its `a` operand and step f's accumulate adds the terms in
        // place); the displacement buffer is re-gathered per step. Same
        // UnsafeCell + stage-counter protocol as `fold_batch`.
        let mut acc4 = vec![U::default(); LANES * d];
        for (l, lane) in lanes.iter().enumerate() {
            for (w, c) in lane.acc.iter().enumerate().take(d) {
                acc4[inv[w] as usize * LANES + l] = c.clone();
            }
        }
        let acc_cell = UnsafeCell::new(acc4);
        let b_cell = UnsafeCell::new(vec![U::default(); LANES * d]);
        // SAFETY: the cells are not moved for the dispatch's duration, so
        // the heap buffers behind their Vecs are stable; everything below
        // goes through these data pointers (the sweeps read shared slices
        // derived from them, the block stages read and write through
        // them), and the stage counters order every write-before-read.
        let acc_data: *mut U = unsafe { (*acc_cell.get()).as_mut_ptr() };
        let b_data: *mut U = unsafe { (*b_cell.get()).as_mut_ptr() };

        // Degree layout of the class space — verbatim from `fold_batch`.
        let pos_degree: Vec<u8> = order.degree_cls().to_vec();
        let max_deg = pos_degree.iter().copied().max().unwrap_or(0) as usize;
        let mut deg_start: Vec<u32> = vec![d as u32; max_deg + 2];
        for (pos, &dg) in pos_degree.iter().enumerate() {
            let dg = dg as usize;
            if (deg_start[dg] as usize) > pos {
                deg_start[dg] = pos as u32;
            }
        }
        deg_start[max_deg + 1] = d as u32;
        for dd in (0..=max_deg).rev() {
            if deg_start[dd] > deg_start[dd + 1] {
                deg_start[dd] = deg_start[dd + 1];
            }
        }
        let zero_shifts = vec![0u32; max_deg + 2];
        let zero_shifts_ptr: *const u32 = zero_shifts.as_ptr();

        // One job table for the whole cohort: level-0 jobs read the SoA
        // accumulator and displacement buffers with the shared supports;
        // interior jobs read child node buffers with the children's
        // scatter-target lists — exactly the `fold_batch` operand fill, at
        // SoA bases (the SoA sweep multiplies every slot index by `LANES`).
        let a_nz_cls: Vec<usize> = a_nonzero.iter().map(|&i| inv[i] as usize).collect();
        let b_nz_cls: Vec<usize> = b_nonzero.iter().map(|&i| inv[i] as usize).collect();
        let mut levels_jobs: Vec<Vec<KernelJob<U>>> = self
            .structure
            .levels
            .iter()
            .skip(1)
            .map(|level| {
                level
                    .iter()
                    .filter_map(|&k| match self.structure.nodes[k as usize] {
                        DagNode::Binary { .. } => Some(KernelJob {
                            a: acc_data,
                            a_len: LANES * d,
                            a_nonzero: &a_nz_cls,
                            b: b_data,
                            b_len: LANES * d,
                            b_nonzero: &b_nz_cls,
                            result: std::ptr::null_mut(),
                            result_len: 0,
                            a_shift: zero_shifts_ptr,
                            b_shift: zero_shifts_ptr,
                            r_shift: zero_shifts_ptr,
                        }),
                        DagNode::Atom(_) => None,
                    })
                    .collect()
            })
            .collect();
        // SAFETY: the operand slices are derived from the cell data
        // pointers and the node buffers (stable allocations, no resize
        // during the dispatch) and only read (by the sweeps); the block
        // stages write through the same raws between stages — the
        // UnsafeCell makes that interleaving defined, and the stage
        // counters order it. The allocations live to the end of this
        // function (the dispatch joins before they drop).
        let lists_local: SharedNodeLists = Arc::clone(&self.nonzeros);
        for (li, level) in self.structure.levels.iter().enumerate().skip(1) {
            for (jj, &k) in level.iter().enumerate() {
                let (left, right) = match self.structure.nodes[k as usize] {
                    DagNode::Binary { left, right } => (left, right),
                    DagNode::Atom(_) => continue,
                };
                let (a_ptr, a_len, a_sh, a_nz): (*const U, usize, *const u32, &[usize]) = match left
                {
                    0 => (acc_data as *const U, LANES * d, zero_shifts_ptr, &a_nz_cls),
                    1 => (b_data as *const U, LANES * d, zero_shifts_ptr, &b_nz_cls),
                    id => (
                        std::ptr::null(),
                        0,
                        zero_shifts_ptr,
                        &lists_local[id as usize - 2],
                    ),
                };
                let (b_ptr, b_len, b_sh, b_nz): (*const U, usize, *const u32, &[usize]) =
                    match right {
                        0 => (acc_data as *const U, LANES * d, zero_shifts_ptr, &a_nz_cls),
                        1 => (b_data as *const U, LANES * d, zero_shifts_ptr, &b_nz_cls),
                        id => (
                            std::ptr::null(),
                            0,
                            zero_shifts_ptr,
                            &lists_local[id as usize - 2],
                        ),
                    };
                let job = &mut levels_jobs[li - 1][jj];
                job.a = a_ptr;
                job.a_len = a_len;
                job.a_nonzero = a_nz;
                job.b = b_ptr;
                job.b_len = b_len;
                job.b_nonzero = b_nz;
                job.a_shift = a_sh;
                job.b_shift = b_sh;
                job.r_shift = zero_shifts_ptr;
            }
        }

        // Plan once to get the jobs' exact cohort scatter sets, then record
        // them as the node lists (same fixed-point bookkeeping as
        // `fold_batch`).
        let (_, scatter_sets) =
            plan_class_sweep_stages(
                series,
                order,
                &levels_jobs,
                &mut self.gating_cache,
                true,
                true,
            );
        {
            // Copy-on-write re-record (same fixed-point bookkeeping as
            // `fold_batch`): the exact cohort scatter sets replace EVERY
            // binary node's list, so a fresh shell avoids copying the old
            // lists — crucial when they are an adopted cache snapshot
            // shared by `Arc` (a shared snapshot must never be mutated).
            let internal = self.structure.nodes.len() - 2;
            let mut updated: Vec<Vec<usize>> = vec![Vec::new(); internal];
            let mut written = 0_usize;
            for (li, level) in self.structure.levels.iter().enumerate().skip(1) {
                for (jj, &k) in level.iter().enumerate() {
                    if matches!(self.structure.nodes[k as usize], DagNode::Binary { .. }) {
                        updated[k as usize - 2] = scatter_sets[li - 1][jj]
                            .iter()
                            .map(|&p| p as usize)
                            .collect();
                        written += 1;
                    }
                }
            }
            debug_assert_eq!(written, internal, "every internal node must be re-recorded");
            self.nonzeros = Arc::new(updated);
        }

        // Compact per-node scratch buffers, SoA-interleaved (each scalar
        // slot becomes `LANES` contiguous lane elements): sizing, degree
        // slices and shift tables are lane-agnostic — verbatim from
        // `fold_batch`, at 4x element counts.
        struct CompactSlice {
            base: u32,
            rs: u32,
            re: u32,
        }
        let mut compact: Vec<Vec<U>> = Vec::with_capacity(internal);
        let mut shifts: Vec<Vec<u32>> = Vec::with_capacity(internal);
        let mut layouts: Vec<Vec<CompactSlice>> = Vec::with_capacity(internal);
        for _ in 0..internal {
            compact.push(Vec::new());
            shifts.push(vec![0u32; max_deg + 2]);
            layouts.push(Vec::new());
        }
        for (li, level) in self.structure.levels.iter().enumerate().skip(1) {
            for (jj, &k) in level.iter().enumerate() {
                let set: &[u32] =
                    if matches!(self.structure.nodes[k as usize], DagNode::Binary { .. }) {
                        &scatter_sets[li - 1][jj]
                    } else {
                        &[]
                    };
                let ki = k as usize - 2;
                let mut degs: Vec<u8> = set.iter().map(|&p| pos_degree[p as usize]).collect();
                degs.sort_unstable();
                degs.dedup();
                let mut buf: Vec<U> = Vec::new();
                let shift = &mut shifts[ki];
                let slices = &mut layouts[ki];
                for dg in degs {
                    let rs = deg_start[dg as usize];
                    let re = deg_start[dg as usize + 1];
                    // The shift stays in SCALAR slot units (class position
                    // minus compact slot) — the SoA sweep scales slot
                    // indices by `LANES` itself — so `base` counts scalar
                    // slots: the SoA element count divided by the lane
                    // width (always exact — the buffer grows in
                    // `LANES`-element steps).
                    let base = (buf.len() / LANES) as u32;
                    buf.resize(buf.len() + (re - rs) as usize * LANES, U::default());
                    shift[dg as usize] = rs - base;
                    slices.push(CompactSlice { base, rs, re });
                }
                compact[ki] = buf;
            }
        }
        let compact_ptrs: Vec<(*mut U, usize)> = compact
            .iter_mut()
            .map(|b| (b.as_mut_ptr(), b.len()))
            .collect();

        // Wire the jobs' results and the interior jobs' operand views to the
        // compact SoA buffers (the shift tables stay in class positions —
        // the SoA sweep scales slot indices by `LANES`).
        for (li, level) in self.structure.levels.iter().enumerate().skip(1) {
            for (jj, &k) in level.iter().enumerate() {
                let job = &mut levels_jobs[li - 1][jj];
                let (rp, rl) = compact_ptrs[k as usize - 2];
                job.result = rp;
                job.result_len = rl;
                job.r_shift = shifts[k as usize - 2].as_ptr();
                let (left, right) = match self.structure.nodes[k as usize] {
                    DagNode::Binary { left, right } => (left, right),
                    DagNode::Atom(_) => continue,
                };
                if left > 1 {
                    let (p, l) = compact_ptrs[left as usize - 2];
                    job.a = p as *const U;
                    job.a_len = l;
                    job.a_shift = shifts[left as usize - 2].as_ptr();
                }
                if right > 1 {
                    let (p, l) = compact_ptrs[right as usize - 2];
                    job.b = p as *const U;
                    job.b_len = l;
                    job.b_shift = shifts[right as usize - 2].as_ptr();
                }
            }
        }

        // Final plan with the wired jobs (gating cache hit — supports
        // unchanged), then the SoA sweep stages over it: identical tasks,
        // packs and gateways, dispatched to the 4-lane kernel.
        let (sweep_stages, _) =
            plan_class_sweep_stages(
                series,
                order,
                &levels_jobs,
                &mut self.gating_cache,
                true,
                true,
            );
        let sweep_stages_soa: Vec<ClassBatchStage<U>> = sweep_stages
            .iter()
            .map(ClassBatchStage::cohort_variant)
            .collect();

        // Block ranges over class positions — same work-adaptive policy as
        // `fold_batch`, funded by the COHORT's per-step work (scalar sweep
        // tickets × lanes; see `run_class_batch_with_work`).
        let threads = rayon::current_num_threads().max(1);
        let planned_scalar = planned_sweep_entries(&sweep_stages.iter().collect::<Vec<_>>());
        let slots = work_adaptive_slots(threads, d.max(1), planned_scalar * lane_n);
        let blocks = d.min(4 * slots).max(1);
        // The block partition must match `for_each_position_run`'s block
        // assignment (`p * blocks / d`, floor): block `b` is exactly the set
        // of positions with `floor(p * blocks / d) == b`, i.e. the CEIL
        // tiling `[ceil(d*b/blocks), ceil(d*(b+1)/blocks))`. The previous
        // floor tiling (`[d*b/blocks, d*(b+1)/blocks)`) disagreed with the
        // run assignment at boundary positions whenever `d % blocks != 0`
        // (e.g. d=14, blocks=4: floor tiling puts p=10 in block 3, the run
        // formula in block 2), splitting one position's displacement-tail
        // run and node-term runs across two concurrently running block
        // tasks — an unordered add race on that accumulator position. The
        // ceil tiling tiles [0, d) identically and agrees with the formula
        // everywhere, so every position's runs land in ONE task and the
        // per-position accumulation order is the deterministic term order.
        // Quantize the block partition to even position counts: the block
        // arrays are 4-lane interleaved (2 f64 positions per 64-byte cache
        // line), so odd boundaries make adjacent block tasks share lines.
        // With even `step` every task's segment starts and ends on a line
        // boundary (only the final block may end short at `d`). This also
        // fixes the partition ONCE for the gather ranges, the runs below
        // (`for_each_position_run`, `first / step`) and the block-level
        // zero runs — one consistent ownership rule, no cross-task adds.
        let step = d.div_ceil(blocks).div_ceil(2) * 2;
        let blocks_q = d.div_ceil(step);
        let ranges: Vec<(u32, u32)> = (0..blocks_q)
            .map(|b| {
                (
                    (b * step).min(d) as u32,
                    ((b + 1) * step).min(d) as u32,
                )
            })
            .collect();

        // Block work lists — same construction as `fold_batch` (term order,
        // live slots only, fused zeroing), with run sources pointing at SoA
        // slots (lane `l` of slot `s` at element `s*LANES + l`).
        let mut block_work: Vec<BlockWork<U>> = (0..blocks_q)
            .map(|_| BlockWork {
                accum: Vec::new(),
                zero: Vec::new(),
            })
            .collect();
        let mut node_terms: Vec<Vec<u32>> = vec![Vec::new(); internal];
        for (ti, (source, _)) in self.structure.terms.iter().enumerate() {
            if let TermSource::Node(k) = source {
                node_terms[*k as usize - 2].push(ti as u32);
            }
        }
        for (ti, (source, _)) in self.structure.terms.iter().enumerate() {
            match source {
                TermSource::Displacement => {
                    // Same argument as `fold_batch`: every skipped prefix
                    // position is a value-zero add (the sign of zero is
                    // unobservable on accumulators — see the bit-exactness
                    // note in the module docs of the tournament).
                    let tail = deg_start[max_deg];
                    let mut sup: Vec<u32> = b_nz_cls
                        .iter()
                        .copied()
                        .map(|p| p as u32)
                        .filter(|&p| p < tail)
                        .collect();
                    sup.sort_unstable();
                    for_each_position_run(&sup, step, |b, first, len| {
                        block_work[b].accum.push(AccumRun {
                            term: ti as u32,
                            g0: first,
                            len,
                            src: unsafe { b_data.add(first as usize * LANES) },
                            zero: false,
                        });
                    });
                    for (b, &(c0, c1)) in ranges.iter().enumerate() {
                        let lo = c0.max(tail);
                        if lo < c1 {
                            block_work[b].accum.push(AccumRun {
                                term: ti as u32,
                                g0: lo,
                                len: c1 - lo,
                                src: unsafe { b_data.add(lo as usize * LANES) },
                                zero: false,
                            });
                        }
                    }
                }
                TermSource::Node(k) => {
                    let ki = *k as usize - 2;
                    let (ptr, _) = compact_ptrs[ki];
                    let sup: Vec<u32> = self.nonzeros[ki].iter().map(|&p| p as u32).collect();
                    debug_assert!(sup.windows(2).all(|w| w[0] < w[1]));
                    let fused = node_terms[ki].len() == 1;
                    debug_assert!(!fused || node_terms[ki][0] == ti as u32);
                    for &CompactSlice { base, rs, re } in &layouts[ki] {
                        let s0 = sup.partition_point(|&p| p < rs);
                        let s1 = sup.partition_point(|&p| p < re);
                        for_each_position_run(&sup[s0..s1], step, |b, first, len| {
                            block_work[b].accum.push(AccumRun {
                                term: ti as u32,
                                g0: first,
                                len,
                                src: unsafe { ptr.add((base + first - rs) as usize * LANES) },
                                zero: fused,
                            });
                        });
                    }
                }
            }
        }
        // Zero runs for the slots no fused accum run zeroes (interior /
        // multi-rooted nodes) — SoA-addressed like the accum runs.
        #[cfg(debug_assertions)]
        let mut zero_covered: Vec<u64> = vec![0; internal];
        for (ni, &(ptr, _)) in compact_ptrs.iter().enumerate() {
            if node_terms[ni].len() == 1 {
                continue;
            }
            let sup: Vec<u32> = self.nonzeros[ni].iter().map(|&p| p as u32).collect();
            debug_assert!(sup.windows(2).all(|w| w[0] < w[1]));
            for &CompactSlice { base, rs, re } in &layouts[ni] {
                let s0 = sup.partition_point(|&p| p < rs);
                let s1 = sup.partition_point(|&p| p < re);
                #[cfg(debug_assertions)]
                let covered = &mut zero_covered[ni];
                for_each_position_run(&sup[s0..s1], step, |b, first, len| {
                    #[cfg(debug_assertions)]
                    {
                        *covered += len as u64;
                    }
                    block_work[b].zero.push(ZeroRun {
                        ptr,
                        s0: (base + first - rs) * LANES as u32,
                        len: len * LANES as u32,
                    });
                });
            }
        }
        #[cfg(debug_assertions)]
        for (ni, sup) in self.nonzeros.iter().enumerate() {
            let fused = node_terms[ni].len() == 1;
            debug_assert_eq!(
                zero_covered[ni] as usize + if fused { sup.len() } else { 0 },
                sup.len(),
                "node {ni}: live-slot zero coverage incomplete"
            );
        }
        let work_shared = SendBlockWork(Arc::new(block_work));
        let ranges: Arc<[(u32, u32)]> = ranges.into();

        // Per-(step, lane) displacement views for the gather task. Raw
        // pointers targeting the lanes' displacement storage, read-only
        // for the whole dispatch (same Send/Sync argument as
        // `SendRhsTable`). Inactive (step, lane) pairs are null.
        let rhs_views_shared = SendRhsViews(Arc::new(
            (0..steps)
                .map(|f| {
                    let mut row: [RhsView<U>; 4] = [
                        RhsView { ptr: std::ptr::null(), len: 0 },
                        RhsView { ptr: std::ptr::null(), len: 0 },
                        RhsView { ptr: std::ptr::null(), len: 0 },
                        RhsView { ptr: std::ptr::null(), len: 0 },
                    ];
                    for (l, lane) in lanes.iter().enumerate() {
                        if f < lane.rhss.len() {
                            row[l] = RhsView {
                                ptr: lane.rhss[f].as_ptr(),
                                len: lane.rhss[f].len(),
                            };
                        }
                    }
                    row
                })
                .collect(),
        ));
        // The per-step active masks, shared into the accumulate task.
        let active_masks_shared: Arc<[u8]> = active_masks.into();
        // "All lanes active" at the cohort's own width (a 2- or 3-lane
        // cohort never flips the vector accumulate on at 0b1111).
        let all_mask = (1u8 << lane_n) - 1;

        let gather_task: Box<dyn Fn(usize, usize) + Send + Sync> = {
            let work = SendBlockWork(Arc::clone(&work_shared.0));
            let b = SendRaw(b_data);
            let ranges = Arc::clone(&ranges);
            let views = SendRhsViews(Arc::clone(&rhs_views_shared.0));
            Box::new(move |bi: usize, fold: usize| {
                // Binding references first forces whole-value captures (see
                // `fold_batch`'s gather task for the precise-capture note).
                let (work, b, ranges, views) = (&work, &b, &ranges, &views);
                // The lead stage (fold 0) also zeroes this block's slot
                // runs of the compact buffers — the accumulate stages
                // maintain the zeroing from there on.
                if fold == 0 {
                    for run in &work.0[bi].zero {
                        unsafe {
                            for si in run.s0 as usize..(run.s0 + run.len) as usize {
                                *run.ptr.add(si) = U::default();
                            }
                        }
                    }
                }
                let (c0, c1) = (ranges[bi].0 as usize, ranges[bi].1 as usize);
                let row = &views.0[fold];
                #[allow(clippy::needless_range_loop)]
                unsafe {
                    for g in c0..c1 {
                        let pw = perm[g] as usize;
                        for l in 0..lane_n {
                            if row[l].ptr.is_null() {
                                continue;
                            }
                            let v = if pw < row[l].len {
                                (*row[l].ptr.add(pw)).clone()
                            } else {
                                U::default()
                            };
                            *b.0.add(g * LANES + l) = v;
                        }
                    }
                }
            })
        };
        let accum_task: Box<dyn Fn(usize, usize) + Send + Sync> = {
            let structure = Arc::clone(&self.structure);
            let acc = SendRaw(acc_data);
            let work = SendBlockWork(Arc::clone(&work_shared.0));
            let masks = Arc::clone(&active_masks_shared);
            Box::new(move |bi: usize, fold: usize| {
                // Whole-value captures (see the gather task above).
                let (work, acc, structure, masks) = (&work, &acc, &structure, &masks);
                let terms = &structure.terms;
                let mask = masks[fold];
                unsafe {
                    // This block's intersecting runs only, in the original
                    // term order: each position's add sequence matches the
                    // full walk bit for bit. All-lanes-active steps run the
                    // 4-lane vector accumulate; masked steps fall back to
                    // per-lane scalar adds for the active lanes (a
                    // zero-multiplier mask could turn an overflowing stale
                    // lane value into a NaN add, diverging from the scalar
                    // engine's skip).
                    for run in &work.0[bi].accum {
                        let weight = &terms[run.term as usize].1;
                        let g0 = run.g0 as usize;
                        if mask == all_mask {
                            for (g, si) in
                                (g0..g0 + run.len as usize).zip(0..run.len as usize)
                            {
                                // `cohort::add_mul4` replicates the scalar
                                // `acc[g] += raw_mul(src[si], weight)` per
                                // lane, bit for bit.
                                lie_rs::add_mul4(
                                    acc.0.add(g * LANES),
                                    run.src.add(si * LANES),
                                    weight,
                                );
                            }
                        } else {
                            for (g, si) in
                                (g0..g0 + run.len as usize).zip(0..run.len as usize)
                            {
                                for l in 0..lane_n {
                                    if mask >> l & 1 == 0 {
                                        continue;
                                    }
                                    lie_rs::raw_add_assign_ptr(
                                        acc.0.add(g * LANES + l),
                                        &lie_rs::raw_mul(&*run.src.add(si * LANES + l), weight),
                                    );
                                }
                            }
                        }
                        // Fused zero-after-add: the run's source slots are
                        // dead once the (single referencing) term has read
                        // them — zero all four lanes' elements (masked
                        // lanes' values are dead too: this term is the
                        // slot's only reader).
                        if run.zero {
                            for si in 0..run.len as usize {
                                let sp = run.src.add(si * LANES);
                                for l in 0..LANES {
                                    *sp.add(l) = U::default();
                                }
                            }
                        }
                    }
                    // zero this block's slot runs of the compact buffers
                    for run in &work.0[bi].zero {
                        for si in run.s0 as usize..(run.s0 + run.len) as usize {
                            *run.ptr.add(si) = U::default();
                        }
                    }
                }
            })
        };

        let shape = BlockShape::new(blocks_q);
        let mut block_stages: Vec<ClassBatchStage<U>> = Vec::with_capacity(2 * steps + 1);
        block_stages.push(ClassBatchStage::block_stage(&shape, &*gather_task, 0));
        for f in 0..steps {
            block_stages.push(ClassBatchStage::block_stage(&shape, &*gather_task, f));
            block_stages.push(ClassBatchStage::block_stage(&shape, &*accum_task, f));
        }
        let mut stages: Vec<&ClassBatchStage<U>> =
            Vec::with_capacity(block_stages.len() + sweep_stages_soa.len() * steps);
        stages.push(&block_stages[0]);
        let sweep_refs_soa: Vec<&ClassBatchStage<U>> = sweep_stages_soa.iter().collect();
        for f in 0..steps {
            stages.push(&block_stages[1 + 2 * f]);
            for s in &sweep_refs_soa {
                stages.push(s);
            }
            stages.push(&block_stages[2 + 2 * f]);
        }

        run_class_batch_with_work(series, order, &stages, steps, planned_scalar * lane_n);

        // Epilogue: each lane's class-space accumulator back to public basis
        // order (assignment — it holds the full sum).
        let acc4 = acc_cell.into_inner();
        for (l, lane) in lanes.iter_mut().enumerate() {
            for (k, &src) in inv.iter().enumerate().take(lane.acc.len()) {
                lane.acc[k] = acc4[src as usize * LANES + l].clone();
            }
        }

        // The lists were used with (common accumulator support, common
        // displacement support) throughout; record those as the atom
        // supports so the next per-fold comparison is against reality.
        self.atom_a.clear();
        self.atom_a.extend_from_slice(a_nonzero);
        self.atom_b.clear();
        self.atom_b.extend_from_slice(b_nonzero);
        self.lists_built = true;
    }
}
