//! The scalar batch fold: one stage chain per batch that gathers each
//! displacement into class space, sweeps the DAG's levels and accumulates
//! the terms.

use lie_rs::{
    BlockShape, ClassBatchStage, KernelJob, LieSeries, plan_class_sweep_stages,
    planned_sweep_entries, run_class_batch, work_adaptive_slots,
};
use lyndon_rs::generators::Generator;
use num_traits::{One, Zero};
use std::cell::UnsafeCell;
use std::hash::Hash;
use std::ops::{AddAssign, MulAssign, Neg};
use std::sync::Arc;

use super::block_work::{
    AccumRun, BlockWork, SendBlockWork, SendRaw, SendRhsTable, ZeroRun, for_each_position_run,
};
use super::{CommutatorDag, DagNode, SharedNodeLists, TermSource};

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
    /// Folds `rhss` (public-order displacement slices) into the class-space
    /// accumulator `acc` (public order, in-out) as ONE continuous stage
    /// chain: per fold, a gather stage brings the displacement into class
    /// space, the DAG's levels sweep (reading the class-space accumulator
    /// and the shared displacement buffer), and an accumulate stage adds
    /// every term into the class-space accumulator, zeroing the node
    /// buffers afterwards for the next fold. One plan and one job table
    /// serve the whole batch; the accumulator never leaves class space, so
    /// the per-fold public→class→public round-trip disappears and the
    /// plan/gating cost is paid once, not per fold.
    ///
    /// Contract (checked by the driver): the accumulator's support equals
    /// the reachable set (the union of the node lists and the displacement
    /// support, see [`Self::batch_eligible`]) and all displacements share
    /// the support `b_nonzero`. Level-0 gating then uses masks that cover
    /// every position any later fold can touch, which stays sound even if
    /// mid-batch accumulation cancels values to zero — the node lists are
    /// scatter-target supersets, so a shrinking support only loses
    /// tightness, never targets.
    pub(crate) fn fold_batch<T>(
        &mut self,
        series: &mut LieSeries<T, U>,
        rhss: &[&[U]],
        a_nonzero: &[usize],
        b_nonzero: &[usize],
    ) where
        T: Clone + Ord + Generator + Hash + Eq,
        U: Send + Sync,
    {
        let internal = self.structure.nodes.len() - 2;
        self.ensure_buffers(series.coefficients.len(), internal);
        let order = self
            .class_order
            .get()
            .expect("batch fold requires the class order of a prior steady fold");
        let inv = order.inv();
        let perm = order.perm();
        let d = series.coefficients.len();

        // The batch's class-space working state. The accumulator stays
        // class-ordered across folds: fold f's sweep reads it as its `a`
        // operand and fold f's accumulate adds the terms in place. Both
        // buffers are written by block stages and read by the sweeps'
        // operand slices; all access goes through the cells' raw pointers
        // (the UnsafeCell makes that interleaving defined) and the stage
        // counters order every write-before-read across stages.
        let mut acc_cls = vec![U::default(); d];
        for (w, c) in series.coefficients.iter().enumerate().take(d) {
            acc_cls[inv[w] as usize] = c.clone();
        }
        let acc_cell = UnsafeCell::new(acc_cls);
        let b_cell = UnsafeCell::new(vec![U::default(); d]);
        // SAFETY: the cells are not moved for the dispatch's duration, so
        // the heap buffers behind their Vecs are stable; everything below
        // goes through these data pointers (the sweeps read shared slices
        // derived from them, the block stages read and write through
        // them), and the stage counters order every write-before-read.
        let acc_data: *mut U = unsafe { (*acc_cell.get()).as_mut_ptr() };
        let b_data: *mut U = unsafe { (*b_cell.get()).as_mut_ptr() };

        // Compact per-node scratch buffers. A node's sweep only ever
        // scatters into the degree slices its recorded support touches
        // (averaging ~1.5 slices per node), and the class space is
        // degree-grouped, so storing whole slices costs ~4-6x less memory
        // than one full-d buffer per node — deep DAGs (hundreds of nodes)
        // drop from megabytes to tens/hundreds of KB and the sweeps' scatter
        // RMWs become L2-resident. Layout per node: its active degree
        // slices concatenated in degree order; class position `x` of degree
        // `δ` lives at `x - shift[δ]` (the per-degree shift table the jobs
        // carry; a full-d buffer uses the all-zero identity table).
        // Per-class-position Lyndon degrees (the class space is
        // degree-grouped, so this is non-decreasing). NOTE: class positions,
        // not the public basis order — the batch works entirely in the
        // class-contiguous space.
        let pos_degree: Vec<u8> = order.degree_cls().to_vec();
        let max_deg = pos_degree.iter().copied().max().unwrap_or(0) as usize;
        // Degree-slice starts: the class space is degree-grouped, so slice
        // `dd` = [first(dd), first(dd+1)) where first(x) is the first class
        // position of degree x (d when degree x is empty; backfilled to the
        // next populated degree's start so empty slices are zero-length).
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

        // Identity shifts for the full-d atom buffers (acc_cls / b_cls).
        // The compact per-node buffers are built after the sweep planner
        // returns the jobs' exact batch scatter sets (the recorded node
        // lists are only a union-level fixed point — the batch eligibility
        // checks the union, not each list — so a level-0 job's list,
        // recorded under an earlier, smaller accumulator support, can miss
        // positions the batch's first fold scatters; sizing from the
        // recorded lists would be unsound).
        let zero_shifts = vec![0u32; max_deg + 2];
        let zero_shifts_ptr: *const u32 = zero_shifts.as_ptr();

        // One job table for the whole batch: level-0 jobs read the
        // class-space accumulator and displacement buffers (fixed
        // allocations, values refreshed in-graph by the gather stages)
        // with the reachable/common supports; interior jobs read child
        // node buffers with the children's scatter-target lists — exactly
        // the per-fold operand fill. The jobs' `result`/`r_shift` are
        // wired after the planner returns the exact per-node scatter sets.
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
                            a_len: d,
                            a_nonzero: &a_nz_cls,
                            b: b_data,
                            b_len: d,
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
        // Local copies of the node lists for the jobs' operand views (the
        // lists are re-recorded from the planner's exact scatter sets
        // below, so the jobs must not borrow `self.nonzeros`). `Arc`-shared
        // read-only: the re-record below installs a NEW Arc wholesale, so
        // the jobs' view of the old lists stays valid without a copy.
        let lists_local: SharedNodeLists = Arc::clone(&self.nonzeros);
        for (li, level) in self.structure.levels.iter().enumerate().skip(1) {
            for (jj, &k) in level.iter().enumerate() {
                let (left, right) = match self.structure.nodes[k as usize] {
                    DagNode::Binary { left, right } => (left, right),
                    DagNode::Atom(_) => continue,
                };
                let (a_ptr, a_len, a_sh, a_nz): (*const U, usize, *const u32, &[usize]) = match left
                {
                    0 => (acc_data as *const U, d, zero_shifts_ptr, &a_nz_cls),
                    1 => (b_data as *const U, d, zero_shifts_ptr, &b_nz_cls),
                    id => (
                        std::ptr::null(),
                        0,
                        zero_shifts_ptr,
                        &lists_local[id as usize - 2],
                    ),
                };
                let (b_ptr, b_len, b_sh, b_nz): (*const U, usize, *const u32, &[usize]) =
                    match right {
                        0 => (acc_data as *const U, d, zero_shifts_ptr, &a_nz_cls),
                        1 => (b_data as *const U, d, zero_shifts_ptr, &b_nz_cls),
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

        // Plan once to get the jobs' exact batch scatter sets (the gating
        // depends only on the fixed support lists, so this plan's stages
        // are throwaway — they borrow the placeholder jobs and are dropped
        // before the jobs are rewired; the final plan below hits the
        // gating cache).
        let (_, scatter_sets) =
            plan_class_sweep_stages(
                series,
                order,
                &levels_jobs,
                &mut self.gating_cache,
                true,
                false,
            );

        // Record the true scatter sets as the node lists: they bound this
        // batch's sweeps exactly, keep the union-level eligibility fixed
        // point for later batches, and stay sound supersets for the per-
        // fold path's gating.
        {
            // Copy-on-write re-record: the exact scatter sets replace EVERY
            // binary node's list (interning only creates binary internal
            // nodes — asserted below), so a fresh shell avoids copying the
            // old lists — crucial when they are an adopted cache snapshot
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

        // Compact per-node scratch buffers, sized from the exact scatter
        // sets: per node, its support's degree slices concatenated in
        // degree order; class position `x` of degree `δ` lives at
        // `x - shift[δ]`.
        struct CompactSlice {
            base: u32,
            rs: u32,
            re: u32,
        }
        // Indexed BY NODE ID (`k - 2`), not by (level, job) position — the
        // DAG's levels are not id-ordered (a late-interned shallow node can
        // share a level with early deep nodes), so (level, job) iteration
        // order must not be conflated with node identity.
        let internal = self.buffers.len();
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
                    let base = buf.len() as u32;
                    buf.resize(buf.len() + (re - rs) as usize, U::default());
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

        // Wire the jobs' results and the interior jobs' operand views to
        // the compact buffers (the nodes' slices cover their recorded
        // supports, which is exactly what the presence tests load).
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

        // Final plan with the wired jobs: the gating cache makes this cheap
        // (the support lists are unchanged); the stages reference the final
        // job table the sweeps read through.
        let (sweep_stages, _) =
            plan_class_sweep_stages(
                series,
                order,
                &levels_jobs,
                &mut self.gating_cache,
                true,
                false,
            );

        // Block ranges over class positions. The block count derives from
        // the SAME work-adaptive policy the walk will use (below): blocks
        // exist to spread the gather/accumulate phases across slots, so
        // cutting more blocks than the policy's slot count only multiplies
        // the walk's per-pack claim/publish protocol (78 micro-blocks per
        // stage × ~1000 folds of atomic RMWs dwarfed the gather itself and
        // kept the tiny-grid e2e 2.3x above its 1t floor even at slots=1).
        // The policy sees the per-fold sweep work; the two block stages a
        // fold adds are counted at their post-cut element count, which for
        // the slot decision changes nothing (both terms are « QUANTUM
        // exactly when the sweep term is).
        let threads = rayon::current_num_threads().max(1);
        let sweep_refs: Vec<&ClassBatchStage<U>> = sweep_stages.iter().collect();
        let slots = work_adaptive_slots(threads, d.max(1), planned_sweep_entries(&sweep_refs));
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

        // Block closures, kept alive for the whole dispatch. Z zeroes the

        // compact node buffers before the first in-batch sweep (the
        // accumulate stages then maintain the zeroing); per fold, BG gathers
        // the displacement into class space and C accumulates the terms.
        //
        // Those runs are materialized once per batch, here at plan time —
        // and only over each source's LIVE slots. A node's sweep scatters
        // exactly onto its recorded scatter set (the planner's exact write
        // set, re-recorded into `self.nonzeros` above), so compact slots
        // outside it are padding that no sweep ever touches — permanently
        // `U::default()`. Their `weight × (+0.0)` adds are value-zero and
        // dropped: the only observable difference is the sign bit of an
        // exact-zero accumulator slot (adding +0.0 turns a −0.0 into +0.0;
        // adding −0.0 is a bitwise identity), never the `==` value. The
        // runs keep the full walk's per-position add order (ascending term
        // index; within a term, ascending degree slice and ascending
        // position), so every position's add sequence — value, weight, term
        // order — is unchanged. Zeroing is fused into the single
        // referencing term's runs (a node is the canonical root of at most
        // one BCH term — distinct Lyndon words have distinct canonical
        // bracketings — so no other run reads the slot after that term's
        // add); nodes no term accumulates (shared interior brackets, or
        // defensively a multi-rooted node) keep explicit zero runs over
        // their live slots.
        let mut block_work: Vec<BlockWork<U>> = (0..blocks_q)
            .map(|_| BlockWork {
                accum: Vec::new(),
                zero: Vec::new(),
            })
            .collect();
        // Referencing terms per node: count 1 → that term's accum runs fuse
        // the slot zeroing; count 0 (interior node) or ≥ 2 (defensive) →
        // explicit zero runs.
        let mut node_terms: Vec<Vec<u32>> = vec![Vec::new(); internal];
        for (ti, (source, _)) in self.structure.terms.iter().enumerate() {
            if let TermSource::Node(k) = source {
                node_terms[*k as usize - 2].push(ti as u32);
            }
        }
        for (ti, (source, _)) in self.structure.terms.iter().enumerate() {
            match source {
                TermSource::Displacement => {
                    // The gather rewrites every b_cls slot each fold with
                    // `rhs[perm[g]]` (`U::default()` beyond the rhs
                    // length), and the batch driver guarantees every rhs is
                    // value-zero outside the common support over the
                    // kernel-reachable positions
                    // `[0, degree_start(max_degree))` — exactly the prefix
                    // below `tail` (class and public layouts share degree-
                    // slice boundaries). The degree-`max_degree` tail lies
                    // beyond that check, so its whole contiguous range
                    // stays covered: its gather-written values fold
                    // unconditionally, exactly as before. Everything
                    // skipped in the prefix is a value-zero add (see the
                    // sign-of-zero note above).
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
                            // b_cls is a full-d buffer: class position `p`
                            // lives at slot `p`.
                            src: unsafe { b_data.add(first as usize) },
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
                                src: unsafe { b_data.add(lo as usize) },
                                zero: false,
                            });
                        }
                    }
                }
                TermSource::Node(k) => {
                    let ki = *k as usize - 2;
                    let (ptr, _) = compact_ptrs[ki];
                    // The node's exact batch scatter set (sorted class
                    // positions): the positions its sweep writes, hence
                    // the only slots that can be non-zero here.
                    let sup: Vec<u32> = self.nonzeros[ki].iter().map(|&p| p as u32).collect();
                    debug_assert!(sup.windows(2).all(|w| w[0] < w[1]));
                    let fused = node_terms[ki].len() == 1;
                    debug_assert!(!fused || node_terms[ki][0] == ti as u32);
                    for &CompactSlice { base, rs, re } in &layouts[ki] {
                        // The support positions inside the slice (both
                        // lists ascending); position `p` lives at slot
                        // `base + p - rs` and belongs to block
                        // `p * blocks / d`.
                        let s0 = sup.partition_point(|&p| p < rs);
                        let s1 = sup.partition_point(|&p| p < re);
                        for_each_position_run(&sup[s0..s1], step, |b, first, len| {
                            block_work[b].accum.push(AccumRun {
                                term: ti as u32,
                                g0: first,
                                len,
                                src: unsafe { ptr.add((base + first - rs) as usize) },
                                zero: fused,
                            });
                        });
                    }
                }
            }
        }
        // Zero runs for the slots no fused accum run zeroes: interior
        // nodes (never accumulated — shared brackets read only as sweep
        // operands) and, defensively, multi-rooted nodes.
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
                        s0: base + first - rs,
                        len,
                    });
                });
            }
        }
        // Plan-time coverage invariant: every live slot of every node is
        // zeroed exactly once per fold — by its single referencing term's
        // fused runs (which cover the node's whole scatter set) or by the
        // explicit zero runs above.
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
        // The batch's TWO block tasks (gather and accumulate), created once
        // per batch and shared by every fold's stages: each stage binds its
        // fold via `ClassBatchStage::block_stage`'s fold index, and the
        // tasks read that fold's inputs through the shared tables. This
        // keeps the per-fold stage setup allocation-free — the per-fold
        // closures (one Box + a ranges clone + a debug Arc each, plus the
        // per-stage task/pack vectors) this replaces dominated the
        // small-grid allocator traffic.
        let ranges: Arc<[(u32, u32)]> = ranges.into();
        let rhs_table = SendRhsTable(
            rhss.iter()
                .map(|r| (r.as_ptr(), r.len()))
                .collect(),
        );
        let gather_task: Box<dyn Fn(usize, usize) + Send + Sync> = {
            let work = SendBlockWork(Arc::clone(&work_shared.0));
            let b = SendRaw(b_data);
            let ranges = Arc::clone(&ranges);
            let rhs_table = rhs_table;
            Box::new(move |bi: usize, fold: usize| {
                // Binding references first forces whole-value captures: a
                // direct `rhs.add(..)` place access would let the precise
                // capture split off the bare raw-pointer field, losing the
                // Send/Sync wrapper's promise.
                let (work, b, ranges, rhs_table) = (&work, &b, &ranges, &rhs_table);
                // The lead stage (fold 0) also zeroes this block's slot
                // runs of the compact buffers — the accumulate stages
                // maintain the zeroing from there on (see the stage-chain
                // comment above).
                if fold == 0 {
                    for run in &work.0[bi].zero {
                        unsafe {
                            for si in run.s0 as usize..(run.s0 + run.len) as usize {
                                *run.ptr.add(si) = U::default();
                            }
                        }
                    }
                }
                let (rhs, rhs_len) = rhs_table.0[fold];
                let (c0, c1) = (ranges[bi].0 as usize, ranges[bi].1 as usize);
                #[allow(clippy::needless_range_loop)]
                unsafe {
                    for g in c0..c1 {
                        let pw = perm[g] as usize;
                        *b.0.add(g) = if pw < rhs_len {
                            (*rhs.add(pw)).clone()
                        } else {
                            U::default()
                        };
                    }
                }
            })
        };
        let accum_task: Box<dyn Fn(usize, usize) + Send + Sync> = {
            let structure = Arc::clone(&self.structure);
            let acc = SendRaw(acc_data);
            let work = SendBlockWork(Arc::clone(&work_shared.0));
            Box::new(move |bi: usize, _fold: usize| {
                // Whole-value captures (see the gather task above).
                let (work, acc, structure) = (&work, &acc, &structure);
                let terms = &structure.terms;
                unsafe {
                    // This block's intersecting runs only, in the
                    // original term order: each position's add sequence
                    // matches the full walk bit for bit. Fused runs
                    // zero their source slot right after consuming it
                    // — the node's slots are dead once its single
                    // referencing term has read them — so the fold's
                    // zeroing touches the same cache lines the adds
                    // already do.
                    for run in &work.0[bi].accum {
                        let weight = &terms[run.term as usize].1;
                        let g0 = run.g0 as usize;
                        // SAFETY: the run's positions stay inside the
                        // block's disjoint range of `acc` and inside the
                        // source buffer's run of live slots; the fused
                        // zero-back writes only slots this run just read
                        // (owned by this block, dead afterwards).
                        if run.zero {
                            for (g, si) in
                                (g0..g0 + run.len as usize).zip(0..run.len as usize)
                            {
                                // `raw_mul`/`raw_add_assign_ptr` skip the
                                // float wrappers' per-op NaN checks
                                // (raw-float fast path, see lie-rs
                                // `raw_mul`'s NaN policy).
                                lie_rs::raw_add_assign_ptr(
                                    acc.0.add(g),
                                    &lie_rs::raw_mul(&*run.src.add(si), weight),
                                );
                                *run.src.add(si) = U::default();
                            }
                        } else {
                            for (g, si) in
                                (g0..g0 + run.len as usize).zip(0..run.len as usize)
                            {
                                lie_rs::raw_add_assign_ptr(
                                    acc.0.add(g),
                                    &lie_rs::raw_mul(&*run.src.add(si), weight),
                                );
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
        let mut block_stages: Vec<ClassBatchStage<U>> = Vec::with_capacity(2 * rhss.len() + 1);
        // Lead stage: fold 0's gather carries the compact buffers'
        // initial zeroing (the gather task zeroes when fold == 0); fold
        // 0's own gather stage then re-runs idempotently, exactly as the
        // per-fold lead closure it replaces.
        block_stages.push(ClassBatchStage::block_stage(&shape, &*gather_task, 0));
        for f in 0..rhss.len() {
            block_stages.push(ClassBatchStage::block_stage(&shape, &*gather_task, f));
            block_stages.push(ClassBatchStage::block_stage(&shape, &*accum_task, f));
        }
        let mut stages: Vec<&ClassBatchStage<U>> =
            Vec::with_capacity(block_stages.len() + sweep_stages.len() * rhss.len());
        stages.push(&block_stages[0]);
        for f in 0..rhss.len() {
            stages.push(&block_stages[1 + 2 * f]);
            for s in &sweep_stages {
                stages.push(s);
            }
            stages.push(&block_stages[2 + 2 * f]);
        }

        // Fold units for the slot policy: the stage chain repeats its
        // gather/sweep/accumulate sub-chain once per displacement (plus
        // one leading zero stage), so each unit's barrier cost is paid
        // per fold — the policy normalizes by the unit count.
        run_class_batch(series, order, &stages, rhss.len().max(1));

        // Epilogue: the class-space accumulator back to public basis order
        // (assignment — it holds the full sum).
        let acc_cls = acc_cell.into_inner();
        for (k, &src) in inv.iter().enumerate().take(series.coefficients.len()) {
            series.coefficients[k] = acc_cls[src as usize].clone();
        }

        // The lists were used with (dense accumulator, common displacement
        // support) throughout; record those as the atom supports so the
        // next per-fold comparison is against reality.
        self.atom_a.clear();
        self.atom_a.extend_from_slice(a_nonzero);
        self.atom_b.clear();
        self.atom_b.extend_from_slice(b_nonzero);
        self.lists_built = true;
    }
}
