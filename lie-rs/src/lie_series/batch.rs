use super::*;
use super::gating::KernelGating;
use super::kernel::NoSink;

/// Target size of a parallel bundle, in visited entries.
const BUNDLE_TARGET_ENTRIES: usize = 2048;
/// Floor for the thread-adaptive bundle target, in visited entries: below
/// this a bundle's work is smaller than its dispatch cost, and unit
/// integrity (bundles never split a unit) provides the natural minimum.
const MIN_BUNDLE_ENTRIES: usize = 16;

/// Shared handle for the parallel sweep's result writes.
///
/// SAFETY (of the `Send`/`Sync` impls): bundles of units own disjoint
/// word sets — a basis word's (degree, content) determines the single
/// unit whose entries write it, and bundles never split a unit — so
/// concurrent read-modify-writes through `ptr` touch disjoint addresses.
/// `U: Send` covers the values moving across threads.
struct RawResult<'a, U> {
    ptr: *mut U,
    _marker: PhantomData<&'a mut [U]>,
}

unsafe impl<U: Send> Send for RawResult<'_, U> {}
unsafe impl<U: Send> Sync for RawResult<'_, U> {}

/// Chooses the entry table for a sweep. Monomorphized per layout, so each
/// path's inner loop compiles to exactly one code shape.
trait ScatterLayout {
    fn entries<'a>(class: Option<&'a ClassOrder>, direct: &'a [Entry]) -> &'a [Entry];
}

/// Public basis order: public entries, writes straight into the job's
/// result buffer.
struct DirectLayout;

/// Fully class-contiguous order: class-ordered operands (read through the
/// relabeled entry table) and a class-ordered result buffer — no scratch,
/// no permutation.
struct ClassInternalLayout;

impl ScatterLayout for DirectLayout {
    #[inline(always)]
    fn entries<'a>(_class: Option<&'a ClassOrder>, direct: &'a [Entry]) -> &'a [Entry] {
        direct
    }
}

impl ScatterLayout for ClassInternalLayout {
    #[inline(always)]
    fn entries<'a>(class: Option<&'a ClassOrder>, _direct: &'a [Entry]) -> &'a [Entry] {
        class
            .expect("class layout without a class order")
            .entries_cls()
    }
}

/// Parallel single-phase per-unit sweep: each bundle's units are visited
/// independently — one term per active entry, its row contributions added
/// straight into the result buffer (no intermediate term buffers, no phase
/// barrier). Units are atomic and their word sets partition the write
/// slots, so concurrent bundles never touch the same word; each word's
/// adds happen inside its one unit in table-entry order, which is the
/// per-word float summation sequence the serial sweep produces, so the
/// bits are unchanged. `L` selects the entry table layout — public
/// order into the job's result buffer, or class-contiguous order into a
/// class-ordered result buffer.
fn sweep_bundles_parallel<L: ScatterLayout, T, U>(
    jobs: &[KernelJob<'_, U>],
    writers: &[RawResult<U>],
    gateways: &[KernelGating],
    entries_slice: &[Entry],
    decomp_rels: &[u32],
    class_order: Option<&ClassOrder>,
    decomp_coeffs_slice: &[U],
    unit_bundles: &[Vec<(usize, u32)>],
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
    use rayon::prelude::*;

    #[cfg(feature = "tracing")]
    let sweep_span = tracing::debug_span!(
        "kernel_sweep_parallel",
        jobs = jobs.len(),
        bundles = unit_bundles.len(),
        threads = rayon::current_num_threads(),
    )
    .entered();
    let entries_tbl = <L as ScatterLayout>::entries(class_order, entries_slice);

    // Single-phase per-unit sweep: each bundle's units are visited
    // independently — one term per active entry (computed inline, exactly
    // the old phase-1 formula) and its row contributions added straight
    // into the result buffer. Units are atomic and their word sets
    // partition the write slots — bundles never split one, so concurrent
    // bundles never touch the same word; each word's adds happen inside
    // its one unit in table-entry order, which is exactly the old
    // two-phase sweep's per-word float summation sequence (and the serial
    // sweep's), so the bits are unchanged. Accumulate-into semantics are
    // preserved (the buffer's current value starts every sum).
    unit_bundles
        .par_iter()
        .enumerate()
        .for_each(|(_bundle_index, bundle)| {
            #[cfg(feature = "tracing")]
            let _bundle_span = tracing::debug_span!(
                "kernel_bundle",
                bundle = _bundle_index,
                units = bundle.len(),
            )
            .entered();
            for &(ji, ui) in bundle {
                let job = &jobs[ji];
                let writer = &writers[ji];
                let gating = &gateways[ji];
                let unit = &gating.units[ui as usize];
                let base = unit.rs as usize;
                for tp in unit.ticket_start as usize..unit.ticket_end as usize {
                    let ticket = gating.tickets[tp];
                    let e = (ticket & TICKET_INDEX_MASK) as usize;
                    let entry = entries_tbl[e];
                    let p_active = ticket & TICKET_P_ACTIVE != 0;
                    let q_active = ticket & TICKET_Q_ACTIVE != 0;
                    let (i, j) = (entry.i as usize, entry.j as usize);
                    // SAFETY: i and j are positions < the operand lengths
                    // (the tickets resolved presence against the operand
                    // supports). `raw_mul` skips the float wrappers'
                    // per-op NaN checks (raw-float fast path); `-` never
                    // checks.
                    let term = unsafe {
                        if p_active {
                            let mut t = raw_mul(&*job.a.add(i), &*job.b.add(j));
                            if q_active {
                                raw_add_assign(
                                    &mut t,
                                    &-raw_mul(&*job.a.add(j), &*job.b.add(i)),
                                );
                            }
                            t
                        } else {
                            -raw_mul(&*job.a.add(j), &*job.b.add(i))
                        }
                    };
                    let from = entry.decomp_start as usize;
                    let to = entries_tbl[e + 1].decomp_start as usize;
                    for (k, &rel) in decomp_rels[from..to].iter().enumerate() {
                        // SAFETY: the unit's word set partitions the
                        // result positions (the ownership invariant) and
                        // this bundle owns the unit whole — single-writer
                        // RMW, in bounds (the gating addresses the working
                        // layout's result space).
                        unsafe {
                            raw_add_assign(
                                &mut *writer.ptr.add(base + rel as usize),
                                &raw_mul(&decomp_coeffs_slice[from + k], &term),
                            );
                        }
                    }
                }
            }
        });
    #[cfg(feature = "tracing")]
    drop(sweep_span);
}

impl<U> RawResult<'_, U>
where
    U: Clone + AddAssign + 'static,
{
}

/// One independent commutation: operand slices plus the destination buffer.
///
/// SAFETY contract for `LieSeries::commutator_coefficients_batch`: the
/// `result` buffers of the jobs passed to one batch call must be pairwise
/// disjoint (in a fold these are distinct DAG-node buffers). Within a job,
/// the anagram partition makes the units conflict-free.
pub struct KernelJob<'a, U> {
    /// The left operand's coefficients, as a raw pointer because the
    /// DAG-fold batch mutates the operand buffers through `UnsafeCell`
    /// between stages while these slices stay live — a shared reference
    /// would be frozen over concurrently-mutated memory (undefined
    /// behavior, and it lets the compiler cache operand reads across
    /// stage boundaries). Valid for `a_len` elements for the duration of
    /// the call.
    pub a: *const U,
    /// The left operand's length.
    pub a_len: usize,
    // SAFETY: `Sync` is sound because the batch's tasks only *read* the
    // operand buffers (through the raw pointers, ordered by the stage
    // counters) and write through `result` at indices owned by exactly
    // one task (disjoint buffers across jobs, anagram-disjoint units
    // within a job); `U: Send` covers the written values crossing threads.
    //
    // (The unsafe impls follow the struct definition.)
    /// The left operand's non-zero indices (superset of its support).
    pub a_nonzero: &'a [usize],
    /// The right operand's coefficients (see `a`).
    pub b: *const U,
    /// The right operand's length.
    pub b_len: usize,
    /// The right operand's non-zero indices.
    pub b_nonzero: &'a [usize],
    /// The destination buffer, as a raw pointer because the batch's parallel
    /// tasks write disjoint regions of it concurrently. Must remain valid
    /// for `result_len` elements for the duration of the call.
    pub result: *mut U,
    /// The destination buffer's length.
    pub result_len: usize,
    /// Compact-layout address shifts, indexed by Lyndon degree (tables must
    /// have at least `max_degree + 1` entries). Class position `x` of degree
    /// `d` lives at `x - shift[d]` in the corresponding buffer; a full-d
    /// buffer uses the all-zero `IDENTITY_SHIFTS` table. The batch
    /// fold's per-node buffers store only the degree slices the node's
    /// sweep writes (4-6x smaller than full-d for deep DAGs), so every
    /// operand/result access subtracts the hoisted per-degree shift.
    pub a_shift: *const u32,
    pub b_shift: *const u32,
    pub r_shift: *const u32,
}

/// The identity shift table for full-d buffers (see [`KernelJob`]).
pub const IDENTITY_SHIFTS: [u32; 128] = [0u32; 128];

// SAFETY: see the comment inside `KernelJob` — tasks write only indices
// owned by their job, and only values of `U: Send`.
unsafe impl<U: Send> Send for KernelJob<'_, U> {}
unsafe impl<U: Send> Sync for KernelJob<'_, U> {}

/// Batch kernel: evaluates several independent commutation jobs in one
/// parallel dispatch. Jobs must have disjoint `result` buffers; per-word
/// accumulation order matches the serial sweep exactly.
pub fn commutator_coefficients_batch<T, U>(
    a_series: &LieSeries<T, U>,
    jobs: &mut [KernelJob<'_, U>],
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
    let mut cache = GatingCache::default();
    commutator_coefficients_batch_with_cache(a_series, jobs, &mut cache);
}

/// [`commutator_coefficients_batch`] with a caller-owned [`GatingCache`]
/// that persists across calls, amortizing the gating scan over repeated
/// jobs with equal degree support.
pub fn commutator_coefficients_batch_with_cache<T, U>(
    a_series: &LieSeries<T, U>,
    jobs: &mut [KernelJob<'_, U>],
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
    // Prologue per job (serial): presence bitsets, degree masks, active
    // units (memoized). Cheap relative to the sweep.
    #[cfg(feature = "tracing")]
    let prologue_span = tracing::debug_span!("kernel_prologue", jobs = jobs.len()).entered();
    let gateways: Vec<KernelGating> = jobs
        .iter()
        .map(|j| {
            LieSeries::<T, U>::kernel_prologue_cached(a_series, j.a_nonzero, j.b_nonzero, cache)
        })
        .collect();
    let total: usize = gateways.iter().map(|g| g.total_pairs).sum();
    #[cfg(feature = "tracing")]
    {
        drop(prologue_span);
        tracing::debug!(
            jobs = jobs.len(),
            total_pairs = total,
            "kernel_prologue done"
        );
    }

    // Hoisted flat-table views: every contribution's ticket indexes the
    // entry table directly, so neither sweep touches the unit table again.
    let table = &a_series.feasible_decompositions;
    let entries_slice = table.entries();
    let decomp_coeffs_slice = table.decomp_coeffs();

    // With more than one thread available the batch always dispatches to
    // the parallel sweep; rayon's work stealing balances the pieces. Only a
    // single-threaded pool runs the serial sweep.
    let threads = rayon::current_num_threads();
    if threads == 1 {
        #[cfg(feature = "tracing")]
        let sweep_span = tracing::debug_span!(
            "kernel_sweep_serial",
            jobs = jobs.len(),
            total_pairs = total,
            threads
        )
        .entered();
        for (job, gating) in jobs.iter_mut().zip(&gateways) {
            // SAFETY: the job's result buffer is valid for `result_len`
            // elements (the struct's contract) and is exclusively ours here.
            let result = unsafe { std::slice::from_raw_parts_mut(job.result, job.result_len) };
            // SAFETY: these serial paths are exclusive per job (no
            // concurrent mutation), so shared slices over the raw operand
            // pointers are sound here.
            let (a_slice, b_slice) = unsafe {
                (
                    std::slice::from_raw_parts(job.a, job.a_len),
                    std::slice::from_raw_parts(job.b, job.b_len),
                )
            };
            // The gating's word classes carry public positions: one
            // single-writer accumulation + store per word, straight into
            // the public result buffer. (The old class-contiguous scratch
            // densified a scatter of repeated `+=`s; one store per word
            // needs no densifying, and the scratch + permutation epilogue
            // it paid is gone.)
            LieSeries::sweep_words_serial(
                a_series,
                a_slice,
                b_slice,
                gating,
                result,
                entries_slice,
                table.decomp_indices_rel(),
                &mut NoSink,
            );
        }
        #[cfg(feature = "tracing")]
        drop(sweep_span);
        return;
    }

    // SAFETY: each job's sweep target is a distinct buffer (disjoint across
    // jobs by the caller's contract), and within a job each word class is
    // written by exactly one task (bundles never split a class), so the
    // concurrent single-writer stores never alias. The buffers are not
    // otherwise accessed during the sweep.
    let writers: Vec<RawResult<U>> = jobs
        .iter()
        .map(|j| RawResult {
            ptr: j.result,
            _marker: PhantomData,
        })
        .collect();

    // Flatten (job, unit) pairs into bundles of roughly
    // `BUNDLE_TARGET_ENTRIES` contributions, weighted by the unit's pair
    // count. Units stay whole within a bundle (target/degree order,
    // preserving each word's accumulation context); a unit is never split.
    let mut bundles: Vec<Vec<(usize, u32)>> = vec![Vec::new()];
    let mut cur = 0usize;
    // Enough bundles that every thread can hold a piece, without dropping
    // below the per-task break-even size. The cut happens between units —
    // always safe: unit word sets never overlap.
    //
    // Balancing by decomposition volume (the alternative: bundles of equal
    // summed decomposition length) was measured and rejected in the
    // unit-atomic engine: per-bundle sweep time is ~uniform per work unit
    // across degree mixes, so counts are already the right balance metric.
    let bundle_target = (total / (2 * threads)).clamp(MIN_BUNDLE_ENTRIES, BUNDLE_TARGET_ENTRIES);
    for (ji, gating) in gateways.iter().enumerate() {
        for (ui, unit) in gating.units.iter().enumerate() {
            let pairs = unit.pairs as usize;
            if cur >= bundle_target {
                bundles.push(Vec::new());
                cur = 0;
            }
            bundles.last_mut().unwrap().push((ji, ui as u32));
            cur += pairs;
        }
    }

    sweep_bundles_parallel::<DirectLayout, T, U>(
        jobs,
        &writers,
        &gateways,
        entries_slice,
        table.decomp_indices_rel(),
        None,
        decomp_coeffs_slice,
        &bundles,
    )
}

/// Batch kernel in the class-contiguous working mode: every job's operand
/// slices are class-ordered, its support lists class-indexed, and its
/// result buffer receives the class-ordered sum directly — no scratch, no
/// permutation. Requires `order` to describe `a_series`' basis.
pub fn commutator_coefficients_class_batch_with_cache<T, U>(
    a_series: &LieSeries<T, U>,
    order: &ClassOrder,
    jobs: &mut [KernelJob<'_, U>],
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

    #[cfg(feature = "tracing")]
    let prologue_span = tracing::debug_span!("kernel_prologue_class", jobs = jobs.len()).entered();
    let gateways: Vec<KernelGating> = jobs
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
        .collect();
    let total: usize = gateways.iter().map(|g| g.total_pairs).sum();
    #[cfg(feature = "tracing")]
    drop(prologue_span);

    let table = &a_series.feasible_decompositions;
    let entries_slice = table.entries();
    let decomp_coeffs_slice = table.decomp_coeffs();
    let co: &ClassOrder = order;

    // Single-threaded pools sweep serially; the result buffer is already
    // class-ordered, so the sweep writes it directly.
    let threads = rayon::current_num_threads();
    if threads == 1 {
        for (job, gating) in jobs.iter_mut().zip(&gateways) {
            // SAFETY: the job's result buffer is valid for `result_len`
            // elements (the struct's contract) and is exclusively ours here;
            // the serial path is exclusive per job, so shared slices over
            // the raw operand pointers are sound.
            let result = unsafe { std::slice::from_raw_parts_mut(job.result, job.result_len) };
            let (a_slice, b_slice) = unsafe {
                (
                    std::slice::from_raw_parts(job.a, job.a_len),
                    std::slice::from_raw_parts(job.b, job.b_len),
                )
            };
            LieSeries::sweep_words_serial(
                a_series,
                a_slice,
                b_slice,
                gating,
                result,
                co.entries_cls(),
                co.decomp_cls(),
                &mut NoSink,
            );
        }
        return;
    }

    // SAFETY: each job's result is a distinct buffer (disjoint across jobs
    // by the caller's contract), and within a job each word class is
    // written by exactly one task (bundles never split a class), so the
    // concurrent single-writer stores never alias.
    let writers: Vec<RawResult<U>> = jobs
        .iter()
        .map(|j| RawResult {
            ptr: j.result,
            _marker: PhantomData,
        })
        .collect();

    // Flatten (job, unit) pairs into bundles — same cut rule as the
    // public kernel: pair-count balanced, cuts only between units (always
    // safe: unit word sets never overlap).
    let mut bundles: Vec<Vec<(usize, u32)>> = vec![Vec::new()];
    let mut cur = 0usize;
    let bundle_target = (total / (2 * threads)).clamp(MIN_BUNDLE_ENTRIES, BUNDLE_TARGET_ENTRIES);
    for (ji, gating) in gateways.iter().enumerate() {
        for (ui, unit) in gating.units.iter().enumerate() {
            let pairs = unit.pairs as usize;
            if cur >= bundle_target {
                bundles.push(Vec::new());
                cur = 0;
            }
            bundles.last_mut().unwrap().push((ji, ui as u32));
            cur += pairs;
        }
    }

    // The class-ordered result buffer needs no epilogue: the sweep's writes
    // are already final.
    sweep_bundles_parallel::<ClassInternalLayout, T, U>(
        jobs,
        &writers,
        &gateways,
        entries_slice,
        co.decomp_cls(),
        Some(co),
        decomp_coeffs_slice,
        &bundles,
    );
}

#[cfg(test)]
mod tests {
    use super::*;
    use ordered_float::NotNan;

    /// Bundling decides only which thread runs which anagram unit; the
    /// per-word accumulation order is unchanged, so the result must be
    /// bit-identical for any thread count. Guards the work-balanced bundle
    /// builder: the parallel sweep must never reorder, duplicate, or drop
    /// an accumulation.
    #[test]
    fn commutator_is_bit_identical_across_thread_counts() {
        use lyndon_rs::lyndon::{LyndonBasis, Sort};

        for (d, m) in [(2usize, 12usize), (3, 8), (4, 6)] {
            let words: Vec<LyndonWord<u8>> =
                LyndonBasis::<u8>::new(d, Sort::Lexicographical).generate_basis(m);
            let a_coefficients: Vec<_> = (0..words.len())
                .map(|i: usize| NotNan::new(((i % 11 + 1) as f64) / 7.0 - 0.9).unwrap())
                .collect();
            let b_coefficients: Vec<_> = (0..words.len())
                .map(|i: usize| NotNan::new(((i * 5 + 3) % 13) as f64 / 6.0 - 1.3).unwrap())
                .collect();
            let a: LieSeries<u8, NotNan<f64>> = LieSeries::new(words.clone(), a_coefficients);
            let b: LieSeries<u8, NotNan<f64>> = LieSeries::new(words, b_coefficients);

            let serial = {
                let pool = rayon::ThreadPoolBuilder::new()
                    .num_threads(1)
                    .build()
                    .expect("1-thread pool");
                let mut out = vec![NotNan::<f64>::default(); a.coefficients.len()];
                pool.install(|| {
                    LieSeries::commutator_coefficients(
                        &a,
                        &a.coefficients,
                        &b.coefficients,
                        &mut out,
                    )
                });
                out
            };
            for threads in [2usize, 4usize, 8usize] {
                let pool = rayon::ThreadPoolBuilder::new()
                    .num_threads(threads)
                    .build()
                    .expect("thread pool");
                let mut out = vec![NotNan::<f64>::default(); a.coefficients.len()];
                pool.install(|| {
                    LieSeries::commutator_coefficients(
                        &a,
                        &a.coefficients,
                        &b.coefficients,
                        &mut out,
                    )
                });
                assert_eq!(
                    serial, out,
                    "parallel result differs from serial for d={d}, m={m}, threads={threads}"
                );
            }
        }
    }
}
