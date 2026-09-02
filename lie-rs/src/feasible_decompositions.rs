//! Feasible-pair decomposition table: a compact, CSR-style replacement for
//! the O(D²) `commutator_basis_map` / `commutator_basis_map_coefficients`
//! tables.
//!
//! Only *canonical* pairs `(i, j)` with `i < j` and
//! `deg(i) + deg(j) <= max_degree` are stored. Entries are grouped into
//! *anagram units*: all pairs whose brackets land in the same letter
//! multiset are contiguous, sorted by `(target degree, content class)`.
//! Two different units never scatter onto the same basis word, so units
//! are conflict-free write regions for a parallel kernel.
//!
//! Decomposition indices are stored relative to the start of the
//! target-degree slice, with an absolute copy for [`FeasibleDecompositions::get`].

use std::sync::{Arc, OnceLock};

/// A canonical pair `(i, j)` with `i < j`. Its decomposition slice is
/// `decomp_start ..` up to the next flat entry's `decomp_start`.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Entry {
    pub(crate) i: u32,
    pub(crate) j: u32,
    pub(crate) decomp_start: u32,
}

/// A contiguous run of canonical pairs whose brackets land in the same
/// letter multiset: the unit's decompositions only hit degree-`target`
/// words of that content, so two units never write the same basis word.
#[derive(Clone, Debug)]
pub(crate) struct UnitMeta {
    /// Bracket degree `p + q` (equal to `|gamma|`).
    pub(crate) target: u8,
    /// The unit's letter multiset (counts per alphabet letter). Two units
    /// with the same `target` differ in `gamma` and vice versa.
    pub(crate) gamma: Vec<u8>,
    /// Bit `p` set iff the unit contains pairs whose smaller index has
    /// degree `p`.
    pub(crate) p_mask: [u64; 2],
    /// Entry range `entries[start .. end]`, sorted by `(p, q, i, j)`.
    /// The decomposition slice of entry `k` is
    /// `entries[k].decomp_start .. entries[k + 1].decomp_start`.
    pub(crate) start: u32,
    pub(crate) end: u32,
}

/// Packed active-entry ticket (see [`ActiveSegment::ticket_start`]): bits
/// 0..30 hold the entry's flat table index, bit 31 the `p_active`
/// orientation (`a = i, b = j`) and bit 30 the mirrored `q_active`. Table
/// sizes are debug-asserted below 2^30 entries at gating build time, so the
/// index never reaches the flag bits.
pub(crate) const TICKET_INDEX_MASK: u32 = 0x3fff_ffff;
pub(crate) const TICKET_P_ACTIVE: u32 = 1 << 31;
pub(crate) const TICKET_Q_ACTIVE: u32 = 1 << 30;

/// One maximal active run of a unit's entries for the commutation kernel's
/// sweep: a contiguous entry span gated on by degree support, scattering onto
/// the degree-`target` result slice.
#[derive(Clone, Copy)]
pub(crate) struct ActiveSegment {
    /// The run's entry span in the flat entry array, *including* each
    /// entry's flat successor: `entries[span_start .. span_end]` zipped with
    /// itself shifted by one gives entry `k`'s decomposition range
    /// `entries[k].decomp_start .. entries[k + 1].decomp_start`.
    pub(crate) span_start: u32,
    pub(crate) span_end: u32,
    /// The result scatter region `[rs, re)` — the degree-`target` slice.
    /// Decomposition indices are relative to `rs`.
    pub(crate) rs: u32,
    pub(crate) re: u32,
    /// Whether this is the last active segment of its unit. Segments of one
    /// unit scatter onto the same target words, so the parallel sweep must
    /// keep them in one bundle; it may only cut a bundle after this flag.
    /// (Per-entry orientation now rides in the tickets' flag bits — the
    /// segment-level `o1`/`o2` gating gates only which tickets exist.)
    pub(crate) last_of_unit: bool,
    /// The run's ticket range in the gating's flat ticket list: only the
    /// span's entries that passed the owning job's per-entry presence
    /// tests, in table order, with the resolved orientation flags packed
    /// in the top two bits ([`TICKET_INDEX_MASK`], [`TICKET_P_ACTIVE`],
    /// [`TICKET_Q_ACTIVE`]). The sweeps iterate this range instead of
    /// re-testing every table entry; built in entry order, the ticket list
    /// is a subsequence of the span, so per-word float summation order is
    /// unchanged.
    pub(crate) ticket_start: u32,
    pub(crate) ticket_end: u32,
    /// The left operand degree `p = deg(basis[i])` of this run — constant
    /// within the run (segments are maximal equal-`p` entry runs). The
    /// batch sweep hoists its compact-buffer address shifts per segment
    /// from this and [`Self::td`].
    pub(crate) p: u8,
    /// The unit's target (bracket) degree: the scatter words of this
    /// segment's decompositions all live in `degree_range(td)`.
    pub(crate) td: u8,
}

/// Coefficient type is generic; decomposition indices are stored as `u32`.
#[derive(Clone)]
pub(crate) struct FeasibleDecompositions<U> {
    /// Per basis word: its Lyndon degree. The basis is degree-grouped, so
    /// this array is non-decreasing.
    index_degrees: Vec<u8>,
    /// Per basis word: its letter content (counts per alphabet letter).
    index_contents: Vec<Vec<u8>>,
    /// For each degree `d`, the basis index where degree-`d` words start.
    /// `degree_range(target)` slices the result vector into the per-degree
    /// write regions used by the commutation kernel.
    degree_starts: Vec<u32>,
    /// Anagram units, sorted by `(target, content class)`; entry ranges
    /// are `entries[start .. end]`.
    units: Vec<UnitMeta>,
    /// Per-block contiguous entries, sorted by `(i, j)` within a block.
    entries: Vec<Entry>,
    /// Flat decomposition: basis indices relative to the owning unit's
    /// target-degree slice (kernel scatter path).
    decomp_indices: Vec<u32>,
    /// Flat decomposition: absolute basis indices (`get` path).
    decomp_indices_abs: Vec<u32>,
    /// Flat decomposition: parallel non-zero coefficients.
    decomp_coefficients: Vec<U>,
    /// Number of real entries (the trailing sentinel entry is excluded).
    num_entries: usize,
    /// Gating for full-support operands: one active segment per unit (every
    /// p-run of every unit is active, both orientations on), so a kernel
    /// call with `a_nonzero.len() == full-support cutoff` skips the gating
    /// walk entirely.
    full_support_segments: Arc<[ActiveSegment]>,
    /// The full-support segments' flat ticket list: presence is all-ones,
    /// so every table entry is active in both orientations and the ticket
    /// list is the entry stream itself, in segment order (each segment's
    /// ticket range mirrors its span minus the trailing successor).
    full_support_tickets: Arc<[u32]>,
    /// Within-degree class-contiguous scatter layout, populated at build
    /// time only when a degree slice outgrows L1 (see
    /// [`CLASS_ORDER_MIN_SLICE_WORDS`]) and otherwise created on demand by
    /// [`FeasibleDecompositions::ensure_class_order`]. Held behind an
    /// `Arc`: every kernel call and every series over the same basis
    /// reuses one ordering.
    class_order: OnceLock<Arc<ClassOrder>>,
}

/// Degree-slice size above which the commutation kernel's scatter writes
/// leave L1 and the class-contiguous layout pays for its permutation
/// epilogue. Below it the public-order scatter is already cache-resident.
const CLASS_ORDER_MIN_SLICE_WORDS: usize = 4096;

/// Within-degree class-contiguous scatter layout: the basis words of each
/// degree are regrouped so that words of the same letter content (anagram
/// class) are contiguous, with degree-slice boundaries preserved. A unit's
/// decompositions only ever hit its own class (content homogeneity), so in
/// this layout every unit's stores are dense and consecutive units write
/// consecutive blocks.
///
/// The ordering depends only on the basis — every series over the same
/// basis rearranges coefficients identically — so one handle (cheap to
/// clone: `Arc` internally) amortizes across operand series, across every
/// kernel call of a fold, and across batches of folds. Obtain it via
/// [`ClassOrderedCommutation::class_order`](crate::ClassOrderedCommutation::class_order).
#[derive(Clone)]
pub struct ClassOrder {
    /// Internal (class-contiguous) position -> public basis index.
    perm: Vec<u32>,
    /// Public basis index -> internal (class-contiguous) position. The
    /// epilogue walks this in public order so its result writes are
    /// sequential (full cache lines) while its scratch reads — internal
    /// positions, just written by the sweep — stay cache-hot.
    inv: Vec<u32>,
    /// Decomposition scatter indices relative to the owning unit's degree
    /// slice, expressed in the internal layout: a class-mode sweep writes
    /// `scratch[rs + rel]` with `rel` dense over the unit's class block.
    decomp_cls: Vec<u32>,
    /// Entry table with `i`/`j` relabeled to internal positions, for
    /// sweeps whose operands are class-ordered. Same order and
    /// `decomp_start`s as the public table.
    entries_cls: Vec<Entry>,
    /// Lyndon degrees indexed by internal position.
    degree_cls: Vec<u8>,
}

impl ClassOrder {
    /// Internal (class-contiguous) position -> public basis index: the
    /// final-epilogue map. `public[k] = class[inv()[k]]` translates a
    /// class-ordered slice back to public order with sequential writes.
    #[inline]
    pub fn perm(&self) -> &[u32] {
        &self.perm
    }

    /// Public basis index -> internal (class-contiguous) position: the
    /// input-conversion map. `class[inv()[k]] = public[k]` gathers a
    /// public-order slice into class order.
    #[inline]
    pub fn inv(&self) -> &[u32] {
        &self.inv
    }

    /// Decomposition scatter indices in the internal layout.
    #[inline]
    pub(crate) fn decomp_cls(&self) -> &[u32] {
        &self.decomp_cls
    }

    /// Entry table with internal-position endpoints.
    #[inline]
    pub(crate) fn entries_cls(&self) -> &[Entry] {
        &self.entries_cls
    }

    /// Lyndon degrees indexed by internal (class) position. The positions
    /// are degree-grouped, so this slice is non-decreasing.
    #[inline]
    pub fn degree_cls(&self) -> &[u8] {
        &self.degree_cls
    }
}

/// Incremental builder used during `LieSeries::new`: canonical pairs arrive
/// row-major (`i` ascending, `j > i` ascending within a row) as the structure
/// constants are computed.
pub(crate) struct Builder<U> {
    /// Per row `i`: `(j, decomposition indices, decomposition coefficients)`.
    rows: Vec<Vec<(u32, Vec<u32>, Vec<U>)>>,
    /// Letter content of each basis word (counts per alphabet letter), in
    /// basis order.
    contents: Vec<Vec<u8>>,
    /// Degrees of the basis words (content sums), in basis order
    /// (non-decreasing).
    degrees: Vec<u8>,
}

impl<U: Clone> Builder<U> {
    /// `contents[i]` is the letter multiset of basis word `i`; the basis
    /// must be degree-grouped (degrees non-decreasing along basis order).
    pub(crate) fn new(contents: &[Vec<u8>]) -> Self {
        let degrees: Vec<u8> = contents
            .iter()
            .map(|c| c.iter().copied().sum::<u8>())
            .collect();
        Self {
            rows: (0..contents.len()).map(|_| Vec::new()).collect(),
            contents: contents.to_vec(),
            degrees,
        }
    }

    /// Records the decomposition of the canonical pair `(i, j)`. Pairs
    /// must arrive row-major; decompositions must only hit
    /// degree-`deg(i) + deg(j)` basis words.
    pub(crate) fn push(&mut self, i: usize, j: usize, indices: &[usize], coefficients: &[U]) {
        self.rows[i].push((
            u32::try_from(j).expect("basis index exceeds u32"),
            indices
                .iter()
                .map(|&x| u32::try_from(x).expect("basis index exceeds u32"))
                .collect(),
            coefficients.to_vec(),
        ));
    }

    pub(crate) fn finish(self) -> FeasibleDecompositions<U> {
        let degrees = &self.degrees;
        assert!(
            degrees.len() <= u32::MAX as usize,
            "basis size exceeds u32 indices"
        );

        // Degree starts: for each degree, the first basis index carrying it.
        // Absent degrees point at the next present degree (their slice is
        // empty). The basis is degree-grouped, so a single forward scan
        // suffices.
        let max_degree = degrees.last().copied().unwrap_or(0) as usize;
        let mut degree_starts = vec![0u32; max_degree + 2];
        let mut d = 0usize;
        for (index, &g) in degrees.iter().enumerate() {
            while d < g as usize {
                d += 1;
                degree_starts[d] = index as u32;
            }
        }
        degree_starts[max_degree + 1] = degrees.len() as u32;

        // Flatten and sort by (target, content class, p, q, i, j): anagram
        // units in kernel iteration order, entries sorted within their
        // unit.
        let contents = &self.contents;
        let mut flat: Vec<(u8, Vec<u8>, u8, u8, u32, u32, Vec<u32>, Vec<U>)> = Vec::new();
        for (i, row) in self.rows.iter().enumerate() {
            for (j, indices, coefficients) in row {
                let (p, q) = (degrees[i], degrees[*j as usize]);
                debug_assert!(
                    p <= q,
                    "non-degree-grouped basis: i < j but deg(i) > deg(j)"
                );
                let mut gamma = contents[i].clone();
                for k in 0..gamma.len() {
                    gamma[k] += contents[*j as usize][k];
                }
                flat.push((
                    p + q,
                    gamma,
                    p,
                    q,
                    i as u32,
                    *j,
                    indices.clone(),
                    coefficients.clone(),
                ));
            }
        }
        flat.sort_unstable_by(|a, b| {
            (a.0, &a.1, a.2, a.3, a.4, a.5).cmp(&(b.0, &b.1, b.2, b.3, b.4, b.5))
        });

        let total: usize = flat
            .iter()
            .map(|(_, _, _, _, _, _, idx, _)| idx.len())
            .sum();
        let mut units: Vec<UnitMeta> = Vec::new();
        let mut entries = Vec::with_capacity(flat.len());
        let mut decomp_indices = Vec::with_capacity(total);
        let mut decomp_indices_abs = Vec::with_capacity(total);
        let mut decomp_coefficients = Vec::with_capacity(total);

        for (target, gamma, p, _q, i, j, indices, coefficients) in flat {
            debug_assert!(
                target as usize <= max_degree,
                "infeasible pair ({i}, {j}) pushed: target degree {target} > {max_degree}"
            );
            let target_start = degree_starts[target as usize];
            if units
                .last()
                .is_none_or(|u| (u.target, &u.gamma) != (target, &gamma))
            {
                units.push(UnitMeta {
                    target,
                    gamma: gamma.clone(),
                    p_mask: [0, 0],
                    start: entries.len() as u32,
                    end: 0,
                });
            }
            let unit = units.last_mut().unwrap();
            let mask_word = &mut unit.p_mask[(p / 64) as usize];
            *mask_word |= 1u64 << (p % 64);
            debug_assert!(
                indices.iter().all(|&x| contents[x as usize] == gamma),
                "decomposition of pair ({i}, {j}) is not content-homogeneous"
            );
            debug_assert!(
                indices.iter().all(|&x| degrees[x as usize] == target),
                "decomposition of pair ({i}, {j}) is not degree-homogeneous"
            );
            entries.push(Entry {
                i,
                j,
                decomp_start: decomp_indices.len() as u32,
            });
            for &x in &indices {
                debug_assert!(
                    x >= target_start,
                    "decomposition index below target-degree slice start"
                );
                decomp_indices.push(x - target_start);
            }
            decomp_indices_abs.extend_from_slice(&indices);
            decomp_coefficients.extend_from_slice(&coefficients);
        }
        // Close the units' entry ranges (the sentinel entry appended below
        // bounds the last real entry's decomposition slice).
        let mut end = entries.len() as u32;
        for unit in units.iter_mut().rev() {
            unit.end = end;
            end = unit.start;
        }
        let total = decomp_indices.len() as u32;

        let num_entries = entries.len();
        entries.push(Entry {
            i: 0,
            j: 0,
            decomp_start: total,
        });
        // Class-contiguous scatter layout: only worth building when some
        // degree slice outgrows L1 — below that the public-order scatter is
        // already cache-resident and the layout's permutation epilogue
        // would only add work (12×2's whole kernel is smaller than one
        // epilogue pass).
        let max_slice = (0..=max_degree)
            .map(|d| degree_starts[d + 1] - degree_starts[d])
            .max()
            .unwrap_or(0);
        let class_order = OnceLock::new();
        if max_slice as usize >= CLASS_ORDER_MIN_SLICE_WORDS {
            let _ = class_order.set(Arc::new(build_class_order(
                degrees,
                contents,
                &degree_starts,
                &units,
                &entries,
                &decomp_indices,
            )));
        }

        // Full-support gating: every run of every unit is active with both
        // orientations. Segments are maximal equal-`p` entry runs (the same
        // shape the degree-gated prologue produces), so a segment's left
        // operand degree — and with it the batch sweep's hoisted compact
        // address shifts — is constant within the segment. Presence is
        // all-ones, so the flat ticket list is the entry stream itself
        // (every entry active, both flag bits packed; see `TICKET_*`).
        debug_assert!(
            num_entries < (1 << 30),
            "entry indices must fit a ticket's 30 bits"
        );
        let mut full_support_segments: Vec<ActiveSegment> = Vec::new();
        let mut full_support_tickets: Vec<u32> = Vec::new();
        for unit in units.iter() {
            let td = unit.target;
            let (rs, re) = (
                degree_starts[td as usize],
                degree_starts[td as usize + 1],
            );
            let mut run_start = unit.start;
            let mut cur_p = u8::MAX;
            for ei in unit.start..unit.end {
                let p = degrees[entries[ei as usize].i as usize];
                if p == cur_p {
                    continue;
                }
                if cur_p != u8::MAX {
                    let ticket_start = full_support_tickets.len() as u32;
                    for e in run_start..ei {
                        full_support_tickets.push(e | TICKET_P_ACTIVE | TICKET_Q_ACTIVE);
                    }
                    full_support_segments.push(ActiveSegment {
                        span_start: run_start,
                        span_end: ei + 1,
                        ticket_start,
                        ticket_end: full_support_tickets.len() as u32,
                        rs,
                        re,
                        last_of_unit: false,
                        p: cur_p,
                        td,
                    });
                }
                cur_p = p;
                run_start = ei;
            }
            if cur_p != u8::MAX {
                let ticket_start = full_support_tickets.len() as u32;
                for e in run_start..unit.end {
                    full_support_tickets.push(e | TICKET_P_ACTIVE | TICKET_Q_ACTIVE);
                }
                full_support_segments.push(ActiveSegment {
                    span_start: run_start,
                    span_end: unit.end + 1,
                    ticket_start,
                    ticket_end: full_support_tickets.len() as u32,
                    rs,
                    re,
                    last_of_unit: true,
                    p: cur_p,
                    td,
                });
            }
        }

        FeasibleDecompositions {
            index_degrees: degrees.clone(),
            index_contents: contents.clone(),
            degree_starts,
            units,
            entries,
            decomp_indices,
            decomp_indices_abs,
            decomp_coefficients,
            num_entries,
            full_support_segments: Arc::from(full_support_segments),
            full_support_tickets: Arc::from(full_support_tickets),
            class_order,
        }
    }
}

/// Builds the within-degree class-contiguous permutation and the
/// re-indexed decomposition scatter array.
///
/// Within each degree, words are grouped by letter content (classes ordered
/// by content bytes, words in public order inside a class); degree-slice
/// boundaries are preserved, so a degree slice's start offset is valid in
/// both layouts and only the positions inside the slice move.
fn build_class_order(
    degrees: &[u8],
    contents: &[Vec<u8>],
    degree_starts: &[u32],
    units: &[UnitMeta],
    entries: &[Entry],
    decomp_rel: &[u32],
) -> ClassOrder {
    let len = degrees.len();
    let mut inv: Vec<u32> = (0..len as u32).collect();
    let max_degree = degree_starts.len() - 2;
    for d in 0..=max_degree {
        let lo = degree_starts[d] as usize;
        let hi = degree_starts[d + 1] as usize;
        if hi <= lo + 1 {
            continue;
        }
        // Stable sort: equal contents (one class) keep their public order.
        let mut order: Vec<usize> = (lo..hi).collect();
        order.sort_by(|&a, &b| contents[a].cmp(&contents[b]));
        for (pos, &w) in order.iter().enumerate() {
            inv[w] = (lo + pos) as u32;
        }
    }
    let mut perm: Vec<u32> = vec![0; len];
    for (w, &p) in inv.iter().enumerate() {
        perm[p as usize] = w as u32;
    }
    let inv = inv;
    // Re-index the scatter array: `rel_cls = inv[rs + rel] - rs` per unit.
    // The entries' sentinel-closed successor chains give each entry's
    // decomposition range, exactly as the sweep reads it.
    let mut decomp_cls = vec![0u32; decomp_rel.len()];
    for unit in units {
        let rs = degree_starts[unit.target as usize] as usize;
        for ei in unit.start..unit.end {
            let from = entries[ei as usize].decomp_start as usize;
            let to = entries[ei as usize + 1].decomp_start as usize;
            for k in from..to {
                decomp_cls[k] = inv[rs + decomp_rel[k] as usize] - rs as u32;
            }
        }
    }
    // Entry table relabeled to internal positions (same order, same
    // decomp_starts; the sentinel's endpoints are never read).
    let entries_cls: Vec<Entry> = entries
        .iter()
        .map(|e| Entry {
            i: inv[e.i as usize],
            j: inv[e.j as usize],
            decomp_start: e.decomp_start,
        })
        .collect();
    let mut degree_cls = vec![0u8; degrees.len()];
    for (w, &p) in inv.iter().enumerate() {
        degree_cls[p as usize] = degrees[w];
    }
    ClassOrder {
        perm,
        inv,
        decomp_cls,
        entries_cls,
        degree_cls,
    }
}

impl<U> FeasibleDecompositions<U> {
    /// Total number of stored feasible pairs.
    pub(crate) fn len(&self) -> usize {
        self.num_entries
    }

    /// The class-contiguous layout when it was prebuilt at table
    /// construction (degree slice above the L1 threshold); `None` means
    /// the default kernel paths run direct.
    #[inline]
    pub(crate) fn cached_class_order(&self) -> Option<&Arc<ClassOrder>> {
        self.class_order.get()
    }

    /// The class-contiguous layout, building it on first request when the
    /// table did not prebuild one. Cheap after the first call, and shared
    /// as an `Arc` across every consumer of this table.
    #[inline]
    pub(crate) fn ensure_class_order(&self) -> Arc<ClassOrder> {
        self.class_order
            .get_or_init(|| {
                let units = self.units.clone();
                let entries = self.entries.clone();
                let decomp = self.decomp_indices.clone();
                Arc::new(build_class_order(
                    &self.index_degrees,
                    &self.index_contents,
                    &self.degree_starts,
                    &units,
                    &entries,
                    &decomp,
                ))
            })
            .clone()
    }

    /// Test hook: builds the class-contiguous layout regardless of the
    /// slice-size threshold, so correctness tests can exercise the class
    /// sweep on small (fast-to-build) shapes.
    #[cfg(test)]
    pub(crate) fn force_class_order(&mut self) {
        let units = self.units.clone();
        let entries = self.entries.clone();
        let decomp = self.decomp_indices.clone();
        let built = Arc::new(build_class_order(
            &self.index_degrees,
            &self.index_contents,
            &self.degree_starts,
            &units,
            &entries,
            &decomp,
        ));
        let _ = self.class_order.set(built);
    }

    /// Test hook: drops the class layout so a threshold-qualifying table
    /// can serve as the direct-layout reference.
    #[cfg(test)]
    pub(crate) fn clear_class_order(&mut self) {
        let _ = self.class_order.take();
    }

    /// The maximum basis degree.
    #[cfg_attr(not(test), allow(dead_code))]
    pub(crate) fn max_degree(&self) -> usize {
        self.degree_starts.len() - 2
    }

    /// The Lyndon degree of basis word `i` (O(1)).
    #[inline]
    pub(crate) fn degree_of(&self, i: usize) -> usize {
        self.index_degrees[i] as usize
    }

    /// Anagram units in kernel iteration order.
    #[inline]
    pub(crate) fn units(&self) -> &[UnitMeta] {
        &self.units
    }

    /// The flat entry array in storage order.
    #[inline]
    pub(crate) fn entries(&self) -> &[Entry] {
        &self.entries
    }

    /// A unit's entry slice, extended by one: `slice[k + 1]` is the flat
    /// successor of `slice[k]`, so entry `k`'s decomposition range is
    /// `slice[k].decomp_start .. slice[k + 1].decomp_start`.
    #[cfg_attr(not(test), allow(dead_code))]
    #[inline]
    pub(crate) fn entry_span(&self, unit: &UnitMeta) -> &[Entry] {
        &self.entries[unit.start as usize..=unit.end as usize]
    }

    /// Active segments for full-support operands (every run of every unit
    /// active, both orientations on) plus their flat ticket list: shared by
    /// all full-support kernel calls, so the gating walk is skipped
    /// entirely. Presence is all-ones, so the tickets cover every table
    /// entry with both flag bits packed, and each segment's ticket range
    /// mirrors its span minus the trailing successor.
    #[inline]
    pub(crate) fn full_support_gating(&self) -> (Arc<[ActiveSegment]>, Arc<[u32]>) {
        (
            Arc::clone(&self.full_support_segments),
            Arc::clone(&self.full_support_tickets),
        )
    }

    /// The flat relative decomposition index array (kernel scatter path).
    #[inline]
    pub(crate) fn decomp_indices_rel(&self) -> &[u32] {
        &self.decomp_indices
    }

    /// The flat decomposition coefficient array.
    #[inline]
    pub(crate) fn decomp_coeffs(&self) -> &[U] {
        &self.decomp_coefficients
    }

    /// A unit's entries as `(entry, decomposition indices, decomposition
    /// coefficients)`. Decomposition indices are relative to the start of
    /// the unit's target-degree slice.
    #[cfg(test)]
    pub(crate) fn unit_iter<'a>(
        &'a self,
        unit: &'a UnitMeta,
    ) -> impl Iterator<Item = (&'a Entry, &'a [u32], &'a [U])> + 'a {
        let span = self.entry_span(unit);
        span[..span.len() - 1]
            .iter()
            .zip(span[1..].iter())
            .map(move |(e, next)| {
                let from = e.decomp_start as usize;
                let to = next.decomp_start as usize;
                (
                    e,
                    &self.decomp_indices[from..to],
                    &self.decomp_coefficients[from..to],
                )
            })
    }

    /// The basis index where degree-`degree` words start (see
    /// [`FeasibleDecompositions::degree_range`]).
    #[inline]
    pub(crate) fn degree_start(&self, degree: usize) -> usize {
        self.degree_starts[degree] as usize
    }

    /// The result-vector range `[start, end)` holding the degree-`target`
    /// basis words.
    #[inline]
    /// DEBUG/telemetry: per-position Lyndon degree + degree-slice starts.
    #[doc(hidden)]
    pub fn debug_degree_layout(&self) -> (Vec<u8>, Vec<u32>) {
        (self.index_degrees.clone(), self.degree_starts.clone())
    }

    pub(crate) fn degree_range(&self, target: u8) -> (usize, usize) {
        let t = target as usize;
        (
            self.degree_starts[t] as usize,
            self.degree_starts[t + 1] as usize,
        )
    }

    /// All entries in storage order as
    /// `(i, j, absolute decomposition indices, decomposition coefficients)`.
    #[cfg(test)]
    pub(crate) fn iter_entries(&self) -> impl Iterator<Item = (usize, usize, &[u32], &[U])> + '_ {
        self.units.iter().flat_map(move |unit| {
            let span = self.entry_span(unit);
            span[..span.len() - 1]
                .iter()
                .zip(span[1..].iter())
                .map(move |(e, next)| {
                    let from = e.decomp_start as usize;
                    let to = next.decomp_start as usize;
                    (
                        e.i as usize,
                        e.j as usize,
                        &self.decomp_indices_abs[from..to],
                        &self.decomp_coefficients[from..to],
                    )
                })
        })
    }

    /// The decomposition of the bracket `[basis[i], basis[j]]`, if stored.
    ///
    /// Querying with `j < i` returns the canonical data with
    /// `swapped = true`; the caller must negate the coefficients. `i == j`
    /// is never stored.
    pub(crate) fn get(&self, i: usize, j: usize) -> Option<(&[u32], &[U], bool)> {
        let (min, max, swapped) = match i.cmp(&j) {
            std::cmp::Ordering::Less => (i, j, false),
            std::cmp::Ordering::Greater => (j, i, true),
            std::cmp::Ordering::Equal => return None,
        };
        let (min, max) = (u32::try_from(min).ok()?, u32::try_from(max).ok()?);
        let (c_min, c_max) = (
            &self.index_contents[min as usize],
            &self.index_contents[max as usize],
        );
        let mut gamma = c_min.clone();
        for (g, c) in gamma.iter_mut().zip(c_max) {
            *g += c;
        }
        let target = gamma.iter().copied().sum::<u8>();
        // Units are sorted by (target, content class) and the class
        // determines the target, so the unit is found with one binary
        // search. The target is widened: arbitrary (infeasible) queries may
        // overflow `u8`.
        let pos = self
            .units
            .partition_point(|u| (u.target as u16, &u.gamma) < (target as u16, &gamma));
        let unit = self.units.get(pos)?;
        if unit.gamma != gamma {
            return None; // no unit for this content pair: infeasible pair
        }
        let entries = &self.entries[unit.start as usize..unit.end as usize];
        let pos = entries.partition_point(|e| (e.i, e.j) < (min, max));
        let entry = entries.get(pos)?;
        if (entry.i, entry.j) != (min, max) {
            return None;
        }
        let from = entry.decomp_start as usize;
        // The flat successor bounds the slice (next unit's first entry or
        // the trailing sentinel; both carry the right `decomp_start`).
        let next = self.entries[unit.start as usize + pos + 1].decomp_start as usize;
        Some((
            &self.decomp_indices_abs[from..next],
            &self.decomp_coefficients[from..next],
            swapped,
        ))
    }
}

#[cfg(test)]
mod test {
    use super::*;

    /// INVARIANT CHECK: within one degree-target slice, two different
    /// units must never scatter onto the same basis word — the batch's
    /// packs cut at unit boundaries and rely on per-word single-writer
    /// accumulation for bit-identical results.
    #[test]
    fn debug_unit_word_ownership_disjoint() {
        // Real 3-letter basis through degree 5: word contents from the
        // Lyndon enumeration, in basis order.
        let contents: Vec<Vec<u8>> = vec![
            vec![1, 0, 0], // a
            vec![0, 1, 0], // b
            vec![0, 0, 1], // c
            vec![1, 1, 0], // ab
            vec![1, 0, 1], // ac
            vec![0, 1, 1], // bc
            vec![2, 1, 0], // aab
            vec![1, 2, 0], // abb
            vec![2, 0, 1], // aac
            vec![1, 0, 2], // acc
            vec![0, 2, 1], // bbc
            vec![0, 1, 2], // bcc
            vec![1, 1, 1], // abc
            vec![3, 1, 0], // aaab
            vec![1, 3, 0], // abbb
            vec![1, 1, 2], // abcc
            vec![2, 2, 1], // aabb
        ];
        let mut b = Builder::<f64>::new(&contents);
        // Push every canonical pair with a decomposition into every word of
        // the pair's content class at the target degree (row-major, as the
        // kernel's table build does).
        let n = contents.len();
        for i in 0..n {
            for j in (i + 1)..n {
                let mut gamma = vec![0u8; 3];
                for k in 0..3 {
                    gamma[k] = contents[i][k] + contents[j][k];
                }
                let _target: u8 = gamma.iter().sum();
                // all words of degree `target` with this content
                let targets: Vec<usize> = (0..n).filter(|w| contents[*w] == gamma).collect();
                if targets.is_empty() {
                    continue; // no word of this content at this degree
                }
                b.push(i, j, &targets, &vec![1.0; targets.len()]);
            }
        }
        let t = b.finish();

        // Walk every unit; collect (target, rel) -> unit index.
        let mut owner: std::collections::HashMap<(u8, u32), usize> =
            std::collections::HashMap::new();
        let mut violations = 0usize;
        for (uid, unit) in t.units().iter().enumerate() {
            let span = t.entry_span(unit);
            for (entry, next) in span[..span.len() - 1].iter().zip(span[1..].iter()) {
                let from = entry.decomp_start as usize;
                let to = next.decomp_start as usize;
                for &rel in &t.decomp_indices_rel()[from..to] {
                    let key = (unit.target, rel);
                    match owner.get(&key) {
                        Some(prev) if *prev != uid => {
                            violations += 1;
                            if violations <= 8 {
                                println!(
                                    "VIOLATION: degree {} rel {} claimed by units {} and {}",
                                    unit.target, rel, prev, uid
                                );
                            }
                        }
                        _ => {
                            owner.insert(key, uid);
                        }
                    }
                }
            }
        }
        println!(
            "units={} rels={} (target,rel) keys={} violations={}",
            t.units().len(),
            t.decomp_indices_rel().len(),
            owner.len(),
            violations
        );
        assert_eq!(
            violations, 0,
            "units share scatter words — pack-level bit-identical accumulation is unsound"
        );
    }

    /// Five fake basis words over a 2-letter alphabet: `a`, `b` (degree 1),
    /// `ab` (degree 2), `aab`, `abb` (degree 3). Degree-`t` decompositions
    /// may only hit degree-`t` words whose content is the multiset union of
    /// the pair's contents — e.g. `[a, b]` -> `ab`, `[a, ab]` -> `aab`.
    #[test]
    fn builder_round_trip() {
        let contents = vec![
            vec![1, 0], // 0: a
            vec![0, 1], // 1: b
            vec![1, 1], // 2: ab
            vec![2, 1], // 3: aab
            vec![1, 2], // 4: abb
        ];
        let mut b = Builder::<f64>::new(&contents);
        // Row-major pushes of canonical pairs, including an empty
        // decomposition and a row with no feasible pairs.
        // (0,1): a+b, target 2, class (1,1,0) -> word 2 (ab)
        b.push(0, 1, &[2], &[1.5]);
        // (0,2): a+ab, target 3, class (2,1,0) -> word 3 (aab)
        b.push(0, 2, &[3], &[0.25]);
        // (1,2): b+ab, target 3, class (1,2,0): word 4 (abb) exists but the
        // pair's decomposition is empty here; the unit still exists.
        b.push(1, 2, &[], &[]);
        let t = b.finish();

        assert_eq!(t.len(), 3);
        assert_eq!(t.max_degree(), 3);
        assert_eq!(t.degree_range(1), (0, 2));
        assert_eq!(t.degree_range(2), (2, 3));
        assert_eq!(t.degree_range(3), (3, 5));

        // Units sorted by (target, content class), entries sorted by (i, j).
        let units: Vec<_> = t
            .units()
            .iter()
            .map(|u| (u.target, u.gamma.clone()))
            .collect();
        assert_eq!(units, [(2, vec![1, 1]), (3, vec![1, 2]), (3, vec![2, 1])]);
        let unit12: &UnitMeta = t
            .units()
            .iter()
            .find(|u| u.target == 3 && u.gamma[0] == 1)
            .unwrap();
        let entries12: Vec<_> = t
            .unit_iter(unit12)
            .map(|(e, indices, coefficients)| (e.i, e.j, indices, coefficients))
            .collect();
        assert_eq!(entries12.len(), 1);
        assert_eq!((entries12[0].0, entries12[0].1), (1, 2));
        assert_eq!(entries12[0].2, &[] as &[u32]);
        assert_eq!(entries12[0].3, &[] as &[f64]);

        // Canonical lookups, including the mirrored (negated) query.
        assert_eq!(t.get(0, 1), Some((&[2u32][..], &[1.5][..], false)));
        assert_eq!(t.get(1, 0), Some((&[2u32][..], &[1.5][..], true)));
        assert_eq!(t.get(2, 1), Some((&[] as &[u32], &[] as &[f64], true)));
        assert_eq!(t.get(0, 0), None);
        assert_eq!(t.get(1, 1), None);
        assert_eq!(t.get(0, 3), None); // feasible degrees, but never pushed
    }

    #[test]
    fn empty_table() {
        let t = Builder::<f64>::new(&[vec![1, 0], vec![0, 1], vec![2, 0]]).finish();
        assert_eq!(t.len(), 0);
        assert_eq!(t.max_degree(), 2);
        assert_eq!(t.units().len(), 0);
        assert_eq!(t.get(0, 1), None);
    }

    /// The kernel scatters with indices relative to the target-degree slice;
    /// the absolute copy used by `get` must reconstruct the same indices.
    #[test]
    fn relative_indices_match_absolute() {
        let contents = vec![
            vec![1, 0], // 0: a
            vec![0, 1], // 1: b
            vec![1, 1], // 2: ab   (degree 2, class (1,1,0))
            vec![2, 1], // 3: aab  (degree 3, class (2,1,0))
        ];
        let mut b = Builder::<f64>::new(&contents);
        // target 2, class (1,1,0): word 2 is the only degree-2 word of that
        // content.
        b.push(0, 1, &[2], &[1.0]);
        // target 3, class (2,1,0): word 3 is the only degree-3 word of that
        // content.
        b.push(0, 2, &[3], &[-1.0]);
        let t = b.finish();

        for unit in t.units() {
            let (start, _) = t.degree_range(unit.target);
            for (e, rel, _) in t.unit_iter(unit) {
                let (abs, _, swapped) = t.get(e.i as usize, e.j as usize).expect("stored");
                assert!(!swapped);
                for (rel, abs) in rel.iter().zip(abs) {
                    assert_eq!(start + *rel as usize, *abs as usize);
                }
            }
        }
    }
}
