use super::*;

/// Memo for the commutation kernel's prologues: maps a job's operand
/// support lists to its resolved gating — the active per-word write classes
/// with their `(ticket, decomp position)` contribution lists, presence
/// results baked in.
///
/// The key is a 128-bit fingerprint per side of the *exact* support lists
/// (not just their degree masks): the cached tickets bake per-entry
/// presence results, which degree masks do not determine — two supports
/// can cover the same degrees through different words (exactly what the
/// old per-entry bit tests re-resolved on every kernel call). The
/// fingerprints mix the list lengths in and avalanche-finish both streams,
/// so distinct lists reuse an entry only with negligible probability; an
/// entry's gating is reused solely for identical supports.
///
/// In a log-signature fold the DAG's node support lists are a value-
/// independent fixed point, so steady-state folds (and whole batches)
/// hit the cache for every job: the per-entry resolution runs once per
/// distinct (job, operand supports), not per fold.
///
/// The cache is valid only for the decomposition table of the series whose
/// prologues populated it.
#[derive(Clone, Default)]
pub struct GatingCache {
    /// Open-addressed linear-probe table keyed by the two support
    /// fingerprints. Distinct pairs per configuration are few (the DAG's
    /// distinct node-support pairs), so a fixed-capacity table with a
    /// cheap multiplicative hash beats a full hash map: the lookup runs per
    /// kernel call, and at small grids the SipHash + bucket walk was a
    /// measurable share of the fold.
    slots: Vec<Slot>,
}

/// Key + value for one [`GatingCache`] slot. `None` key = empty; the all-
/// zero fingerprint pair is a legitimate key (two empty supports), so
/// emptiness is tracked out of band.
#[derive(Clone, Default)]
struct Slot {
    key: Option<([u64; 2], [u64; 2])>,
    value: (Arc<[UnitRange]>, Arc<[u32]>, Arc<[u32]>, usize),
}

impl GatingCache {
    /// Looks up the support fingerprint pair, returning the memoized
    /// `(active units, flat tickets, active word positions, pair count)`.
    #[inline]
    fn get(
        &self,
        key: ([u64; 2], [u64; 2]),
    ) -> Option<&(Arc<[UnitRange]>, Arc<[u32]>, Arc<[u32]>, usize)> {
        let cap = self.slots.len();
        if cap == 0 {
            return None;
        }
        let mut idx = Self::hash(&key) & (cap - 1);
        loop {
            let slot = &self.slots[idx];
            match &slot.key {
                Some(k) if *k == key => return Some(&slot.value),
                Some(_) => idx = (idx + 1) & (cap - 1),
                None => return None,
            }
        }
    }

    /// Inserts, growing the table to the next power of two when full.
    fn insert(
        &mut self,
        key: ([u64; 2], [u64; 2]),
        value: (Arc<[UnitRange]>, Arc<[u32]>, Arc<[u32]>, usize),
    ) {
        if self.slots.len() * 3 / 4 <= self.len() {
            let old = std::mem::take(&mut self.slots);
            self.slots = vec![Slot::default(); (old.len() * 2).max(8)];
            for slot in old {
                if let Some(k) = slot.key {
                    self.insert_unchecked(k, slot.value);
                }
            }
        }
        self.insert_unchecked(key, value);
    }

    fn insert_unchecked(
        &mut self,
        key: ([u64; 2], [u64; 2]),
        value: (Arc<[UnitRange]>, Arc<[u32]>, Arc<[u32]>, usize),
    ) {
        let cap = self.slots.len();
        let mut idx = Self::hash(&key) & (cap - 1);
        while self.slots[idx].key.is_some() {
            idx = (idx + 1) & (cap - 1);
        }
        self.slots[idx] = Slot {
            key: Some(key),
            value,
        };
    }

    fn len(&self) -> usize {
        self.slots.iter().filter(|s| s.key.is_some()).count()
    }

    /// Multiplicative fold of the four mask words — cheaper than SipHash
    /// and collision-free enough for a table holding a handful of entries.
    #[inline]
    fn hash(key: &([u64; 2], [u64; 2])) -> usize {
        let mut h = 0xcbf29ce484222325u64;
        for w in [&key.0, &key.1] {
            for x in w {
                h = (h ^ x).wrapping_mul(0x100000001b3);
            }
        }
        h as usize
    }
}

impl<
    T: Clone + Ord + Generator + Hash + Eq,
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
> LieSeries<T, U>
{
    /// Indices of non-zero coefficients that the commutation kernel will
    /// actually use: non-zero values on basis elements below `max_degree`.
    pub(super) fn kernel_prologue(
        a_series: &LieSeries<T, U>,
        a_nonzero: &[usize],
        b_nonzero: &[usize],
    ) -> KernelGating {
        let mut cache = GatingCache::default();
        Self::kernel_prologue_cached(a_series, a_nonzero, b_nonzero, &mut cache)
    }
    pub(super) fn kernel_prologue_cached(
        a_series: &LieSeries<T, U>,
        a_nonzero: &[usize],
        b_nonzero: &[usize],
        cache: &mut GatingCache,
    ) -> KernelGating {
        let table = &a_series.feasible_decompositions;

        // Full-support fast path: the nonzero lists must be EXACTLY the
        // kernel-visible prefix `0..cutoff` — every index the sweep's
        // presence tests can ever reference (entries only pair positions
        // strictly below `max_degree`). Length alone is NOT sufficient:
        // node support lists recorded by a batch fold legitimately contain
        // degree-`max_degree` positions, and a same-length list that
        // includes them would skip the presence tests and read compact
        // slots the operand layout does not cover (wrong values — see the
        // pooled-dag leaf-reuse regression test).
        let cutoff = table.degree_start(table.max_degree());
        let full_prefix = |s: &[usize]| s.len() == cutoff && s.iter().all(|&p| p < cutoff);
        if full_prefix(a_nonzero) && full_prefix(b_nonzero) {
            let (units, tickets, unit_words) = table.full_support_gating_public();
            return KernelGating {
                total_pairs: units.iter().map(|u| u.pairs as usize).sum(),
                units,
                tickets,
                unit_words,
            };
        }

        // Memoized gating keyed by the exact support lists (see
        // `GatingCache`): on a hit neither the bitsets, the per-entry
        // presence tests, nor the per-unit word-set collection run at all.
        let key = (
            Self::support_fingerprint(a_nonzero),
            Self::support_fingerprint(b_nonzero),
        );
        if let Some((units, tickets, unit_words, total)) = cache.get(key) {
            return KernelGating {
                units: units.clone(),
                tickets: tickets.clone(),
                unit_words: unit_words.clone(),
                total_pairs: *total,
            };
        }
        // Fresh gating: presence bitsets drive the per-entry ticket
        // flags, degree-support masks the run-level gating. A unit's
        // entries are grouped into contiguous p-runs (entries are sorted
        // `(p, q, i, j)` within the unit and `q = target - p` is forced),
        // and a run is kept only when its own `(p, target - p)` degree
        // pair is supported — unit-level gating would drag every other p's
        // entries through tests that always fail. Surviving entries carry
        // pre-packed orientation flags, and the run walk simultaneously
        // collects each active entry's row targets into the owning unit's
        // word set, so no sweep re-derives anything.
        let words = a_series.basis.len().div_ceil(64);
        let mut presence = vec![0u64; 2 * words];
        let (a_present, b_present) = presence.split_at_mut(words);
        let mut a_deg = [0u64; 2];
        let mut b_deg = [0u64; 2];
        for &i in a_nonzero {
            a_present[i / 64] |= 1u64 << (i % 64);
            let d = table.degree_of(i);
            debug_assert!(d < 128, "degree masks cover degrees 0..127");
            a_deg[d / 64] |= 1u64 << (d % 64);
        }
        for &j in b_nonzero {
            b_present[j / 64] |= 1u64 << (j % 64);
            let d = table.degree_of(j);
            debug_assert!(d < 128, "degree masks cover degrees 0..127");
            b_deg[d / 64] |= 1u64 << (d % 64);
        }
        let value = Self::build_gating(
            table,
            table.entries(),
            table.entries(),
            table.decomp_indices_rel(),
            a_present,
            b_present,
            a_deg,
            b_deg,
        );
        let KernelGating {
            ref units,
            ref tickets,
            ref unit_words,
            total_pairs,
        } = value;
        cache.insert(
            key,
            (
                Arc::clone(units),
                Arc::clone(tickets),
                Arc::clone(unit_words),
                total_pairs,
            ),
        );
        value
    }

    /// 128-bit fingerprint of a support list, the gating cache's key half.
    /// Content-addressed: an entry's baked per-entry presence results are
    /// only valid for the exact supports they were resolved from, so two
    /// calls may share a cache entry only when their support lists are
    /// identical. Length is mixed in and both streams are avalanche-
    /// finished (murmur3 fmix64), so lists of different lengths or content
    /// collide with negligible probability.
    /// The gating cache's key fingerprint for one support list (murmur3-
    /// finish FNV pairs). Length is mixed in and both streams are
    /// avalanche-finished, so distinct lists collide with negligible
    /// probability.
    #[inline]
    fn support_fingerprint(list: &[usize]) -> [u64; 2] {
        let mut h1 = 0xcbf2_9ce4_8422_2325u64 ^ (list.len() as u64).rotate_left(32);
        let mut h2 = 0x9ae1_6a3b_2f90_404fu64 ^ (list.len() as u64);
        for &x in list {
            let v = x as u64;
            h1 = (h1 ^ v).wrapping_mul(0x1000_0000_01b3);
            h2 = (h2 ^ v).wrapping_mul(0x94d0_49bb_1331_11eb);
        }
        h1 ^= h1 >> 33;
        h1 = h1.wrapping_mul(0xff51_afd7_ed55_8ccd);
        h1 ^= h1 >> 33;
        h2 ^= h2 >> 33;
        h2 = h2.wrapping_mul(0xff51_afd7_ed55_8ccd);
        h2 ^= h2 >> 33;
        [h1, h2]
    }

    /// Walks the decomposition table's units once and resolves the gating
    /// for the given presence bitsets and degree-support masks: the flat
    /// per-entry ticket list (orientation bits packed) plus, in the same
    /// walk, the per-word transposition of every active contribution into
    /// its target word's write class. Shared by both prologue variants.
    ///
    /// `entries` (public i/j) drives the p-run walk's degree lookups —
    /// degrees are layout-independent — while `presence_entries` must
    /// carry the i/j relabeling that matches the presence bitsets' index
    /// space (public entries for the public prologue, the class-order's
    /// relabeled table for the class one). The two tables share order,
    /// index space, and `decomp_start`s, so the ticket's entry index is
    /// layout-independent either way. `decomp_tbl`/`pos_degree` are the
    /// working layout's scatter indices and position→degree map, and
    /// `space` is the result buffer's position count.
    fn build_gating(
        table: &FeasibleDecompositions<U>,
        entries: &[Entry],
        presence_entries: &[Entry],
        decomp_tbl: &[u32],
        a_present: &[u64],
        b_present: &[u64],
        a_deg: [u64; 2],
        b_deg: [u64; 2],
    ) -> KernelGating {
        debug_assert!(
            table.len() < (1 << 30),
            "entry indices must fit a ticket's 30 bits"
        );
        let mut tickets: Vec<u32> = Vec::new();
        let mut unit_words: Vec<u32> = Vec::new();
        let mut units: Vec<UnitRange> = Vec::with_capacity(table.units().len());
        // Transient per-unit rel bitset over the unit's degree slice,
        // cleared by iterating its set bits after extraction.
        let mut rel_bits: Vec<u64> = Vec::new();
        for unit in table.units().iter() {
            let t = unit.target as usize;
            let rs = table.degree_start(t) as u32;
            let ticket_start = tickets.len() as u32;
            let mut pairs = 0u32;
            rel_bits.clear();
            rel_bits.resize((table.degree_start(t + 1) as u32 - rs).div_ceil(64) as usize, 0);
            let mut cur_p = u8::MAX;
            let mut run_start = unit.start;
            // Real entries only: `unit.end` is the trailing sentinel's
            // slot (its decomp_start closes the last run's last
            // decomposition range via the +1 span).
            for ei in unit.start..unit.end {
                let p = table.degree_of(entries[ei as usize].i as usize) as u8;
                if p == cur_p {
                    continue;
                }
                if cur_p != u8::MAX {
                    Self::push_run(
                        presence_entries,
                        decomp_tbl,
                        a_present,
                        b_present,
                        &mut tickets,
                        &mut rel_bits,
                        &mut pairs,
                        a_deg,
                        b_deg,
                        cur_p,
                        t,
                        run_start,
                        ei,
                    );
                }
                cur_p = p;
                run_start = ei;
            }
            if cur_p != u8::MAX {
                Self::push_run(
                    presence_entries,
                    decomp_tbl,
                    a_present,
                    b_present,
                    &mut tickets,
                    &mut rel_bits,
                    &mut pairs,
                    a_deg,
                    b_deg,
                    cur_p,
                    t,
                    run_start,
                    unit.end,
                );
            }
            // Extract the unit's ascending word positions (pos = rs + rel)
            // and clear the bitset for the next unit.
            let mut emitted = 0u32;
            for (w, bits) in rel_bits.iter_mut().enumerate() {
                let mut b = *bits;
                while b != 0 {
                    let bit = b.trailing_zeros();
                    b &= b - 1;
                    unit_words.push(rs + (w as u32) * 64 + bit);
                    emitted += 1;
                }
                *bits = 0;
            }
            // A unit with no active contributions (every run gated out, or
            // every active entry's row empty) is omitted: no tickets, no
            // words, no work.
            if emitted > 0 {
                units.push(UnitRange {
                    rs,
                    td: unit.target,
                    ticket_start,
                    ticket_end: tickets.len() as u32,
                    pairs,
                });
            } else {
                debug_assert_eq!(pairs, 0);
                tickets.truncate(ticket_start as usize);
            }
        }
        // Units of one degree are ordered by CONTENT BYTES (the table's
        // unit order), not by position — the per-unit concatenation is not
        // ascending. Sort the union so the flat list is globally ascending:
        // the sink/scatter-set walk order must stay byte-identical to the
        // per-word fan-in's (which emitted active positions in ascending
        // order).
        unit_words.sort_unstable();
        let total_pairs = units.iter().map(|u| u.pairs as usize).sum();
        KernelGating {
            units: Arc::from(units),
            tickets: Arc::from(tickets),
            unit_words: Arc::from(unit_words),
            total_pairs,
        }
    }

    /// Class-space variant of [`Self::kernel_prologue_cached`]: the
    /// support lists are class-indexed, so the presence bitsets are
    /// class-positioned, the degree masks read through the ordering's
    /// relabeled degree table, and the per-word transposition targets
    /// class positions. The memo key — the exact class-indexed support
    /// lists — is the class-space image of the public variant's, so each
    /// working mode resolves its own gating entries (the class order is a
    /// pure function of the basis, so the two modes' transpositions are
    /// the same permutation apart).
    pub(super) fn kernel_prologue_cached_class(
        a_series: &LieSeries<T, U>,
        a_nonzero_cls: &[usize],
        b_nonzero_cls: &[usize],
        order: &ClassOrder,
        cache: &mut GatingCache,
    ) -> KernelGating {
        Self::kernel_prologue_cached_class_core(
            &a_series.feasible_decompositions,
            a_series.basis.len(),
            a_nonzero_cls,
            b_nonzero_cls,
            order,
            cache,
        )
    }

    /// T-free core of [`Self::kernel_prologue_cached_class`] (see
    /// `sweep_words_serial_core` for why the letter type must not reach
    /// the parallel kernels).
    pub(super) fn kernel_prologue_cached_class_core(
        table: &FeasibleDecompositions<U>,
        basis_len: usize,
        a_nonzero_cls: &[usize],
        b_nonzero_cls: &[usize],
        order: &ClassOrder,
        cache: &mut GatingCache,
    ) -> KernelGating {

        // Full-support fast path (class space): same predicate and cutoff
        // logic as the public prologue — the supports must be EXACTLY the
        // prefix `0..cutoff` (length alone misfires on batch-recorded node
        // lists containing degree-`max_degree` positions; see the public
        // path's comment and the pooled-dag leaf-reuse regression test),
        // per-unit gating cached on the ordering.
        let cutoff = table.degree_start(table.max_degree());
        let full_prefix = |s: &[usize]| s.len() == cutoff && s.iter().all(|&p| p < cutoff);
        if full_prefix(a_nonzero_cls) && full_prefix(b_nonzero_cls) {
            let (units, tickets, unit_words) = order.full_support_gating_class(table);
            return KernelGating {
                total_pairs: units.iter().map(|u| u.pairs as usize).sum(),
                units,
                tickets,
                unit_words,
            };
        }

        let key = (
            Self::support_fingerprint(a_nonzero_cls),
            Self::support_fingerprint(b_nonzero_cls),
        );
        if let Some((units, tickets, unit_words, total)) = cache.get(key) {
            return KernelGating {
                units: units.clone(),
                tickets: tickets.clone(),
                unit_words: unit_words.clone(),
                total_pairs: *total,
            };
        }
        // Class-positioned presence bitsets (indexed by class positions)
        // and relabeled degrees; the rest of the gating walk is shared
        // with the public prologue (see `build_gating`).
        let words = basis_len.div_ceil(64);
        let mut presence = vec![0u64; 2 * words];
        let (a_present, b_present) = presence.split_at_mut(words);
        let mut a_deg = [0u64; 2];
        let mut b_deg = [0u64; 2];
        for &i in a_nonzero_cls {
            a_present[i / 64] |= 1u64 << (i % 64);
            let d = order.degree_cls()[i] as usize;
            debug_assert!(d < 128, "degree masks cover degrees 0..127");
            a_deg[d / 64] |= 1u64 << (d % 64);
        }
        for &j in b_nonzero_cls {
            b_present[j / 64] |= 1u64 << (j % 64);
            let d = order.degree_cls()[j] as usize;
            debug_assert!(d < 128, "degree masks cover degrees 0..127");
            b_deg[d / 64] |= 1u64 << (d % 64);
        }
        let value = Self::build_gating(
            table,
            table.entries(),
            order.entries_cls(),
            order.decomp_cls(),
            a_present,
            b_present,
            a_deg,
            b_deg,
        );
        let KernelGating {
            ref units,
            ref tickets,
            ref unit_words,
            total_pairs,
        } = value;
        cache.insert(
            key,
            (
                Arc::clone(units),
                Arc::clone(tickets),
                Arc::clone(unit_words),
                total_pairs,
            ),
        );
        value
    }
    /// Resolves one maximal p-run of a unit: run-level gating on the
    /// degree-support masks, then the per-entry presence tests whose
    /// results become the run's packed tickets — and, in the same walk,
    /// each active entry's decomposition row marks its target rels in the
    /// owning unit's transient bitset and counts the unit's pairs.
    /// Tickets are built in table order, so each word's contribution
    /// sequence is a subsequence of the entry stream in table order:
    /// per-word float summation order is provably unchanged.
    #[allow(clippy::too_many_arguments)]
    fn push_run(
        presence_entries: &[Entry],
        decomp_tbl: &[u32],
        a_present: &[u64],
        b_present: &[u64],
        tickets: &mut Vec<u32>,
        rel_bits: &mut [u64],
        pairs: &mut u32,
        a_deg: [u64; 2],
        b_deg: [u64; 2],
        p: u8,
        t: usize,
        run_start: u32,
        run_end: u32,
    ) {
        let pu = p as usize;
        let qu = t - pu;
        let o1 = a_deg[pu / 64] >> (pu % 64) & 1 != 0 && b_deg[qu / 64] >> (qu % 64) & 1 != 0;
        let o2 = b_deg[pu / 64] >> (pu % 64) & 1 != 0 && a_deg[qu / 64] >> (qu % 64) & 1 != 0;
        if !o1 && !o2 {
            return;
        }
        // The same bitmap expressions the sweeps used to re-evaluate per
        // kernel call, resolved once here: `o1`/`o2` ANDed in identically,
        // presence tested per index.
        for ei in run_start..run_end {
            // i/j in the presence bitsets' index space (public for the
            // public prologue, class positions for the class one).
            let entry = &presence_entries[ei as usize];
            let (i, j) = (entry.i as usize, entry.j as usize);
            let p_active = o1
                && a_present[i / 64] & (1u64 << (i % 64)) != 0
                && b_present[j / 64] & (1u64 << (j % 64)) != 0;
            let q_active = o2
                && a_present[j / 64] & (1u64 << (j % 64)) != 0
                && b_present[i / 64] & (1u64 << (i % 64)) != 0;
            if !p_active && !q_active {
                continue;
            }
            tickets.push(
                ei | if p_active { TICKET_P_ACTIVE } else { 0 }
                    | if q_active { TICKET_Q_ACTIVE } else { 0 },
            );
            // Mark the entry's row targets in the owning unit's bitset
            // (rels are degree-slice relative — the same in both working
            // spaces) and count the pairs.
            let from = entry.decomp_start as usize;
            let to = presence_entries[ei as usize + 1].decomp_start as usize;
            for &rel in &decomp_tbl[from..to] {
                rel_bits[(rel / 64) as usize] |= 1u64 << (rel % 64);
                *pairs += 1;
            }
        }
    }
}
