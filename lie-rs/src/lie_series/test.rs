use ordered_float::NotNan;
use rstest::rstest;

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

/// INVARIANT CHECK on the real tables: within one degree-target slice,
/// two different units must never scatter onto the same basis word (the
/// batch's packs cut at unit boundaries and rely on per-word
/// single-writer accumulation).
#[test]
fn debug_real_table_unit_word_ownership() {
    for (d, m) in [
        (2usize, 4usize),
        (3usize, 4usize),
        (3usize, 5usize),
        (3usize, 8usize),
    ] {
        let basis = lyndon_rs::lyndon::LyndonBasis::<u8>::new(
            d,
            lyndon_rs::lyndon::Sort::Lexicographical,
        )
        .generate_basis(m);
        let len = basis.len();
        let a = LieSeries::<u8, NotNan<f64>>::new(
            basis.clone(),
            (0..len)
                .map(|k| NotNan::new((k % 17) as f64 / 19.0 - 0.4).unwrap())
                .collect(),
        );
        let _b = LieSeries::<u8, NotNan<f64>>::new(
            basis,
            (0..len)
                .map(|k| NotNan::new((k % 13) as f64 / 17.0 - 0.3).unwrap())
                .collect(),
        );
        let fd = &a.feasible_decompositions;
        let mut owner: std::collections::HashMap<(u8, u32), usize> =
            std::collections::HashMap::new();
        let mut violations = 0usize;
        for (uid, unit) in fd.units().iter().enumerate() {
            let span = fd.entry_span(unit);
            for (entry, next) in span[..span.len() - 1].iter().zip(span[1..].iter()) {
                let from = entry.decomp_start as usize;
                let to = next.decomp_start as usize;
                for &rel in &fd.decomp_indices_rel()[from..to] {
                    let key = (unit.target, rel);
                    match owner.get(&key) {
                        Some(prev) if *prev != uid => {
                            violations += 1;
                            if violations <= 6 {
                                println!(
                                    "VIOLATION d={d} m={m}: degree {} rel {} by units {} and {}",
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
            "d={d} m={m}: units={} rels={} violations={}",
            fd.units().len(),
            fd.decomp_indices_rel().len(),
            violations
        );
        assert_eq!(violations, 0);
    }
}

/// Anagram-class scatter locality of the commutation kernel's write set:
/// per anagram unit, how many distinct 64-byte lines its stores touch and
/// how far they spread inside the public degree slice, versus the same
/// stores routed into the class-contiguous (internal) layout.
#[test]
#[ignore = "layout stats: run explicitly with --ignored --nocapture"]
fn dump_scatter_locality_stats() {
    for (d, m) in [(2usize, 12usize), (3, 8), (4, 10), (12, 2)] {
        let basis = lyndon_rs::LyndonBasis::<u8>::new(d, lyndon_rs::Sort::Lexicographical);
        let words = basis.generate_basis(m);
        let len = words.len();
        let a = LieSeries::new(words.clone(), vec![NotNan::new(1.0).unwrap(); len]);
        let table = &a.feasible_decompositions;
        let decomp = table.decomp_indices_rel();
        let entries = table.entries();
        let mut cur_lines = 0usize;
        let mut cls_lines = 0usize;
        let mut total_rmw = 0usize;
        let mut spread_max = 0usize;
        let mut report = String::new();
        for unit in table.units() {
            let (rs, re) = table.degree_range(unit.target);
            let slice = (re - rs) as usize;
            let mut words_seen = std::collections::BTreeSet::new();
            let mut rmw = 0usize;
            for ei in unit.start..unit.end {
                let from = entries[ei as usize].decomp_start as usize;
                let to = entries[ei as usize + 1].decomp_start as usize;
                for &rel in &decomp[from..to] {
                    words_seen.insert(rel as usize);
                }
                rmw += to - from;
            }
            let distinct = words_seen.len();
            if distinct == 0 {
                continue;
            }
            let lo = *words_seen.first().unwrap();
            let hi = *words_seen.last().unwrap();
            let spread = hi - lo + 1;
            let lines_touched = words_seen
                .iter()
                .map(|&w| (rs as usize + w) / 8)
                .collect::<std::collections::HashSet<_>>()
                .len();
            let lines_cls = (distinct + 7) / 8;
            cur_lines += lines_touched;
            cls_lines += lines_cls;
            total_rmw += rmw;
            spread_max = spread_max.max(spread);
            report.push_str(&format!(
                "    unit t={:2} words={:5} rmw={:6} slice={:7} spread={:6} ({:5.1}% of slice) lines_cur={:5} lines_cls={:5}\n",
                unit.target,
                distinct,
                rmw,
                slice,
                spread,
                100.0 * spread as f64 / slice.max(1) as f64,
                lines_touched,
                lines_cls,
            ));
        }
        println!(
            "{d}x{m}: len={len} store_rmw={total_rmw} lines_cur_sum={cur_lines} lines_cls_sum={cls_lines} ({:.2}x reduction) max_unit_spread={spread_max}",
            cur_lines as f64 / cls_lines.max(1) as f64
        );
        print!("{report}");
    }
}

use super::*;

#[rstest]
#[case(2, 2, vec![1, 2, 3], vec![4, 5, 6], vec![0, 0, -3])]
#[case(2, 2, vec![3, 2, 1], vec![1, 2, 3], vec![0, 0, 4])]
#[case(2, 3, vec![1, 2, 3, 4, 5], vec![6, 7, 8, 9, 10], vec![0, 0, -5, -10, 5])]
#[case(3, 3,
    vec![1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    vec![5, 3, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    vec![0, 0, 0, -7, -14, -7, 0, 0, 0, 0, 0, 0, 0, 0])]
#[case(3, 3,
    vec![1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    vec![0, 0, 0, -7, -14, -7, 0, 0, 0, 0, 0, 0, 0, 0],
    vec![0, 0, 0, 0, 0, 0, -7, -14, 14, 14, 49, 42, -14, 21])]
#[case(3, 4,
    vec![
        1, 2, 3, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0
    ],
    vec![
        5, 3, 1, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0
    ],
    vec![
        0, 0, 0, -7, -14, -7, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0
    ],
)]
#[case(3, 4,
    vec![
        1, 2, 3, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0
    ],
    vec![
        0, 0, 0, -7, -14, -7, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0
    ],
    vec![
        0, 0, 0, 0, 0, 0, -7, -14,
        14, 14, 49, 42, -14, 21, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
    ],
)]
fn test_lie_series_commutation(
    #[case] num_generators: usize,
    #[case] basis_depth: usize,
    #[case] a_coefficients: Vec<i128>,
    #[case] b_coefficients: Vec<i128>,
    #[case] expected_coefficients: Vec<i128>,
) {
    use lyndon_rs::lyndon::{LyndonBasis, Sort};
    use num_rational::Ratio;

    let basis = LyndonBasis::<u8>::new(num_generators, Sort::Lexicographical)
        .generate_basis(basis_depth);
    let a_coefficients = a_coefficients
        .into_iter()
        .map(Ratio::<i128>::from_integer)
        .collect::<Vec<_>>();
    let a = LieSeries::new(basis.clone(), a_coefficients);

    let b_coefficients = b_coefficients
        .into_iter()
        .map(Ratio::<i128>::from_integer)
        .collect::<Vec<_>>();
    let b = LieSeries::new(basis.clone(), b_coefficients);
    let expected_coefficients = expected_coefficients
        .into_iter()
        .map(Ratio::<i128>::from_integer)
        .collect::<Vec<_>>();

    let series = comm![a, b];
    assert_eq!(series.coefficients.len(), expected_coefficients.len());
    dbg!(&series.coefficients);
    for (i, c) in series.coefficients.iter().enumerate() {
        assert_eq!(
            *c, expected_coefficients[i],
            "{i}: {c:?} != {:?}",
            expected_coefficients[i]
        );
    }
}

#[test]
fn decompositions_contain_no_zero_coefficients() {
    use lyndon_rs::lyndon::{LyndonBasis, Sort};
    use ordered_float::NotNan;

    for (d, m) in [(2usize, 6usize), (3, 5), (4, 4)] {
        let words = LyndonBasis::<u8>::new(d, Sort::Lexicographical).generate_basis(m);
        let series = LieSeries::<u8, NotNan<f64>>::new(words, Vec::new());
        assert!(
            series
                .feasible_decompositions
                .iter_entries()
                .flat_map(|(_, _, _, coefficients)| coefficients)
                .all(|c| !c.is_zero()),
            "zero coefficient found in decompositions for d={d}, m={m}"
        );
    }
}

/// Independently recomputes the decomposition of every feasible canonical
/// pair and pins the compact table's contents (indices, coefficients,
/// ordering, feasibility, and canonicality) against it, plus the mirrored
/// (negated) query path.
#[test]
fn feasible_decompositions_match_independent_recomputation() {
    use lyndon_rs::lyndon::{LyndonBasis, Sort};
    use ordered_float::NotNan;

    for (d, m) in [(2usize, 6usize), (3, 5), (4, 4)] {
        let words = LyndonBasis::<u8>::new(d, Sort::Lexicographical).generate_basis(m);
        let series = LieSeries::<u8, NotNan<f64>>::new(words, Vec::new());
        let degree = |x: usize| series.commutator_basis[x].degree();
        let index_of = |term: &CommutatorTerm<NotNan<f64>, u8>| {
            series.commutator_basis_index_map[&term.unit_hash()]
        };
        let basis_set: HashSet<_> = series
            .commutator_basis
            .iter()
            .map(CommutatorTerm::unit_hash)
            .collect();

        let mut feasible = 0;
        for unit in series.feasible_decompositions.units() {
            let (p, _q, target) = (unit.p_mask, unit.target as usize, unit.target as usize);
            assert_eq!(
                p.iter().map(|w| w.count_ones()).sum::<u32>() as usize >= 1,
                true,
                "empty unit"
            );
            assert_eq!(
                unit.gamma.iter().map(|&x| x as usize).sum::<usize>(),
                target,
                "unit target degree mismatch"
            );
            let entries: Vec<_> = series
                .feasible_decompositions
                .unit_iter(unit)
                .map(|(e, _, _)| (e.i as usize, e.j as usize))
                .collect();
            assert!(
                entries.windows(2).all(|w| w[0] < w[1]),
                "unit {:?} entries not sorted by (i, j)",
                unit.gamma
            );

            for (i, j) in entries {
                assert!(i < j, "non-canonical pair stored");
                assert!(
                    p[(degree(i) / 64) as usize] >> (degree(i) % 64) & 1 == 1,
                    "pair ({i}, {j}) in a unit of the wrong degrees"
                );
                assert_eq!(
                    degree(i) + degree(j),
                    target,
                    "pair ({i}, {j}) in a unit of the wrong target degree"
                );
                assert!(degree(i) + degree(j) <= m, "infeasible pair stored");

                let mut term = comm![&series.commutator_basis[i], &series.commutator_basis[j]];
                term.lyndon_sort();
                let expected: Vec<_> = term
                    .lyndon_basis_decomposition(&basis_set)
                    .into_iter()
                    .filter(|x| !x.coefficient().is_zero())
                    .collect();

                let expected_indices: Vec<u32> =
                    expected.iter().map(|x| index_of(x) as u32).collect();
                let expected_coefficients: Vec<_> =
                    expected.iter().map(|x| x.coefficient().clone()).collect();

                // The canonical query returns the stored data unflagged;
                // comparing through `decomposition` exercises the
                // degree-block lookup as well as the contents.
                let (canonical_indices, canonical_coefficients, swapped) =
                    series.decomposition(i, j).expect("canonical pair query");
                assert!(!swapped, "(i={i}, j={j}) canonical query flagged");
                assert_eq!(
                    canonical_indices,
                    &expected_indices[..],
                    "(i={i}, j={j}) indices"
                );
                assert_eq!(
                    canonical_coefficients,
                    &expected_coefficients[..],
                    "(i={i}, j={j}) coeffs"
                );
                feasible += 1;

                // The mirrored query returns the same (canonical) data,
                // flagged so the caller negates it into
                // [basis[j], basis[i]] orientation.
                let (mirrored_indices, mirrored_coefficients, swapped) =
                    series.decomposition(j, i).expect("mirrored pair query");
                assert!(swapped, "(i={i}, j={j}) mirrored query not flagged");
                assert_eq!(mirrored_indices, &expected_indices[..]);
                assert_eq!(
                    mirrored_coefficients,
                    &expected_coefficients[..],
                    "(i={i}, j={j}) mirrored data differs from canonical"
                );
            }
        }
        assert_eq!(series.feasible_decompositions.len(), feasible);
    }
}

/// The fused kernel evaluates both bracket orientations of every canonical
/// pair, so `[A, A]` must vanish exactly — every pair contributes
/// `c * (a_min * a_max - a_max * a_min) = 0`.
#[test]
fn commutator_of_series_with_itself_is_zero() {
    use lyndon_rs::lyndon::{LyndonBasis, Sort};
    use ordered_float::NotNan;

    for (d, m) in [(2usize, 6usize), (3, 5)] {
        let words: Vec<LyndonWord<u8>> =
            LyndonBasis::<u8>::new(d, Sort::Lexicographical).generate_basis(m);
        let coefficients: Vec<_> = (0..words.len())
            .map(|i: usize| NotNan::new(((i % 7 + 1) as f64) / 3.0 - 1.1).unwrap())
            .collect();
        let series: LieSeries<u8, NotNan<f64>> = LieSeries::new(words, coefficients);
        let result: LieSeries<u8, NotNan<f64>> = series.commutator(&series);
        assert!(
            result.coefficients.iter().all(|c| c.is_zero()),
            "non-zero coefficient in [A, A] for d={d}, m={m}"
        );
    }
}

/// ADVERSARIAL (per-unit division): the gating's per-unit structure
/// must be LOSSLESS and ORDERED against an independent walk of the
/// table with presence resolved straight from the support lists:
/// (a) the unit ticket ranges partition the flat ticket list without
/// overlap, (b) each unit's orientation flags match the independent
/// presence resolution per entry, (c) the unit word sets PARTITION the
/// active word positions — pairwise disjoint, ascending per unit, and
/// their concatenation ascending (so CollectSink order is unchanged) —
/// and equal the set of positions the independent walk produces, (d)
/// each unit's contribution sequence — its tickets expanded in order,
/// each entry's row (rel, dp) pairs — is exactly the independent
/// walk's subsequence for that unit's words, in table-entry order (the
/// bit-exactness contract), and (e) total_pairs matches. A
/// misbucketed, lost, duplicated, or reordered contribution fails
/// here with a distinct message per invariant.
#[test]
fn write_class_gating_transposition_is_lossless_and_ordered() {
    use lyndon_rs::lyndon::{LyndonBasis, Sort};
    use ordered_float::NotNan;

    use super::ClassOrderedCommutation;

    for (d, m) in [(2usize, 4usize), (3usize, 5usize), (3usize, 8usize)] {
        let words: Vec<LyndonWord<u8>> =
            LyndonBasis::<u8>::new(d, Sort::Lexicographical).generate_basis(m);
        let len = words.len();
        let series: LieSeries<u8, NotNan<f64>> = LieSeries::new(
            words.clone(),
            (0..len)
                .map(|i| NotNan::new(((i % 7 + 1) as f64) / 5.0 - 0.7).unwrap())
                .collect(),
        );
        let table = &series.feasible_decompositions;
        let entries = table.entries();
        let decomp = table.decomp_indices_rel();
        let degrees = table.index_degrees();
        let order = series.class_order();

        // Independent support pairs: full, one-sided, random halves,
        // single letters, empty.
        let cutoff = table.degree_start(table.max_degree());
        let mut support_sets: Vec<Vec<usize>> = vec![Vec::new()];
        support_sets.push((0..cutoff).collect());
        support_sets.push((0..cutoff).filter(|i| i % 2 == 0).collect());
        support_sets.push((0..cutoff).filter(|i| i % 3 != 1).collect());
        support_sets.push((0..cutoff).filter(|i| degrees[*i] == 1).collect());
        support_sets.push(vec![0]);

        for a_support in &support_sets {
            for b_support in &support_sets {
                // Independent active-contribution walk over the table,
                // in table-entry order: presence straight from the
                // support lists (degree-mask gating is implied by
                // membership).
                let present = |s: &[usize]| {
                    let mut p = vec![false; len];
                    for &i in s {
                        p[i] = true;
                    }
                    p
                };
                let (ap, bp) = (present(a_support), present(b_support));
                // (pos, entry, dp) in table-entry order.
                let mut expected: Vec<(usize, usize, usize)> = Vec::new();
                for (ei, entry) in entries[..entries.len() - 1].iter().enumerate() {
                    let (i, j) = (entry.i as usize, entry.j as usize);
                    let p_active = ap[i] && bp[j];
                    let q_active = ap[j] && bp[i];
                    if !p_active && !q_active {
                        continue;
                    }
                    let from = entry.decomp_start as usize;
                    let to = entries[ei + 1].decomp_start as usize;
                    let rs = table.degree_start(degrees[i] as usize + degrees[j] as usize);
                    for (k, &rel) in decomp[from..to].iter().enumerate() {
                        expected.push((rs + rel as usize, ei, from + k));
                    }
                }
                // Active entries with orientation, table order (the
                // tickets must be exactly this list).
                let mut expected_tickets: Vec<(usize, bool, bool)> = Vec::new();
                for (ei, entry) in entries[..entries.len() - 1].iter().enumerate() {
                    let (i, j) = (entry.i as usize, entry.j as usize);
                    let p_active = ap[i] && bp[j];
                    let q_active = ap[j] && bp[i];
                    if p_active || q_active {
                        expected_tickets.push((ei, p_active, q_active));
                    }
                }

                // Public-mode gating: positions are public basis
                // indices.
                let mut cache = GatingCache::default();
                let gating = LieSeries::<u8, NotNan<f64>>::kernel_prologue_cached(
                    &series,
                    a_support,
                    b_support,
                    &mut cache,
                );
                let ctx = "public";

                // (a) Ticket-range partition: consecutive, ordered,
                // non-overlapping, covering the ticket list.
                let mut cursor = 0usize;
                for (ui, u) in gating.units.iter().enumerate() {
                    assert_eq!(
                        u.ticket_start as usize, cursor,
                        "{ctx}: unit {ui} ticket range not contiguous"
                    );
                    assert!(
                        u.ticket_end > u.ticket_start,
                        "{ctx}: unit {ui} emitted with no tickets"
                    );
                    assert_eq!(
                        degrees[u.rs as usize], u.td,
                        "{ctx}: unit {ui} rs/td degree mismatch"
                    );
                    cursor = u.ticket_end as usize;
                }
                assert_eq!(
                    cursor,
                    gating.tickets.len(),
                    "{ctx}: ticket ranges do not cover the ticket list"
                );

                // (b) Orientation flags match the independent walk,
                // ticket list = expected active entries in table order.
                assert_eq!(
                    gating.tickets.len(),
                    expected_tickets.len(),
                    "{ctx}: active-entry count mismatch"
                );
                for (tp, &(ei, want_p, want_q)) in
                    expected_tickets.iter().enumerate()
                {
                    let ticket = gating.tickets[tp];
                    assert_eq!(
                        (ticket & TICKET_INDEX_MASK) as usize, ei,
                        "{ctx}: ticket {tp} entry mismatch"
                    );
                    assert_eq!(
                        ticket & TICKET_P_ACTIVE != 0, want_p,
                        "{ctx}: p_active flag mismatch at entry {ei}"
                    );
                    assert_eq!(
                        ticket & TICKET_Q_ACTIVE != 0, want_q,
                        "{ctx}: q_active flag mismatch at entry {ei}"
                    );
                }

                // (c)+(d) Word-set partition AND per-unit sequence:
                // reconstruct each unit's word set and contribution
                // sequence from its rows; sets must be pairwise
                // disjoint, within the unit's degree slice, and their
                // union must equal both the gating's flat `unit_words`
                // list and the independent walk's positions; each
                // unit's (entry, dp) sequence must be exactly the
                // independent walk's subsequence for the unit's words,
                // in table-entry order (the bit-exactness contract).
                let mut all_positions: Vec<usize> = Vec::new();
                let mut seen_pairs: Vec<(usize, usize)> = Vec::new();
                for (ui, u) in gating.units.iter().enumerate() {
                    let mut set: Vec<usize> = Vec::new();
                    let mut unit_pairs: Vec<(usize, usize)> = Vec::new();
                    for tp in u.ticket_start as usize..u.ticket_end as usize {
                        let ei = (gating.tickets[tp] & TICKET_INDEX_MASK) as usize;
                        let entry = entries[ei];
                        let from = entry.decomp_start as usize;
                        let to = entries[ei + 1].decomp_start as usize;
                        for (k, &rel) in decomp[from..to].iter().enumerate() {
                            let pos = u.rs as usize + rel as usize;
                            set.push(pos);
                            unit_pairs.push((ei, from + k));
                        }
                    }
                    set.sort_unstable();
                    set.dedup();
                    assert!(
                        set.iter().all(|&p| {
                            degrees[p] == u.td
                                && (u.rs as usize..table.degree_start(u.td as usize + 1))
                                    .contains(&p)
                        }),
                        "{ctx}: unit {ui} word set outside its degree slice"
                    );
                    for &p in &set {
                        assert!(
                            !all_positions.contains(&p),
                            "{ctx}: position {p} written by two units"
                        );
                    }
                    all_positions.extend_from_slice(&set);
                    seen_pairs.extend_from_slice(&unit_pairs);
                    // The unit's sequence must contain only its own
                    // words, in table-entry order.
                    let mut last_ei = None::<usize>;
                    for &(ei, dp) in &unit_pairs {
                        let rel = decomp[dp];
                        assert!(
                            set.binary_search(&(u.rs as usize + rel as usize)).is_ok(),
                            "{ctx}: unit {ui} contributes to a word outside its set"
                        );
                        assert!(
                            last_ei.is_none_or(|prev| prev <= ei),
                            "{ctx}: unit {ui} contributions out of table order"
                        );
                        last_ei = Some(ei);
                    }
                }
                let mut sorted_all = all_positions.clone();
                sorted_all.sort_unstable();
                let got_flat: Vec<usize> =
                    gating.unit_words.iter().map(|&p| p as usize).collect();
                assert_eq!(
                    sorted_all, got_flat,
                    "{ctx}: flat unit_words list differs from the union of the rows' positions"
                );
                let mut want_positions: Vec<usize> =
                    expected.iter().map(|&(p, _, _)| p).collect();
                want_positions.sort_unstable();
                want_positions.dedup();
                assert_eq!(
                    got_flat, want_positions,
                    "{ctx}: unit word-set partition differs from the independent walk's positions"
                );
                // (d) global sequence: every unit's contributions
                // concatenated = the independent walk's (entry, dp)
                // pairs in table-entry order.
                let mut want_pairs: Vec<(usize, usize)> = expected
                    .iter()
                    .map(|&(_, ei, dp)| (ei, dp))
                    .collect();
                assert_eq!(
                    seen_pairs, want_pairs,
                    "{ctx}: unit contribution sequences differ from the independent walk"
                );
                let _ = &mut want_pairs;

                // (e) Pair count.
                assert_eq!(
                    gating.total_pairs,
                    expected.len(),
                    "{ctx}: contribution count mismatch (supports {a_support:?}/{b_support:?})"
                );

                // Class-mode gating: positions are class positions; the
                // same invariants must hold against the class-space
                // image of the expected list.
                let inv = order.inv();
                let a_cls: Vec<usize> = a_support.iter().map(|&i| inv[i] as usize).collect();
                let b_cls: Vec<usize> = b_support.iter().map(|&i| inv[i] as usize).collect();
                let gating_cls = LieSeries::<u8, NotNan<f64>>::kernel_prologue_cached_class(
                    &series,
                    &a_cls,
                    &b_cls,
                    &order,
                    &mut cache,
                );
                let degree_cls = order.degree_cls();
                let decomp_cls = order.decomp_cls();
                let entries_cls = order.entries_cls();
                // Class-space image of the expected walk: map public
                // positions through inv, keep table-entry order.
                let mut expected_cls: Vec<(usize, usize, usize)> = Vec::new();
                let perm = order.perm();
                for (ei, entry) in entries_cls[..entries_cls.len() - 1].iter().enumerate() {
                    let (i, j) = (entry.i as usize, entry.j as usize);
                    // entries_cls endpoints are CLASS positions; the
                    // presence vectors are indexed by PUBLIC position
                    // (perm: class -> public).
                    let p_active = ap[perm[i] as usize] && bp[perm[j] as usize];
                    let q_active = ap[perm[j] as usize] && bp[perm[i] as usize];
                    if !p_active && !q_active {
                        continue;
                    }
                    let from = entry.decomp_start as usize;
                    let to = entries_cls[ei + 1].decomp_start as usize;
                    let rs = table.degree_start(
                        degree_cls[i] as usize + degree_cls[j] as usize,
                    );
                    for &rel in decomp_cls[from..to].iter() {
                        expected_cls.push((rs + rel as usize, ei, from as usize));
                    }
                }
                // (a)+(c) for class space: ticket ranges cover the
                // list; the flat word list is globally ascending and
                // equals the independent class-space walk's positions.
                let mut cursor_cls = 0usize;
                for (ui, u) in gating_cls.units.iter().enumerate() {
                    assert_eq!(
                        u.ticket_start as usize, cursor_cls,
                        "class: unit {ui} ticket range not contiguous"
                    );
                    cursor_cls = u.ticket_end as usize;
                }
                assert_eq!(cursor_cls, gating_cls.tickets.len(), "class: ticket ranges do not cover the list");
                let mut want_cls: Vec<usize> =
                    expected_cls.iter().map(|&(p, _, _)| p).collect();
                want_cls.sort_unstable();
                want_cls.dedup();
                let got_cls: Vec<usize> = gating_cls
                    .unit_words
                    .iter()
                    .map(|&p| p as usize)
                    .collect();
                let mut last_cls = None::<usize>;
                for &p in &got_cls {
                    assert!(
                        last_cls.is_none_or(|prev| prev < p),
                        "class: unit word list not globally ascending at {p}"
                    );
                    last_cls = Some(p);
                }
                assert_eq!(
                    got_cls, want_cls,
                    "class: unit word-set partition differs from the independent walk"
                );
                // (e) Pair count.
                assert_eq!(
                    gating_cls.total_pairs,
                    expected_cls.len(),
                    "class: class-mode contribution count mismatch"
                );
            }
        }
    }
}

/// ADVERSARIAL (write-access division): the kernel sweep must equal an
/// INDEPENDENT per-word fan-in reference — built by walking the table's
/// entries in table order with presence resolved straight from the
/// support lists — bit for bit, across coefficient types (raw-float
/// fast path and exact rationals), layouts (public direct and
/// class-contiguous), thread counts, and adversarial inputs (all-zero,
/// single-hot, planted exact cancellations). Distinct failure: this
/// guards arithmetic and per-word accumulation order, not gating
/// structure (the transposition test above).
#[test]
fn write_class_sweep_matches_independent_fanin_reference() {
    use lyndon_rs::lyndon::{LyndonBasis, Sort};
    use ordered_float::NotNan;
    use std::collections::BTreeMap;

    use super::ClassOrderedCommutation;

    // The reference: result[w] += c * term over the table's entries in
    // table order, term = (p_active ? a_i*b_j : 0) - (q_active ?
    // a_j*b_i : 0), entries with neither orientation skipped. Same add
    // sequence per word as the kernel's write classes, derived without
    // touching the gating.
    fn reference<U>(
        series: &LieSeries<u8, U>,
        a: &[U],
        a_present: &[bool],
        b: &[U],
        b_present: &[bool],
    ) -> Vec<U>
    where
        U: Clone + Default + Mul<Output = U> + AddAssign + Neg<Output = U> + PartialEq,
    {
        let table = &series.feasible_decompositions;
        let entries = table.entries();
        let coeffs = table.decomp_coeffs();
        let mut acc: BTreeMap<usize, U> = BTreeMap::new();
        for (ei, entry) in entries[..entries.len() - 1].iter().enumerate() {
            let (i, j) = (entry.i as usize, entry.j as usize);
            let p_active = a_present[i] && b_present[j];
            let q_active = a_present[j] && b_present[i];
            if !p_active && !q_active {
                continue;
            }
            let term = if p_active {
                let mut t = a[i].clone() * b[j].clone();
                if q_active {
                    t += -(a[j].clone() * b[i].clone());
                }
                t
            } else {
                -(a[j].clone() * b[i].clone())
            };
            let from = entry.decomp_start as usize;
            let to = entries[ei + 1].decomp_start as usize;
            let rs = table.degree_start(
                table.degree_of(i) + table.degree_of(j),
            );
            for dp in from..to {
                // Absolute decomp position -> the row's target: the
                // relative array is stored per entry in row order.
                let rel = table.decomp_indices_rel()[dp];
                let w = rs + rel as usize;
                let contrib = coeffs[dp].clone() * term.clone();
                *acc.entry(w).or_insert_with(U::default) += contrib;
            }
        }
        let mut out = vec![U::default(); series.coefficients.len()];
        for (w, v) in acc {
            out[w] = v;
        }
        out
    }

    for (d, m) in [(2usize, 4usize), (3usize, 6usize)] {
        let words: Vec<LyndonWord<u8>> =
            LyndonBasis::<u8>::new(d, Sort::Lexicographical).generate_basis(m);
        let len = words.len();
        let table_deg = |i: usize| words[i].len();
        let _ = table_deg;

        // Adversarial value sets: zeros (support holes), signed values
        // that plant exact cancellations between entries feeding the
        // same word, and asymmetric magnitudes.
        let mk = |f: &dyn Fn(usize) -> f64| {
            (0..len)
                .map(|i| NotNan::new(f(i)).unwrap())
                .collect::<Vec<_>>()
        };
        let a_sets = [
            mk(&|_| 0.0),
            mk(&|i| if i == 0 { 1.0 } else { 0.0 }),
            mk(&|i| ((i % 5 + 1) as f64) / 3.0 - 1.0),
            mk(&|i| if i % 2 == 0 { 1.0 } else { -1.0 }),
        ];
        let b_sets = [
            mk(&|_| 0.0),
            mk(&|i| if i == 1 { 1.0 } else { 0.0 }),
            mk(&|i| ((i * 3 + 2) % 7) as f64 / 2.0 - 1.5),
            mk(&|i| if i % 3 == 0 { 2.0 } else { -0.5 }),
        ];

        let cutoff_len = len; // full lists; the kernel filters by value
        for a_coefficients in &a_sets {
            for b_coefficients in &b_sets {
                let a: LieSeries<u8, NotNan<f64>> =
                    LieSeries::new(words.clone(), a_coefficients.clone());
                let b: LieSeries<u8, NotNan<f64>> =
                    LieSeries::new(words.clone(), b_coefficients.clone());
                let a_nonzero = a.nonzero_coefficient_indices(a_coefficients);
                let b_nonzero = b.nonzero_coefficient_indices(b_coefficients);
                let ap = {
                    let mut p = vec![false; cutoff_len];
                    for &i in &a_nonzero {
                        p[i] = true;
                    }
                    p
                };
                let bp = {
                    let mut p = vec![false; cutoff_len];
                    for &i in &b_nonzero {
                        p[i] = true;
                    }
                    p
                };
                let expected =
                    reference(&a, a_coefficients, &ap, b_coefficients, &bp);

                for threads in [1usize, 4usize] {
                    let pool = rayon::ThreadPoolBuilder::new()
                        .num_threads(threads)
                        .build()
                        .expect("pool");
                    pool.install(|| {
                        // Public direct kernel.
                        let mut out =
                            vec![NotNan::<f64>::default(); len];
                        LieSeries::commutator_coefficients_with_nonzero(
                            &a,
                            a_coefficients,
                            &a_nonzero,
                            b_coefficients,
                            &b_nonzero,
                            &mut out,
                        );
                        assert_eq!(
                            out, expected,
                            "public kernel differs (d={d}, m={m}, t={threads}, \
                             a={:?}, b={:?})",
                            a_coefficients.iter().map(|x| x.into_inner()).collect::<Vec<_>>(),
                            b_coefficients.iter().map(|x| x.into_inner()).collect::<Vec<_>>()
                        );

                        // Class-contiguous kernel: operands class-
                        // ordered, result mapped back through `inv`.
                        let order = a.class_order();
                        let a_cls = a.class_coefficients(&order, a_coefficients);
                        let b_cls = b.class_coefficients(&order, b_coefficients);
                        let a_nz_cls: Vec<usize> =
                            a_nonzero.iter().map(|&i| order.inv()[i] as usize).collect();
                        let b_nz_cls: Vec<usize> =
                            b_nonzero.iter().map(|&i| order.inv()[i] as usize).collect();
                        let mut out_cls =
                            vec![NotNan::<f64>::default(); len];
                        a.class_commutation(
                            &order,
                            &a_cls,
                            &a_nz_cls,
                            &b_cls,
                            &b_nz_cls,
                            &mut out_cls,
                            &mut GatingCache::default(),
                        );
                        let public = a.public_coefficients(&order, &out_cls);
                        assert_eq!(
                            public, expected,
                            "class kernel differs (d={d}, m={m}, t={threads})"
                        );
                    });
                }
            }
        }
    }
}

/// ADVERSARIAL (write-access division): a word class stays OWNED under
/// total cancellation — `a == b` makes every term exactly zero, yet
/// every reachable word is still swept (one single-writer store of the
/// exact zero) and still reported by the collecting variant: the
/// reported targets are a WRITE-set superset, never a value-filtered
/// set. Distinct failure: guards target ownership and the superset
/// property, not arithmetic (the fan-in test above).
#[test]
fn write_class_canceled_word_stays_owned_and_reported() {
    use lyndon_rs::lyndon::{LyndonBasis, Sort};
    use ordered_float::NotNan;

    for (d, m) in [(2usize, 5usize), (3usize, 4usize)] {
        let words: Vec<LyndonWord<u8>> =
            LyndonBasis::<u8>::new(d, Sort::Lexicographical).generate_basis(m);
        let len = words.len();
        let coeffs: Vec<_> = (0..len)
            .map(|i| NotNan::new(((i % 9 + 1) as f64) / 4.0 - 1.1).unwrap())
            .collect();
        let a: LieSeries<u8, NotNan<f64>> =
            LieSeries::new(words.clone(), coeffs.clone());
        let b: LieSeries<u8, NotNan<f64>> = LieSeries::new(words, coeffs);
        let full: Vec<usize> = (0..len).collect();

        // [A, A] = 0 exactly (float products are commutative), but
        // every word with a feasible decomposition is still written.
        let mut out = vec![NotNan::<f64>::default(); len];
        LieSeries::commutator_coefficients_with_nonzero(
            &a,
            &a.coefficients,
            &full,
            &b.coefficients,
            &full,
            &mut out,
        );
        for (w, v) in out.iter().enumerate() {
            assert_eq!(*v, NotNan::<f64>::default(), "word {w} not exactly zero");
        }

        // The collecting variant must still report every written word:
        // value-filtered targets would be unsound for list reuse (a
        // canceled target can go nonzero again when values change).
        let mut dirty = vec![0u64; len / 64 + 1];
        let mut targets = Vec::new();
        LieSeries::commutator_coefficients_with_nonzero_collecting(
            &a,
            &a.coefficients,
            &full,
            &b.coefficients,
            &full,
            &mut out,
            &mut dirty,
            &mut targets,
        );
        // Independent expectation: every word targeted by any active
        // contribution of any entry (full supports -> every entry).
        let table = &a.feasible_decompositions;
        let entries = table.entries();
        let decomp = table.decomp_indices_rel();
        let degrees = table.index_degrees();
        let mut expected: Vec<usize> = Vec::new();
        for (ei, entry) in entries[..entries.len() - 1].iter().enumerate() {
            let (i, j) = (entry.i as usize, entry.j as usize);
            let from = entry.decomp_start as usize;
            let to = entries[ei + 1].decomp_start as usize;
            let rs = table.degree_start(degrees[i] as usize + degrees[j] as usize);
            for &rel in &decomp[from..to] {
                expected.push(rs + rel as usize);
            }
        }
        expected.sort_unstable();
        expected.dedup();
        // The sink reports kernel-VISIBLE positions only: the degree-
        // `max_degree` tail lies above the support cutoff (nothing can
        // consume those values in a truncated BCH fold), so both sides
        // filter at `degree_start(max_degree)`.
        let cutoff = table.degree_start(table.max_degree());
        expected.retain(|&w| w < cutoff);
        targets.sort_unstable();
        assert_eq!(
            targets, expected,
            "collecting targets lost canceled words (write-set superset violated)"
        );
    }
}

/// Bundling decides only which thread runs which anagram unit; the
/// per-word accumulation order is unchanged, so the result must be
/// bit-identical for any thread count. Guards the work-balanced bundle
/// builder: the parallel sweep must never reorder, duplicate, or drop
/// an accumulation.
#[test]
fn commutator_is_bit_identical_across_thread_counts() {
    use lyndon_rs::lyndon::{LyndonBasis, Sort};
    use ordered_float::NotNan;

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

/// The trait's internal working mode (`class_commutation`) must agree
/// bit-for-bit with the direct kernel after one
/// `public_coefficients` epilogue, at every thread count; the
/// collecting variant must report the class-indexed image of the
/// direct layout's first-touch sequence.
#[test]
fn class_commutation_round_trip_is_bit_identical() {
    use lyndon_rs::lyndon::{LyndonBasis, Sort};
    use ordered_float::NotNan;

    use super::ClassOrderedCommutation;

    for (d, m, force) in [(2usize, 4usize, true), (3usize, 10usize, false)] {
        let words: Vec<LyndonWord<u8>> =
            LyndonBasis::<u8>::new(d, Sort::Lexicographical).generate_basis(m);
        let a_coefficients: Vec<_> = (0..words.len())
            .map(|i: usize| NotNan::new(((i % 11 + 1) as f64) / 7.0 - 0.9).unwrap())
            .collect();
        let b_coefficients: Vec<_> = (0..words.len())
            .map(|i: usize| NotNan::new(((i * 5 + 3) % 13) as f64 / 6.0 - 1.3).unwrap())
            .collect();
        let mut direct = LieSeries::new(words.clone(), a_coefficients.clone());
        let class = LieSeries::new(words, b_coefficients.clone());
        if force {
            Arc::make_mut(&mut direct.feasible_decompositions).clear_class_order();
        }
        let order = class.class_order();
        assert_eq!(order.inv().len(), class.coefficients.len());

        let a_cls = class.class_coefficients(&order, &a_coefficients);
        let b_cls = class.class_coefficients(&order, &b_coefficients);
        let inv_of = |i: usize| order.inv()[i] as usize;
        let a_nonzero = direct.nonzero_coefficient_indices(&a_coefficients);
        let b_nonzero = direct.nonzero_coefficient_indices(&b_coefficients);
        let a_nz_cls: Vec<usize> = a_nonzero.iter().copied().map(inv_of).collect();
        let b_nz_cls: Vec<usize> = b_nonzero.iter().copied().map(inv_of).collect();

        let run = |threads: usize, series: &LieSeries<u8, NotNan<f64>>| {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .expect("thread pool");
            let mut out = vec![NotNan::<f64>::default(); series.coefficients.len()];
            pool.install(|| {
                LieSeries::commutator_coefficients(
                    series,
                    &a_coefficients,
                    &b_coefficients,
                    &mut out,
                )
            });
            out
        };

        let reference = run(1, &direct);
        for threads in [1usize, 2, 4, 8] {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .expect("thread pool");
            let mut result = vec![NotNan::<f64>::default(); class.coefficients.len()];
            pool.install(|| {
                class.class_commutation(
                    &order,
                    &a_cls,
                    &a_nz_cls,
                    &b_cls,
                    &b_nz_cls,
                    &mut result,
                    &mut GatingCache::default(),
                )
            });
            let public = class.public_coefficients(&order, &result);
            assert_eq!(
                reference, public,
                "class round trip differs for d={d}, m={m}, threads={threads}"
            );
        }

        // Collecting variant: class-indexed first-touch sequence.
        let mut direct_result = vec![NotNan::<f64>::default(); direct.coefficients.len()];
        let mut direct_dirty = vec![0u64; direct.coefficients.len() / 64 + 1];
        let mut direct_targets = Vec::new();
        LieSeries::commutator_coefficients_with_nonzero_collecting(
            &direct,
            &a_coefficients,
            &a_nonzero,
            &b_coefficients,
            &b_nonzero,
            &mut direct_result,
            &mut direct_dirty,
            &mut direct_targets,
        );
        let mut class_result = vec![NotNan::<f64>::default(); class.coefficients.len()];
        let mut class_dirty = vec![0u64; class.coefficients.len() / 64 + 1];
        let mut class_targets = Vec::new();
        class.class_commutation_with_nonzero_collecting(
            &order,
            &a_cls,
            &a_nz_cls,
            &b_cls,
            &b_nz_cls,
            &mut class_result,
            &mut class_dirty,
            &mut class_targets,
        );
        for (k, &src) in order.inv().iter().enumerate() {
            assert_eq!(
                direct_result[k], class_result[src as usize],
                "collecting result differs for d={d}, m={m} at public {k}"
            );
        }
        // Class targets are internal positions: map them back with
        // `perm` (internal -> public) before comparing. Under the
        // write-access division each variant emits its targets sorted
        // ascending in ITS OWN layout's position order (one store per
        // word class, classes in position order), so the public
        // sequences differ by the layout permutation — the invariant is
        // the same target SET, each position reported exactly once
        // (per-word accumulation order is guarded by the bit-identical
        // result comparison above).
        let perm = order.perm();
        let mut relabeled: Vec<usize> = class_targets
            .iter()
            .copied()
            .map(|p| perm[p] as usize)
            .collect();
        relabeled.sort_unstable();
        let mut direct_sorted = direct_targets.clone();
        direct_sorted.sort_unstable();
        assert_eq!(direct_sorted, relabeled, "collecting targets differ");
    }
}

/// The class-contiguous scatter layout must produce bit-identical
/// results to the direct layout, at every thread count, on both kernel
/// entry points: the layout is a pure relabeling of write addresses and
/// never reorders the per-word accumulation sequence. `(2, 4)` forces
/// the layout on a small table (fast build); `(3, 10)` qualifies
/// through the real slice-size threshold (its degree-10 slice exceeds
/// 4096 words).
#[test]
fn commutator_is_bit_identical_across_scatter_layouts() {
    use lyndon_rs::lyndon::{LyndonBasis, Sort};
    use ordered_float::NotNan;

    for (d, m, force) in [(2usize, 4usize, true), (3usize, 10usize, false)] {
        let words: Vec<LyndonWord<u8>> =
            LyndonBasis::<u8>::new(d, Sort::Lexicographical).generate_basis(m);
        let a_coefficients: Vec<_> = (0..words.len())
            .map(|i: usize| NotNan::new(((i % 11 + 1) as f64) / 7.0 - 0.9).unwrap())
            .collect();
        let b_coefficients: Vec<_> = (0..words.len())
            .map(|i: usize| NotNan::new(((i * 5 + 3) % 13) as f64 / 6.0 - 1.3).unwrap())
            .collect();
        let mut direct = LieSeries::new(words.clone(), a_coefficients.clone());
        let mut class = LieSeries::new(words, b_coefficients.clone());
        if force {
            Arc::make_mut(&mut class.feasible_decompositions).force_class_order();
        } else {
            // (3, 10) auto-qualifies through the threshold: strip the
            // reference's layout so the pair compares direct vs class
            // on the same table.
            Arc::make_mut(&mut direct.feasible_decompositions).clear_class_order();
        }
        assert!(
            direct
                .feasible_decompositions
                .cached_class_order()
                .is_none(),
            "expected the reference series to run the direct layout"
        );
        assert!(
            class.feasible_decompositions.cached_class_order().is_some(),
            "expected the layout series to carry a class order"
        );
        let len = direct.coefficients.len();

        let run = |series: &LieSeries<u8, NotNan<f64>>, threads: usize| {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .expect("thread pool");
            let mut out = vec![NotNan::<f64>::default(); len];
            pool.install(|| {
                LieSeries::commutator_coefficients(
                    series,
                    &a_coefficients,
                    &b_coefficients,
                    &mut out,
                )
            });
            out
        };

        let reference = run(&direct, 1);
        for threads in [1usize, 2, 4, 8] {
            let out = run(&class, threads);
            assert_eq!(
                reference, out,
                "class layout differs from direct for d={d}, m={m}, threads={threads}"
            );
        }

        // Collecting entry point: result, nonzero set, and first-touch
        // order must all match the direct layout exactly.
        let collect = |series: &LieSeries<u8, NotNan<f64>>| {
            let a_nonzero = series.nonzero_coefficient_indices(&a_coefficients);
            let b_nonzero = series.nonzero_coefficient_indices(&b_coefficients);
            let mut result = vec![NotNan::<f64>::default(); len];
            let mut dirty = vec![0u64; len / 64 + 1];
            let mut targets = Vec::new();
            LieSeries::commutator_coefficients_with_nonzero_collecting(
                series,
                &a_coefficients,
                &a_nonzero,
                &b_coefficients,
                &b_nonzero,
                &mut result,
                &mut dirty,
                &mut targets,
            );
            (result, targets)
        };
        let (direct_result, direct_targets) = collect(&direct);
        let (class_result, class_targets) = collect(&class);
        assert_eq!(direct_result, class_result, "collecting result differs");
        assert_eq!(direct_targets, class_targets, "collecting targets differ");
    }
}
