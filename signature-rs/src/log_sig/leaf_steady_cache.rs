use super::{Arc, ClassOrder, HashMap, Mutex, OnceLock, SharedDirtyLists, SharedNodeLists};

/// Key: (table identity, displacement support). The support is the
/// exact sorted-deduped vector — never a hash — so distinct supports
/// can never collide by construction.
pub(super) type LeafPlanKey = (u64, Vec<usize>);

pub(super) struct LeafSteadyLists {
    /// The converged reachable-support plan (public basis indices).
    pub(super) a_plan: Vec<usize>,
    /// The plan dag's node scatter lists (class-contiguous indices),
    /// exactly as the last collecting rebuild left them. `Arc`-shared:
    /// a hit adopts these lists by reference into the group's dag —
    /// zero copy — and the dag's first list mutation copies
    /// (copy-on-write), so the shared snapshot is never mutated.
    pub(super) nonzeros: SharedNodeLists,
    /// The matching per-node dirty bitsets. `class_collect_kernel`
    /// zeroes the dirty words at entry, so the post-collection
    /// state is a valid pre-state for any later collection. Shared
    /// exactly like [`Self::nonzeros`].
    pub(super) dirty: SharedDirtyLists,
    /// The class order the lists are expressed in (the table's shared
    /// ordering — the same `Arc` every dag over this table installs).
    pub(super) order: Arc<ClassOrder>,
}

static LEAF_STEADY_PLANS: OnceLock<Mutex<HashMap<LeafPlanKey, Arc<LeafSteadyLists>>>> =
    OnceLock::new();

/// Bound on distinct cached supports. Uniform workloads — the common
/// case, every group sharing one displacement support — hold ONE
/// entry; adversarially diverse workloads trade cache memory for
/// planning sweeps, and clearing at the cap keeps the map bounded.
const LEAF_STEADY_PLANS_CAP: usize = 4096;

pub(super) fn get(table_id: u64, b: &[usize]) -> Option<Arc<LeafSteadyLists>> {
    let map = LEAF_STEADY_PLANS.get_or_init(|| Mutex::new(HashMap::new()));
    let guard = map.lock().expect("leaf steady plan cache mutex");
    guard.get(&(table_id, b.to_vec())).map(Arc::clone)
}

pub(super) fn store(
    table_id: u64,
    b: &[usize],
    a_plan: Vec<usize>,
    nonzeros: SharedNodeLists,
    dirty: SharedDirtyLists,
    order: Arc<ClassOrder>,
) {
    let map = LEAF_STEADY_PLANS.get_or_init(|| Mutex::new(HashMap::new()));
    let mut guard = map.lock().expect("leaf steady plan cache mutex");
    if guard.len() >= LEAF_STEADY_PLANS_CAP {
        guard.clear();
    }
    guard.insert(
        (table_id, b.to_vec()),
        Arc::new(LeafSteadyLists {
            a_plan,
            nonzeros,
            dirty,
            order,
        }),
    );
}
#[cfg(test)]
mod tests {
    //! Leaf steady-plan cache tests: fixed-point content equality, support
    //! keying, mutation isolation, Arc sharing and copy-on-write behavior.
    //! Test-only cache resets live here too (`clear_for_tests`).

    use super::*;
    use super::super::engine::{
        SeriesTemplate, leaf_steady_adopt_into, leaf_steady_plan, leaf_steady_plan_cached,
    };
    use crate::LogSignature;
    use crate::LogSignatureBuilder;
    use crate::commutator_dag::CommutatorDag;
    use crate::commutator_dag::set_cohort_off;
    use num_rational::Ratio;
    use ordered_float::NotNan;
    use std::sync::Arc;

    /// Clears the process-wide leaf plan cache between tests: every test
    /// controls whether its first derivation is a cache miss or a hit.
    fn clear_for_tests() {
        if let Some(map) = LEAF_STEADY_PLANS.get() {
            map.lock().expect("leaf steady plan cache mutex").clear();
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


    /// A (template, fresh dag) pair over (d, m) with Ratio coefficients —
    /// exact arithmetic makes any planning divergence a value divergence.

    fn leaf_plan_mk(
        d: usize,
        m: usize,
    ) -> (
        SeriesTemplate<Ratio<i64>, u8>,
        CommutatorDag<Ratio<i64>>,
        usize,
    ) {
        let builder = LogSignatureBuilder::<u8>::new()
            .with_num_dimensions(d)
            .with_max_degree(m);
        let LogSignature {
            series,
            bch_series: _,
            dag,
        } = builder.build::<Ratio<i64>>();
        let basis_len = series.basis.len();
        (SeriesTemplate(series), dag, basis_len)
    }


    /// Displacements with a controlled support: nonzero exactly at
    /// `coords` (LCG values elsewhere zero). Ratio keeps everything exact.

    fn leaf_plan_rhs(
        basis_len: usize,
        coords: &[usize],
        seed: u64,
    ) -> Vec<Ratio<i64>> {
        let mut s = seed;
        let mut lcg = || {
            s = s.wrapping_mul(6364136223846793005).wrapping_add(1);
            ((s >> 33) % 9) as i64 - 4
        };
        (0..basis_len)
            .map(|k| {
                if coords.contains(&k) {
                    let v = lcg();
                    Ratio::from_integer(if v == 0 { 3 } else { v })
                } else {
                    Ratio::from_integer(0)
                }
            })
            .collect()
    }


    /// Fixed-point content equality: a cache-hit adoption must produce the
    /// SAME plan and the SAME node lists a fresh derivation produces (the
    /// lists are a pure function of the supports — this pins that, and it
    /// is what makes every other test's bit-identity non-vacuous).
    #[test]

    fn leaf_plan_cache_adopted_state_equals_fresh_derivation() {
        clear_for_tests();
        let (template, mut dag_cold, basis_len) = leaf_plan_mk(3, 5);
        let first = leaf_plan_rhs(basis_len, &[0, 1, 2], 0xAB1);
        let b: Vec<usize> = {
            let mut v = template.0.nonzero_coefficient_indices(&first);
            v.sort_unstable();
            v
        };
        let zero = vec![Ratio::<i64>::default(); basis_len];

        // MISS: derive + store.
        let a_cold = leaf_steady_plan_cached(&template, &mut dag_cold, &zero, &first, &b);
        let (nz_cold, dirty_cold) = dag_cold.steady_lists_snapshot();
        let snap = get(template.0.feasible_decompositions_handle().table_id(), &b)
            .expect("derivation must populate the cache");
        assert_eq!(snap.a_plan, a_cold, "cached plan must be the derived plan");
        assert_eq!(snap.nonzeros, nz_cold, "cached lists must be the derived lists");
        assert_eq!(snap.dirty, dirty_cold, "cached dirty state must match");

        // HIT: adopt into a fresh dag; the adopted state must equal the
        // derived state content-for-content.
        let (_, mut dag_warm, _) = leaf_plan_mk(3, 5);
        let a_warm = leaf_steady_plan_cached(&template, &mut dag_warm, &zero, &first, &b);
        assert_eq!(a_warm, a_cold, "adopted plan diverged from derived plan");
        let (nz_warm, _) = dag_warm.steady_lists_snapshot();
        assert_eq!(nz_warm, nz_cold, "adopted node lists diverged from derived lists");

        // A second dag derived WITHOUT the cache (cleared) must reach the
        // same fixed point — the cache cannot shift the fixed point.
        clear_for_tests();
        let (_, mut dag_fresh, _) = leaf_plan_mk(3, 5);
        let a_fresh = leaf_steady_plan(&template, &mut dag_fresh, &zero, &first, &b);
        assert_eq!(a_fresh, a_cold, "uncached derivation diverged");
    }


    /// Support-signature correctness: two different displacement supports
    /// must get DIFFERENT cache entries and DIFFERENT plans — no false hit
    /// by key collision or staleness — and each must equal its own fresh
    /// derivation.
    #[test]

    fn leaf_plan_cache_distinguishes_distinct_supports() {
        clear_for_tests();
        let (template, mut dag_a, basis_len) = leaf_plan_mk(3, 5);
        let zero = vec![Ratio::<i64>::default(); basis_len];

        let dense = leaf_plan_rhs(basis_len, &[0, 1, 2], 0xC22);
        let sparse = leaf_plan_rhs(basis_len, &[0], 0xC33);
        let support = |r: &Vec<Ratio<i64>>| {
            let mut v = template.0.nonzero_coefficient_indices(r);
            v.sort_unstable();
            v
        };
        let b_dense = support(&dense);
        let b_sparse = support(&sparse);
        assert!(b_sparse.len() < b_dense.len(), "test needs distinct supports");

        let a_dense = leaf_steady_plan_cached(&template, &mut dag_a, &zero, &dense, &b_dense);
        let a_sparse = leaf_steady_plan_cached(&template, &mut dag_a, &zero, &sparse, &b_sparse);
        assert_ne!(
            a_dense, a_sparse,
            "different supports must not share a plan"
        );

        // Cross-check both against fresh uncached derivations.
        clear_for_tests();
        let (_, mut dag_f, _) = leaf_plan_mk(3, 5);
        assert_eq!(
            leaf_steady_plan(&template, &mut dag_f, &zero, &dense, &b_dense),
            a_dense,
            "dense plan diverged from fresh derivation"
        );
        assert_eq!(
            leaf_steady_plan(&template, &mut dag_f, &zero, &sparse, &b_sparse),
            a_sparse,
            "sparse plan diverged from fresh derivation"
        );
    }


    /// Zero-subset edge: a displacement whose support is a STRICT SUBSET
    /// of its lane-0 neighbor's first slice must fold correctly under the
    /// group's shared (superset) plan — cached or derived — and the
    /// results must match sequential folding bit for bit.
    #[test]

    fn leaf_plan_cache_subset_displacement_matches_sequential() {
        clear_for_tests();
        let (d, m) = (3usize, 5usize);
        let (template, mut dag, basis_len) = leaf_plan_mk(d, m);
        let zero = vec![Ratio::<i64>::default(); basis_len];

        // Lane 0: dense displacement; lane 1: strict subset support.
        let dense = leaf_plan_rhs(basis_len, &[0, 1, 2], 0xD41);
        let subset = leaf_plan_rhs(basis_len, &[1], 0xD42);
        let b0 = {
            let mut v = template.0.nonzero_coefficient_indices(&dense);
            v.sort_unstable();
            v
        };
        // The documented uniformity condition: every slice's support is a
        // subset of lane 0's first slice's support.
        assert!(
            template
                .0
                .nonzero_coefficient_indices(&subset)
                .iter()
                .all(|k| b0.contains(k)),
            "test construction: subset support must be contained"
        );

        // Shared plan via the cache (first call: derive + store).
        let a_plan = leaf_steady_plan_cached(&template, &mut dag, &zero, &dense, &b0);
        let mk_series = || {
            let mut s = template.0.clone();
            s.coefficients = vec![Ratio::<i64>::default(); basis_len];
            s
        };

        // Lane 0 (dense) on the plan dag; lane 1 (subset) on a dag that
        // ADOPTS the cached lists (the fold-loop path).
        let mut series0 = mk_series();
        let dense_slice: Vec<&[Ratio<i64>]> = vec![&dense];
        dag.fold_batch(&mut series0, &dense_slice, &a_plan, &b0);

        let (_, mut dag2, _) = leaf_plan_mk(d, m);
        assert!(
            leaf_steady_adopt_into(&mut dag2, template.0.feasible_decompositions_handle().table_id(), &a_plan, &b0),
            "adoption must succeed for a same-config dag"
        );
        let mut series1 = mk_series();
        let subset_slice: Vec<&[Ratio<i64>]> = vec![&subset];
        dag2.fold_batch(&mut series1, &subset_slice, &a_plan, &b0);

        // Ground truth: sequential folds.
        let mut seq0 = LogSignatureBuilder::<u8>::new()
            .with_num_dimensions(d)
            .with_max_degree(m)
            .build::<Ratio<i64>>();
        seq0.concatenate_assign_coefficients(&dense);
        let mut seq1 = LogSignatureBuilder::<u8>::new()
            .with_num_dimensions(d)
            .with_max_degree(m)
            .build::<Ratio<i64>>();
        seq1.concatenate_assign_coefficients(&subset);

        let bits = |v: &[Ratio<i64>]| -> Vec<(i64, i64)> {
            v.iter().map(|r| (*r.numer(), *r.denom())).collect()
        };
        assert_eq!(bits(&series0.coefficients), bits(&seq0.series.coefficients), "dense lane diverged from sequential");
        assert_eq!(bits(&series1.coefficients), bits(&seq1.series.coefficients), "subset lane diverged from sequential");
    }


    /// Mutation isolation: the snapshot handed out by the cache is never
    /// mutated by folds running on adopting dags (fold_batch re-records
    /// the DAG's OWN lists with its exact scatter sets — the snapshot must
    /// stay the untouched fixed point), and a later group adopting after
    /// earlier groups folded still matches a fresh derivation.
    #[test]

    fn leaf_plan_cache_snapshot_immune_to_fold_mutation() {
        clear_for_tests();
        let (template, mut dag, basis_len) = leaf_plan_mk(3, 5);
        let zero = vec![Ratio::<i64>::default(); basis_len];
        let r1 = leaf_plan_rhs(basis_len, &[0, 1], 0xE51);
        let b1 = {
            let mut v = template.0.nonzero_coefficient_indices(&r1);
            v.sort_unstable();
            v
        };
        let a_plan = leaf_steady_plan_cached(&template, &mut dag, &zero, &r1, &b1);
        let table_id = template.0.feasible_decompositions_handle().table_id();
        let before = get(table_id, &b1).expect("snapshot present");
        let before_nz = before.nonzeros.clone();

        // Fold on the adopting dag (fold_batch re-records ITS lists).
        let mut series = template.0.clone();
        series.coefficients = vec![Ratio::<i64>::default(); basis_len];
        let slices: Vec<&[Ratio<i64>]> = vec![&r1];
        dag.fold_batch(&mut series, &slices, &a_plan, &b1);

        let after = get(table_id, &b1).expect("snapshot survives");
        assert_eq!(
            after.nonzeros, before_nz,
            "fold mutated the shared cached snapshot"
        );
        assert!(Arc::ptr_eq(&before, &after), "cache entry must not be replaced by a fold");

        // A later adoption (after the fold) still reproduces the fresh
        // derivation's fold bit for bit.
        clear_for_tests();
        let (_, mut dag_fresh, _) = leaf_plan_mk(3, 5);
        let a_fresh = leaf_steady_plan(&template, &mut dag_fresh, &zero, &r1, &b1);
        let mut series_fresh = template.0.clone();
        series_fresh.coefficients = vec![Ratio::<i64>::default(); basis_len];
        dag_fresh.fold_batch(&mut series_fresh, &slices, &a_fresh, &b1);
        let bits = |v: &[Ratio<i64>]| -> Vec<(i64, i64)> {
            v.iter().map(|r| (*r.numer(), *r.denom())).collect()
        };
        assert_eq!(bits(&series.coefficients), bits(&series_fresh.coefficients), "post-fold adoption diverged from fresh derivation");
    }


    /// Arc-level mutation isolation: a cache hit installs the snapshot's
    /// lists by REFERENCE (`Arc::ptr_eq` — zero copy, the (d) contract at
    /// 3x5 scale), a fold on the adopting dag replaces the dag's Arc
    /// wholesale (never mutating through the shared one), and the
    /// snapshot survives byte-identical for the next adopting dag — which
    /// folds bit-identically to the first.
    #[test]

    fn leaf_plan_lists_adopted_by_arc_and_isolated_from_folds() {
        clear_for_tests();
        let (template, mut dag_a, basis_len) = leaf_plan_mk(3, 5);
        let zero = vec![Ratio::<i64>::default(); basis_len];
        let r1 = leaf_plan_rhs(basis_len, &[0, 1], 0xF61);
        let b1 = {
            let mut v = template.0.nonzero_coefficient_indices(&r1);
            v.sort_unstable();
            v
        };
        let table_id = template.0.feasible_decompositions_handle().table_id();

        // MISS: derive + store, then drop the deriving dag so the cache is
        // the snapshot's only holder.
        let a_plan = leaf_steady_plan_cached(&template, &mut dag_a, &zero, &r1, &b1);
        let snap = get(table_id, &b1).expect("snapshot present");
        let snapshot_nz: Vec<Vec<usize>> = (*snap.nonzeros).clone();
        drop(dag_a);

        // HIT: adopt into a fresh dag — by reference, not by copy.
        let (_, mut dag_b, _) = leaf_plan_mk(3, 5);
        assert!(leaf_steady_adopt_into(&mut dag_b, table_id, &a_plan, &b1));
        assert!(
            Arc::ptr_eq(&dag_b.nonzeros, &snap.nonzeros),
            "adoption must share the snapshot Arc, not copy it"
        );
        assert_eq!(
            Arc::strong_count(&snap.nonzeros),
            2,
            "exactly the cache and the adopting dag hold the lists"
        );

        // Fold on the adopting dag: fold_batch re-records the dag's OWN
        // lists — through a NEW Arc, never through the shared snapshot.
        let mut series = template.0.clone();
        series.coefficients = vec![Ratio::<i64>::default(); basis_len];
        let slices: Vec<&[Ratio<i64>]> = vec![&r1];
        dag_b.fold_batch(&mut series, &slices, &a_plan, &b1);
        assert!(
            !Arc::ptr_eq(&dag_b.nonzeros, &snap.nonzeros),
            "the fold's re-record must replace the dag's Arc, not mutate the shared one"
        );
        let snap_after = get(table_id, &b1).expect("snapshot survives the fold");
        assert!(Arc::ptr_eq(&snap, &snap_after), "cache entry must not be replaced");
        assert_eq!(*snap_after.nonzeros, snapshot_nz, "fold mutated the shared snapshot");

        // A post-fold adoption still shares the pristine snapshot and
        // folds bit-identically.
        let (_, mut dag_c, _) = leaf_plan_mk(3, 5);
        assert!(leaf_steady_adopt_into(&mut dag_c, table_id, &a_plan, &b1));
        assert!(Arc::ptr_eq(&dag_c.nonzeros, &snap.nonzeros));
        let mut series_c = template.0.clone();
        series_c.coefficients = vec![Ratio::<i64>::default(); basis_len];
        dag_c.fold_batch(&mut series_c, &slices, &a_plan, &b1);
        let bits = |v: &[Ratio<i64>]| -> Vec<(i64, i64)> {
            v.iter().map(|r| (*r.numer(), *r.denom())).collect()
        };
        assert_eq!(
            bits(&series.coefficients),
            bits(&series_c.coefficients),
            "post-fold adoption diverged from the first adopting fold"
        );
    }


    /// Concurrent adoption: 8 threads adopt the SAME snapshot into
    /// distinct dags and fold — all bit-identical to each other and to a
    /// sequential fold, with the shared snapshot's contents untouched
    /// (lockless `Arc` reads; no mutation through the shared Arc).
    #[test]

    fn leaf_plan_lists_concurrent_adoption_bit_identical() {
        clear_for_tests();
        let (template, mut dag0, basis_len) = leaf_plan_mk(3, 5);
        let zero = vec![Ratio::<i64>::default(); basis_len];
        let r1 = leaf_plan_rhs(basis_len, &[0, 2], 0xF71);
        let b1 = {
            let mut v = template.0.nonzero_coefficient_indices(&r1);
            v.sort_unstable();
            v
        };
        let table_id = template.0.feasible_decompositions_handle().table_id();
        let a_plan = leaf_steady_plan_cached(&template, &mut dag0, &zero, &r1, &b1);
        let snap = get(table_id, &b1).expect("snapshot present");
        let snapshot_nz: Vec<Vec<usize>> = (*snap.nonzeros).clone();
        drop(dag0);

        let bits = |v: &[Ratio<i64>]| -> Vec<(i64, i64)> {
            v.iter().map(|r| (*r.numer(), *r.denom())).collect()
        };
        let mut results: Vec<Vec<(i64, i64)>> = Vec::new();
        std::thread::scope(|s| {
            let handles: Vec<_> = (0..8)
                .map(|_| {
                    let template = &template;
                    let a_plan = &a_plan;
                    let r1 = &r1;
                    let b1 = &b1;
                    s.spawn(move || {
                        let (_, mut dag, _) = leaf_plan_mk(3, 5);
                        assert!(leaf_steady_adopt_into(
                            &mut dag,
                            table_id,
                            a_plan,
                            b1
                        ));
                        let mut series = template.0.clone();
                        series.coefficients = vec![Ratio::<i64>::default(); basis_len];
                        let slices: Vec<&[Ratio<i64>]> = vec![r1.as_slice()];
                        dag.fold_batch(&mut series, &slices, a_plan, b1);
                        bits(&series.coefficients)
                    })
                })
                .collect();
            for h in handles {
                results.push(h.join().expect("adoption thread panicked"));
            }
        });
        for r in &results[1..] {
            assert_eq!(r, &results[0], "concurrent adoptions diverged");
        }

        // Ground truth: sequential folding of the same displacement.
        let mut seq = LogSignatureBuilder::<u8>::new()
            .with_num_dimensions(3)
            .with_max_degree(5)
            .build::<Ratio<i64>>();
        seq.concatenate_assign_coefficients(&r1);
        assert_eq!(
            results[0],
            bits(&seq.series.coefficients),
            "concurrent adoption folds diverged from sequential"
        );

        let snap_after = get(table_id, &b1).expect("snapshot survives");
        assert!(Arc::ptr_eq(&snap, &snap_after));
        assert_eq!(*snap_after.nonzeros, snapshot_nz, "snapshot mutated under concurrency");
    }


    /// Copy-on-write trigger: a dag holding ADOPTED (shared) lists that
    /// faces a support change must rebuild onto its OWN Arc — the shared
    /// snapshot stays pristine — and the rebuilt lists must equal a fresh
    /// uncached derivation's fixed point.
    #[test]

    fn leaf_plan_lists_cow_on_support_change_leaves_snapshot_pristine() {
        clear_for_tests();
        let (template, mut dag, basis_len) = leaf_plan_mk(3, 5);
        let zero = vec![Ratio::<i64>::default(); basis_len];
        let r1 = leaf_plan_rhs(basis_len, &[0, 1], 0xF81);
        let r2 = leaf_plan_rhs(basis_len, &[2], 0xF82);
        let support = |r: &Vec<Ratio<i64>>| {
            let mut v = template.0.nonzero_coefficient_indices(r);
            v.sort_unstable();
            v
        };
        let b1 = support(&r1);
        let b2 = support(&r2);
        assert_ne!(b1, b2, "test needs distinct supports");
        let table_id = template.0.feasible_decompositions_handle().table_id();

        // Adopt b1's snapshot into dag.
        let _a1 = leaf_steady_plan_cached(&template, &mut dag, &zero, &r1, &b1);
        let snap1 = get(table_id, &b1).expect("snapshot present");
        assert!(
            Arc::ptr_eq(&dag.nonzeros, &snap1.nonzeros),
            "precondition: the dag must share the snapshot"
        );

        // Support change on the SAME dag: the fixed-point loop rebuilds
        // via cow — onto the dag's own Arc, never through the shared one.
        let a2 = leaf_steady_plan(&template, &mut dag, &zero, &r2, &b2);
        assert!(
            !Arc::ptr_eq(&dag.nonzeros, &snap1.nonzeros),
            "the rebuild must install a new Arc (cow), not mutate the shared snapshot"
        );
        let snap1_after = get(table_id, &b1).expect("snapshot survives");
        assert!(Arc::ptr_eq(&snap1, &snap1_after));

        // The rebuilt lists must equal a fresh uncached derivation.
        let (_, mut dag_f, _) = leaf_plan_mk(3, 5);
        let a2_fresh = leaf_steady_plan(&template, &mut dag_f, &zero, &r2, &b2);
        assert_eq!(a2, a2_fresh, "cow rebuild's plan diverged");
        let (nz_cow, _) = dag.steady_lists_snapshot();
        let (nz_fresh, _) = dag_f.steady_lists_snapshot();
        assert_eq!(*nz_cow, *nz_fresh, "cow rebuild's lists diverged");
    }


    /// Memory discipline at scale: at 4x8 (the ~5MB-snapshot regime) the
    /// adopt must be O(nodes) — the dag's lists ARE the snapshot's
    /// allocation (`Arc::ptr_eq`), with exactly the cache and the adopting
    /// dag holding it. A deep copy would show up as a fresh allocation
    /// and fail the aliasing assert.
    #[test]

    fn leaf_plan_lists_adopt_is_zero_copy_at_scale() {
        clear_for_tests();
        let (template, mut dag, basis_len) = leaf_plan_mk(4, 8);
        let zero = vec![Ratio::<i64>::default(); basis_len];
        let r1 = leaf_plan_rhs(basis_len, &[0, 1, 2, 3], 0xF91);
        let b1 = {
            let mut v = template.0.nonzero_coefficient_indices(&r1);
            v.sort_unstable();
            v
        };
        let table_id = template.0.feasible_decompositions_handle().table_id();
        let a_plan = leaf_steady_plan_cached(&template, &mut dag, &zero, &r1, &b1);
        let snap = get(table_id, &b1).expect("snapshot present");
        let total: usize = snap.nonzeros.iter().map(|v| v.len()).sum();
        assert!(total > 0, "test needs a nontrivial snapshot");
        drop(dag);

        let (_, mut dag2, _) = leaf_plan_mk(4, 8);
        assert!(leaf_steady_adopt_into(&mut dag2, table_id, &a_plan, &b1));
        assert!(
            Arc::ptr_eq(&dag2.nonzeros, &snap.nonzeros),
            "adopt must not copy the lists"
        );
        assert_eq!(
            Arc::strong_count(&snap.nonzeros),
            2,
            "exactly the cache and the adopting dag hold the lists"
        );
    }


    /// Kill-switch bit identity: cached (warm) and uncached (cleared)
    /// tournaments must be bit-identical in BOTH cohort switch states —
    /// the cache swaps planning work, never engine or ordering semantics.
    #[test]

    fn leaf_plan_cache_tournament_bit_identity_both_switch_states() {
        use ordered_float::NotNan;

        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(4)
            .build()
            .expect("pool");
        let (d, m) = (3usize, 5usize);
        let builder = LogSignatureBuilder::<u8>::new()
            .with_num_dimensions(d)
            .with_max_degree(m);
        let basis =
            lyndon_rs::lyndon::LyndonBasis::<u8>::new(d, lyndon_rs::lyndon::Sort::Lexicographical)
                .generate_basis(m);
        let mut seed = 0x1eaf_u64 + d as u64 * 733 + m as u64;
        let lcg = |seed: &mut u64| {
            *seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((*seed >> 33) % 17) as i64 - 8
        };
        // Uniform dense supports (the cache's target workload), 17 leaves.
        let rhss: Vec<Vec<NotNan<f64>>> = (0..17)
            .map(|_| {
                (0..basis.len())
                    .map(|k| {
                        if k < d {
                            let v = lcg(&mut seed);
                            NotNan::new(if v == 0 { 2 } else { v } as f64).unwrap()
                        } else {
                            NotNan::new(0.0).unwrap()
                        }
                    })
                    .collect()
            })
            .collect();
        let slices: Vec<&[NotNan<f64>]> = rhss.iter().map(|r| r.as_slice()).collect();

        let run = |off: bool, warm: bool| -> Vec<NotNan<f64>> {
            if !warm {
                clear_for_tests();
            }
            let sig = builder.build::<NotNan<f64>>();
            set_cohort_off(off);
            let reduced = pool.install(|| sig.tournament_reduce(&slices));
            set_cohort_off(false);
            reduced
        };

        // Warm the cache once so the warm runs exercise the adoption path.
        let warm_cohort = run(false, false);
        let warm_scalar = run(true, false);
        assert_bits_identical(
            &warm_cohort,
            &warm_scalar,
            "warm-cache cohort vs scalar tournament diverged",
        );
        let cold_cohort = run(false, true);
        let cold_scalar = run(true, true);
        assert_bits_identical(
            &warm_cohort,
            &cold_cohort,
            "warm-cache cohort diverged from cold-cache cohort",
        );
        assert_bits_identical(
            &warm_scalar,
            &cold_scalar,
            "warm-cache scalar diverged from cold-cache scalar",
        );

        // Correctness anchor: the sequential fold chain. Float tournament
        // results differ from sequential folding by documented leaf-internal
        // reassociation ulps, so the sequential anchor runs on EXACT
        // rationals (same displacements, same tree) where any real planning
        // divergence — cache-related or not — is a value difference.
        {
            let rhss_r: Vec<Vec<Ratio<i64>>> = rhss
                .iter()
                .map(|r| {
                    r.iter()
                        // All LCG values are small integers (exact in f64),
                        // so the conversion is exact.
                        .map(|c| Ratio::from_integer(c.into_inner() as i64))
                        .collect()
                })
                .collect();
            let mut seq = builder.build::<Ratio<i64>>();
            for r in &rhss_r {
                seq.concatenate_assign_coefficients(r);
            }
            let sig_r = builder.build::<Ratio<i64>>();
            let slices_r: Vec<&[Ratio<i64>]> = rhss_r.iter().map(|r| r.as_slice()).collect();
            let tournament_r = pool.install(|| sig_r.tournament_reduce(&slices_r));
            let bits = |v: &[Ratio<i64>]| -> Vec<(i64, i64)> {
                v.iter().map(|r| (*r.numer(), *r.denom())).collect()
            };
            assert_eq!(
                bits(&tournament_r),
                bits(&seq.series.coefficients),
                "exact-arithmetic tournament (warm cache) diverged from sequential folding"
            );
        }
    }
}
