use super::{
    AddAssign, Any, Arc, BchSeriesGenerator, CommutatorDag, Div, FeasibleDecompositions,
    FromPrimitive, Generator, HashMap, Hash, LieSeries, LieSeriesGenerator, LyndonBasis,
    LyndonWord, MulAssign, Mutex, Neg, OnceLock, One, Sort, SubAssign, TypeId, Zero,
};

/// Everything the basis and its structure-constant table depend on.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
struct SeriesPlanKey {
    generators: TypeId,
    coefficients: TypeId,
    alphabet_size: usize,
    max_degree: usize,
    sort: u8,
}

fn sort_tag(sort: Sort) -> u8 {
    match sort {
        Sort::Lexicographical => 0,
        Sort::Topological => 1,
    }
}

/// The shared plan for one configuration. Both members are read-only
/// after construction; `basis` is shared into every series so repeated
/// builds also skip `generate_basis`.
pub(super) struct SeriesPlan<T, U> {
    pub(super) basis: Arc<Vec<LyndonWord<T>>>,
    pub(super) table: Arc<FeasibleDecompositions<U>>,
}

static SERIES_PLANS: OnceLock<Mutex<HashMap<SeriesPlanKey, Arc<dyn Any + Send + Sync>>>> =
    OnceLock::new();

/// Returns the cached plan for this builder configuration, building it
/// once on first use. Callers must hold `T: 'static + Send + Sync` and
/// `U: 'static + Send + Sync` (the `TypeId` key and the shared cache
/// require it).
pub(super) fn series_plan<T, U>(
    basis_config: &LyndonBasis<T>,
    max_degree: usize,
) -> Arc<SeriesPlan<T, U>>
where
    T: Clone + Ord + Generator + Hash + Eq + Send + Sync + 'static,
    U: Clone
        + One
        + Zero
        + Eq
        + MulAssign
        + Neg<Output = U>
        + Hash
        + AddAssign
        + Ord
        + Send
        + Sync
        + 'static,
{
    let key = SeriesPlanKey {
        generators: TypeId::of::<T>(),
        coefficients: TypeId::of::<U>(),
        alphabet_size: basis_config.alphabet_size,
        max_degree,
        sort: sort_tag(basis_config.sort),
    };
    let map = SERIES_PLANS.get_or_init(|| Mutex::new(HashMap::new()));
    {
        let guard = map.lock().expect("series plan cache mutex");
        if let Some(hit) = guard.get(&key) {
            if let Ok(plan) = Arc::clone(hit).downcast::<SeriesPlan<T, U>>() {
                return plan;
            }
        }
    }
    // Cold path outside the lock: identical inputs produce an identical
    // deterministic plan, so a lost race only wastes one build.
    let basis = Arc::new(basis_config.generate_basis(max_degree));
    let table = Arc::new(LieSeries::<T, U>::build_feasible_decompositions(&basis));
    let plan = Arc::new(SeriesPlan {
        basis: Arc::clone(&basis),
        table,
    });
    let mut guard = map.lock().expect("series plan cache mutex");
    let entry = guard
        .entry(key)
        .or_insert_with(|| plan as Arc<dyn Any + Send + Sync>);
    entry
        .clone()
        .downcast::<SeriesPlan<T, U>>()
        .expect("plan type keyed by its TypeId pair")
}

/// The shared per-`(max_degree, U)` BCH plan: the 2-letter BCH Lie
/// series and a pristine template DAG. Both are pure functions of
/// `(max_degree, U)`; the series is cloned per build (cheap: Arc-shared
/// internals, fresh coefficient vector) and the DAG is cloned via
/// `clone_shallow` (shares the compiled structure and the node lists'
/// fixed point; per-value fold scratch — buffers, gating cache — stays
/// with each clone, and the template itself is never folded on).
pub(super) struct BchPlan<U> {
    pub(super) series: LieSeries<u8, U>,
    pub(super) dag: CommutatorDag<U>,
}

static BCH_PLANS: OnceLock<Mutex<HashMap<(TypeId, usize), Arc<dyn Any + Send + Sync>>>> =
    OnceLock::new();

pub(super) fn bch_plan<U>(max_degree: usize) -> Arc<BchPlan<U>>
where
    U: Clone
        + Default
        + AddAssign
        + Div<Output = U>
        + FromPrimitive
        + Neg<Output = U>
        + One
        + Ord
        + Hash
        + Send
        + Sync
        + Zero
        + SubAssign
        + MulAssign
        + 'static,
{
    let key = (TypeId::of::<U>(), max_degree);
    let map = BCH_PLANS.get_or_init(|| Mutex::new(HashMap::new()));
    {
        let guard = map.lock().expect("bch plan cache mutex");
        if let Some(hit) = guard.get(&key) {
            if let Ok(plan) = Arc::clone(hit).downcast::<BchPlan<U>>() {
                return plan;
            }
        }
    }
    let bch_basis = LyndonBasis::<u8>::new(2, Sort::Lexicographical);
    let series: LieSeries<u8, U> =
        BchSeriesGenerator::new(bch_basis, max_degree).generate_lie_series();
    let dag = CommutatorDag::from_bch_series(&series);
    let plan = Arc::new(BchPlan { series, dag });
    let mut guard = map.lock().expect("bch plan cache mutex");
    let entry = guard
        .entry(key)
        .or_insert_with(|| plan as Arc<dyn Any + Send + Sync>);
    entry
        .clone()
        .downcast::<BchPlan<U>>()
        .expect("plan type keyed by its TypeId")
}
#[cfg(test)]
mod tests {
    //! Shared structure-constant table (series plan cache) adversarial
    //! tests: config-scoped sharing, lazy-init races, read-only sharing.

    use super::*;
    use crate::{LogSignature, LogSignatureBuilder};
    use lyndon_rs::lyndon::{LyndonBasis, Sort};
    use ndarray::Array2;
    use ordered_float::NotNan;

    type F = NotNan<f64>;


    /// Deterministic small displacement path (LCG, no external rng).

    fn lcg_path(d: usize, n: usize, seed: u64) -> Array2<F> {
        let mut s = seed;
        let mut lcg = || {
            s = s.wrapping_mul(6364136223846793005).wrapping_add(1);
            ((s >> 33) % 11) as i64 - 5
        };
        Array2::from_shape_fn((n, d), |_| NotNan::new(lcg() as f64).unwrap())
    }



    fn coeff_bits(sig: &LogSignature<u8, F>) -> Vec<u64> {
        sig.series
            .coefficients
            .iter()
            .map(|c| c.into_inner().to_bits())
            .collect()
    }



    fn dec_bits(
        series: &LieSeries<u8, F>,
        i: usize,
        j: usize,
    ) -> Option<(Vec<u32>, Vec<u64>, bool)> {
        series
            .decomposition(i, j)
            .map(|(idx, c, swapped)| {
                (
                    idx.to_vec(),
                    c.iter().map(|x| x.into_inner().to_bits()).collect(),
                    swapped,
                )
            })
    }


    /// (a) The cache must share one table only within an identical
    /// configuration: same (d, m, sort, generator type, coefficient type)
    /// share one `Arc` (same `table_id`), while any differing axis gets its
    /// own entry — and every cached plan must still match a fresh
    /// `LieSeries::new` over the same basis (table size, here).
    #[test]

    fn plan_cache_shares_within_config_and_never_across_configs() {
        let b2m4 = LogSignatureBuilder::<u8>::new()
            .with_num_dimensions(2)
            .with_max_degree(4);
        let b3m4 = LogSignatureBuilder::<u8>::new()
            .with_num_dimensions(3)
            .with_max_degree(4);
        let b2m3 = LogSignatureBuilder::<u8>::new()
            .with_num_dimensions(2)
            .with_max_degree(3);
        let s2a = b2m4.build::<F>();
        let s2b = b2m4.build::<F>();
        let s3 = b3m4.build::<F>();
        let s2m3 = b2m3.build::<F>();
        assert_eq!(
            s2a.series.table_id(),
            s2b.series.table_id(),
            "same configuration must share one table Arc"
        );
        assert_ne!(
            s2a.series.table_id(),
            s3.series.table_id(),
            "different alphabet size must not share"
        );
        assert_ne!(
            s2a.series.table_id(),
            s2m3.series.table_id(),
            "different max_degree must not share"
        );
        for (built, (d, m)) in [(&s2a, (2usize, 4usize)), (&s3, (3, 4)), (&s2m3, (2, 3))] {
            let basis = LyndonBasis::<u8>::new(d, Sort::Lexicographical).generate_basis(m);
            let fresh =
                LieSeries::<u8, F>::new(basis.clone(), vec![F::default(); basis.len()]);
            assert_eq!(
                fresh.feasible_table_len(),
                built.series.feasible_table_len(),
                "cached table for ({d},{m}) must match a fresh build"
            );
        }
    }


    /// (b) Lazy-init race: 8 threads building/folding concurrently on ONE
    /// builder must all observe the same table contents and produce
    /// bit-identical results to the serial reference build. The folds run
    /// inside a 1-thread pool (parallel pack scheduling legitimately
    /// reorders float accumulation; the cache race under test is the
    /// concurrent `series_plan` init from the 8 outer threads).
    #[test]

    fn plan_cache_concurrent_builds_are_bit_identical() {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(1)
            .build()
            .expect("1-thread pool");
        let b = std::sync::Arc::new(
            LogSignatureBuilder::<u8>::new()
                .with_num_dimensions(3)
                .with_max_degree(5),
        );
        let path = lcg_path(3, 40, 0x51DE);
        let reference = pool.install(|| b.build_from_path(&path.view()));
        let ref_coeffs = coeff_bits(&reference);
        let ref_pair01 = dec_bits(&reference.series, 0, 1);
        let ref_pair12 = dec_bits(&reference.series, 1, 2);
        std::thread::scope(|scope| {
            for t in 0..8u64 {
                let b = &b;
                let path = &path;
                let pool = &pool;
                let ref_coeffs = &ref_coeffs;
                let ref_pair01 = &ref_pair01;
                let ref_pair12 = &ref_pair12;
                scope.spawn(move || {
                    let sig = pool.install(|| b.build_from_path(&path.view()));
                    assert_eq!(
                        coeff_bits(&sig),
                        *ref_coeffs,
                        "thread {t}: e2e coefficients diverged from serial reference"
                    );
                    assert_eq!(
                        dec_bits(&sig.series, 0, 1),
                        ref_pair01.clone(),
                        "thread {t}: pair (0,1) table contents diverged"
                    );
                    assert_eq!(
                        dec_bits(&sig.series, 1, 2),
                        ref_pair12.clone(),
                        "thread {t}: pair (1,2) table contents diverged"
                    );
                    let empty = b.build::<F>();
                    assert_eq!(
                        empty.series.table_id(),
                        sig.series.table_id(),
                        "thread {t}: empty build must hit the same cached table"
                    );
                    let folded = pool.install(|| empty.concatenate(&sig));
                    assert_eq!(
                        coeff_bits(&folded),
                        coeff_bits(&sig),
                        "thread {t}: concatenating onto the empty signature must be the identity"
                    );
                });
            }
        });
    }


    /// (c) The shared table is read-only: interleaved folds through one
    /// builder (different supports, different fold histories) must not
    /// change later folds, and everything must match a builder that never
    /// shared anything. All folds run in one 1-thread pool so parallel
    /// scheduling cannot reorder float accumulation.
    #[test]

    fn plan_cache_shared_table_leaks_no_fold_state() {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(1)
            .build()
            .expect("1-thread pool");
        let b = LogSignatureBuilder::<u8>::new()
            .with_num_dimensions(2)
            .with_max_degree(4);
        let pa = lcg_path(2, 24, 0xA);
        let pb = lcg_path(2, 24, 0xB);
        let a1 = pool.install(|| b.build_from_path(&pa.view()));
        let c = pool.install(|| b.build_from_path(&pb.view()));
        let a2 = pool.install(|| b.build_from_path(&pa.view()));
        assert_eq!(
            coeff_bits(&a1),
            coeff_bits(&a2),
            "re-folding the same path after a different fold diverged"
        );
        let fresh_b = LogSignatureBuilder::<u8>::new()
            .with_num_dimensions(2)
            .with_max_degree(4);
        let fresh = pool.install(|| fresh_b.build_from_path(&pa.view()));
        let fresh_c = pool.install(|| fresh_b.build_from_path(&pb.view()));
        assert_eq!(
            coeff_bits(&a1),
            coeff_bits(&fresh),
            "cached-table fold must match a fold on a never-cached builder"
        );
        // Explicit concatenation folds through the shared-table series.
        let r1 = pool.install(|| b.build::<F>().concatenate(&a1));
        let r1_again = pool.install(|| b.build::<F>().concatenate(&a1));
        assert_eq!(coeff_bits(&r1), coeff_bits(&r1_again));
        let r_cross = pool.install(|| b.build::<F>().concatenate(&c));
        let r_cross_fresh = pool.install(|| fresh_b.build::<F>().concatenate(&fresh_c));
        assert_eq!(
            coeff_bits(&r_cross),
            coeff_bits(&r_cross_fresh),
            "cross-path concatenation must match the fresh-builder fold"
        );
    }


    /// (d) The injected table is the true structure-constant fixed point:
    /// a series built through the normal `LieSeries::new` path (never
    /// touching the cache) must agree with the builder's cached table on
    /// EVERY feasible pair — indices, coefficient bits, and canonical
    /// orientation — while being a distinct instance.
    #[test]

    fn plan_cache_injected_table_matches_fresh_construction() {
        let (d, m) = (3usize, 4usize);
        let basis = LyndonBasis::<u8>::new(d, Sort::Lexicographical).generate_basis(m);
        let fresh = LieSeries::<u8, F>::new(basis.clone(), vec![F::default(); basis.len()]);
        let b = LogSignatureBuilder::<u8>::new()
            .with_num_dimensions(d)
            .with_max_degree(m);
        let built = b.build::<F>();
        assert_ne!(
            fresh.table_id(),
            built.series.table_id(),
            "the fresh path must not consume the builder's cache entry"
        );
        assert_eq!(fresh.feasible_table_len(), built.series.feasible_table_len());
        let degrees: Vec<usize> = basis.iter().map(|w| w.len()).collect();
        let mut checked = 0usize;
        for i in 0..basis.len() {
            for j in (i + 1)..basis.len() {
                if degrees[i] + degrees[j] > m {
                    continue;
                }
                let (fi, fc, fs) = fresh
                    .decomposition(i, j)
                    .unwrap_or_else(|| panic!("pair ({i},{j}) missing in fresh table"));
                let (bi, bc, bs) = built
                    .series
                    .decomposition(i, j)
                    .unwrap_or_else(|| panic!("pair ({i},{j}) missing in cached table"));
                assert!(!fs && !bs, "canonical pairs must not be flagged swapped");
                assert_eq!(fi, bi, "pair ({i},{j}) index mismatch");
                let fb: Vec<u64> = fc.iter().map(|c| c.into_inner().to_bits()).collect();
                let bb: Vec<u64> = bc.iter().map(|c| c.into_inner().to_bits()).collect();
                assert_eq!(fb, bb, "pair ({i},{j}) coefficient bits mismatch");
                checked += 1;
            }
        }
        assert!(checked > 10, "too few feasible pairs checked: {checked}");
    }
}
