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
