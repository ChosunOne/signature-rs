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
pub(super) fn clear_for_tests() {
    if let Some(map) = LEAF_STEADY_PLANS.get() {
        map.lock().expect("leaf steady plan cache mutex").clear();
    }
}
