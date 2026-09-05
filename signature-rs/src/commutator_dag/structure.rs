//! The compiled DAG's node graph: atom and bracket nodes, the shared
//! structure tables, and the constructor that compiles a BCH series into
//! the plan (interning shared subtrees, stripping dead nodes).

use commutator_rs::CommutatorTerm;
use lie_rs::{GatingCache, LieSeries};
use num_traits::{One, Zero};
use std::collections::HashMap;
use std::hash::Hash;
use std::sync::{Arc, OnceLock};

use super::CommutatorDag;

/// One node of the compiled plan.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum DagNode {
    /// Leaf: `0` is the accumulator (`A`), `1` the displacement (`B`).
    Atom(u8),
    /// Internal bracket node; children have strictly smaller ids.
    Binary { left: u32, right: u32 },
}

/// Where a term's coefficient slice lives after evaluation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum TermSource {
    /// The term's root node result, held in the DAG's buffers.
    Node(u32),
    /// The term is the bare displacement atom `B`; its coefficients are the
    /// fold's `rhs` input slice.
    Displacement,
}

pub(crate) struct DagStructure<U> {
    /// Topologically sorted: ids 0 and 1 are the atoms, every internal node's
    /// children have strictly smaller ids. Invariant: every internal node is
    /// reachable from at least one kept term's root (enforced by
    /// [`strip_dead_nodes`] at construction), so the fold loop evaluates no
    /// dead nodes.
    pub(crate) nodes: Vec<DagNode>,
    /// Nodes grouped by dependency height (atoms at height 0, a bracket one
    /// above its tallest child). All nodes of a level depend only on earlier
    /// levels, so a level evaluates as one conflict-free parallel batch.
    pub(crate) levels: Vec<Vec<u32>>,
    /// `(source, BCH weight)` per accumulated term, in the original term
    /// order (float summation order preserved).
    pub(crate) terms: Vec<(TermSource, U)>,
}

/// Strips dead nodes from a freshly compiled DAG: removes every internal
/// node unreachable from some kept term's root, renumbers the survivors to
/// compact ids, and remaps the term sources.
///
/// A term's coefficients are read exactly at its root node's buffer (its
/// weight applies to that root result), and a root's value needs exactly
/// its descendants — so reachability from the kept roots is liveness. The
/// atoms are always live: every bracket evaluation reads them, and the
/// displacement term reads atom `B` through the fold input.
///
/// Survivors keep their relative id order, and a live node's children are
/// live with smaller ids (reachability is closed downward), so ascending
/// renumbering preserves the topological order the level tables rely on.
/// The buffers, scatter-target supersets and dirty bitsets are sized from
/// the node count at fold time and don't exist yet here — compacting the
/// node vector before those allocations is the whole rebuild.
///
/// Returns the number of removed nodes (0 = the mark sweep found every
/// internal node live and nothing was touched).
pub(super) fn strip_dead_nodes<U>(nodes: &mut Vec<DagNode>, terms: &mut [(TermSource, U)]) -> usize {
    // Mark: atoms plus the downward closure of every kept term root.
    let mut live = vec![false; nodes.len()];
    live[0] = true;
    live[1] = true;
    for (source, _) in terms.iter() {
        if let TermSource::Node(root) = *source {
            let mut stack = vec![root];
            while let Some(id) = stack.pop() {
                if live[id as usize] {
                    continue;
                }
                live[id as usize] = true;
                if let DagNode::Binary { left, right } = nodes[id as usize] {
                    stack.push(left);
                    stack.push(right);
                }
            }
        }
    }
    let internal = nodes.len() - 2;
    let live_internal = (2..nodes.len()).filter(|&i| live[i]).count();
    if live_internal == internal {
        return 0; // fast path: the construction produced no dead nodes
    }

    // Compact: ascending renumber (atoms keep 0 and 1 — set explicitly:
    // the zero-initialized vector would alias atom `B` onto atom `A`),
    // then one rebuild pass that drops dead subtrees and remaps live
    // children.
    let mut remap = vec![0u32; nodes.len()];
    remap[1] = 1;
    let mut next = 2u32;
    for (idx, _) in nodes.iter().enumerate().skip(2) {
        if live[idx] {
            remap[idx] = next;
            next += 1;
        }
    }
    let old = std::mem::take(nodes);
    nodes.reserve(next as usize);
    for (idx, node) in old.iter().enumerate() {
        match node {
            DagNode::Atom(atom) => nodes.push(DagNode::Atom(*atom)),
            DagNode::Binary { left, right } => {
                if live[idx] {
                    debug_assert!(live[*left as usize] && live[*right as usize]);
                    nodes.push(DagNode::Binary {
                        left: remap[*left as usize],
                        right: remap[*right as usize],
                    });
                }
            }
        }
    }
    for (source, _) in terms.iter_mut() {
        if let TermSource::Node(id) = source {
            debug_assert!(live[*id as usize]);
            *id = remap[*id as usize];
        }
    }
    internal - live_internal
}

impl<U> CommutatorDag<U>
where
    U: Clone + Zero,
{
    /// Compiles the plan from a BCH series' bracket trees.
    ///
    /// Only terms with non-zero BCH weights are kept; their relative order
    /// (and therefore the accumulation's float summation order) is preserved.
    /// The compiled node set is stripped of dead nodes before the level
    /// tables are built: every internal node of the finished DAG is
    /// reachable from some kept term's root (see [`strip_dead_nodes`]).
    pub(crate) fn from_bch_series(bch_series: &LieSeries<u8, U>) -> Self
    where
        U: Eq + Hash + One,
    {
        let mut nodes = vec![DagNode::Atom(0), DagNode::Atom(1)];
        // Intern map keyed by the structural `unit_hash`; the bucket vector
        // resolves (astronomically rare) hash collisions by structural
        // equality, so interning never clones a tree for a map key.
        let mut interned: HashMap<u64, Vec<(CommutatorTerm<U, u8>, u32)>> = HashMap::new();

        fn intern<U: Eq + Hash + Clone + One>(
            term: &CommutatorTerm<U, u8>,
            nodes: &mut Vec<DagNode>,
            interned: &mut HashMap<u64, Vec<(CommutatorTerm<U, u8>, u32)>>,
        ) -> u32 {
            match term {
                CommutatorTerm::Atom { atom, .. } => *atom as u32,
                CommutatorTerm::Expression { .. } => {
                    let hash = term.unit_hash();
                    if let Some(bucket) = interned.get(&hash) {
                        for (existing, id) in bucket {
                            if existing == term {
                                return *id;
                            }
                        }
                    }
                    let left = intern(
                        term.left().expect("expression has a left child"),
                        nodes,
                        interned,
                    );
                    let right = intern(
                        term.right().expect("expression has a right child"),
                        nodes,
                        interned,
                    );
                    let id = nodes.len() as u32;
                    nodes.push(DagNode::Binary { left, right });
                    interned.entry(hash).or_default().push((term.clone(), id));
                    id
                }
            }
        }

        let mut terms = Vec::new();
        for (i, term) in bch_series.commutator_basis.iter().enumerate() {
            if i == 0 {
                continue; // the degree-1 `A` term: its contribution is the accumulator
            }
            let weight = &bch_series.coefficients[i];
            if weight.is_zero() {
                continue;
            }
            let root = intern(term, &mut nodes, &mut interned);
            let source = match root {
                0 => continue, // defensive: an atom-`A` root adds nothing new
                1 => TermSource::Displacement,
                id => TermSource::Node(id),
            };
            terms.push((source, weight.clone()));
        }

        // Dead-node strip, before any id-indexed table is derived: the
        // height/level pass below, and every downstream structure (node
        // buffers, scatter-target supersets, dirty bitsets — all sized from
        // `nodes.len()` at fold time) must describe live nodes only.
        //
        // The dead set is empty for today's construction: zero-weight terms
        // are filtered *before* interning and `intern` materializes only
        // whole subtrees of kept roots, so every created node is reachable
        // from some kept term root. The strip stays as a cheap invariant
        // guard — one O(nodes + terms) mark-and-sweep with an early exit
        // when everything is live — so a future construction change (e.g.
        // interning before filtering, or partial-tree assembly) cannot
        // silently grow the per-fold evaluation loop with unreachable nodes.
        let _stripped = strip_dead_nodes(&mut nodes, &mut terms);

        // Dependency heights and levels: atoms at height 0, a bracket one
        // above its tallest child. The nodes are topologically sorted, so a
        // single forward pass suffices.
        let mut heights = vec![0u32; nodes.len()];
        let mut max_height = 0u32;
        for (idx, node) in nodes.iter().enumerate() {
            if let DagNode::Binary { left, right } = node {
                let h = 1 + heights[*left as usize].max(heights[*right as usize]);
                heights[idx] = h;
                max_height = max_height.max(h);
            }
        }
        let mut levels: Vec<Vec<u32>> = vec![Vec::new(); max_height as usize + 1];
        for (idx, h) in heights.iter().enumerate() {
            levels[*h as usize].push(idx as u32);
        }

        Self {
            structure: Arc::new(DagStructure {
                nodes,
                terms,
                levels,
            }),
            buffers: Vec::new(),
            nonzeros: Arc::new(Vec::new()),
            dirty: Arc::new(Vec::new()),
            atom_a: Vec::new(),
            atom_b: Vec::new(),
            lists_built: false,
            gating_cache: GatingCache::default(),
            class_order: OnceLock::new(),
        }
    }
}
