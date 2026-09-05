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

#[cfg(test)]
mod tests {
    use super::*;
    use commutator_rs::CommutatorTerm;
    use lie_rs::bch_series_generator::BchSeriesGenerator;
    use lie_rs::LieSeriesGenerator;
    use lyndon_rs::lyndon::{LyndonBasis, Sort};
    use num_rational::Ratio;
    use num_traits::{One, Zero};

    type R = Ratio<i128>;

    fn lcg(seed: &mut u64) -> i128 {
        *seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((*seed >> 33) % 19) as i128 - 9
    }

    fn atom(atom: u8) -> CommutatorTerm<R, u8> {
        CommutatorTerm::Atom {
            coefficient: R::one(),
            atom,
        }
    }

    fn bracket(left: CommutatorTerm<R, u8>, right: CommutatorTerm<R, u8>) -> CommutatorTerm<R, u8> {
        let degree = left.degree() + 1;
        CommutatorTerm::Expression {
            coefficient: R::one(),
            left: Box::new(left),
            right: Box::new(right),
            degree,
        }
    }

    /// Independent mark-and-sweep over a built DAG: the set of internal
    /// nodes reachable from the kept terms' roots. Mirrors
    /// `strip_dead_nodes`' marking phase so the tests do not trust the code
    /// under test to define its own oracle.
    fn reachable_roots(dag: &CommutatorDag<R>) -> Vec<bool> {
        let nodes = &dag.structure.nodes;
        let mut live = vec![false; nodes.len()];
        live[0] = true;
        live[1] = true;
        for (source, _) in dag.structure.terms.iter() {
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
        live
    }

    /// Census across the m grid: the BCH DAG is built from the 2-letter
    /// alphabet {A, B} regardless of the signature's dimension d, so the
    /// m axis spans every harness configuration (2x8/3x8/4x8 all build the
    /// m=8 DAG). Expected: zero dead nodes — construction filters
    /// zero-weight terms before interning and materializes only kept
    /// subtrees. This is the receipts for `strip_dead_nodes`' fast path.
    #[test]
    #[ignore = "expensive BCH builds across the m grid; run with --ignored --nocapture"]
    fn dag_dead_node_census() {
        for m in 3..=8usize {
            let bch_basis = LyndonBasis::<u8>::new(2, Sort::Lexicographical);
            let series: LieSeries<u8, R> =
                BchSeriesGenerator::new(bch_basis, m).generate_lie_series();
            let dag = CommutatorDag::<R>::from_bch_series(&series);
            let live = reachable_roots(&dag);
            let dead: Vec<usize> = (2..dag.structure.nodes.len())
                .filter(|&i| !live[i])
                .collect();
            let level_nodes: usize = dag.structure.levels.iter().map(|l| l.len()).sum();
            eprintln!(
                "m={m}: nodes={} internal={} terms={} levels={} level_nodes={} dead={}",
                dag.structure.nodes.len(),
                dag.structure.nodes.len() - 2,
                dag.structure.terms.len(),
                dag.structure.levels.len(),
                level_nodes,
                dead.len()
            );
            assert!(
                dead.is_empty(),
                "m={m}: {} dead nodes: {:?}",
                dead.len(),
                &dead[..dead.len().min(10)]
            );
        }
    }

    /// A dead node that *references* live children and a live node sitting
    /// after a hole: the strip must remove exactly the unreachable nodes,
    /// keep live children alive (they may be shared with live terms),
    /// renumber survivors in ascending id order (shifting ids past removed
    /// holes), remap every term source, and be a no-op on the second call.
    #[test]
    fn strip_dead_nodes_removes_unreachable_subtree_and_renumbers() {
        // 2: [A,B] live root | 3: [[A,B],B] dead (references live 2)
        // 4: [[A,B],A] live root | 5: [4-node, 3-node] dead (references
        // live 4 and dead 3).
        let mut nodes = vec![
            DagNode::Atom(0),
            DagNode::Atom(1),
            DagNode::Binary {
                left: 0,
                right: 1,
            },
            DagNode::Binary {
                left: 2,
                right: 1,
            },
            DagNode::Binary {
                left: 2,
                right: 0,
            },
            DagNode::Binary {
                left: 4,
                right: 3,
            },
        ];
        let mut terms = vec![
            (TermSource::Node(2), R::from_integer(3)),
            (TermSource::Node(4), R::from_integer(-7)),
        ];

        let removed = strip_dead_nodes(&mut nodes, &mut terms);
        assert_eq!(removed, 2, "nodes 3 and 5 are unreachable");
        assert_eq!(
            nodes,
            vec![
                DagNode::Atom(0),
                DagNode::Atom(1),
                DagNode::Binary {
                    left: 0,
                    right: 1
                },
                DagNode::Binary {
                    left: 2,
                    right: 0
                },
            ],
            "survivors renumbered ascending: old 4 became id 3"
        );
        assert_eq!(
            terms,
            vec![
                (TermSource::Node(2), R::from_integer(3)),
                (TermSource::Node(3), R::from_integer(-7)),
            ],
            "term sources remapped to the compact ids"
        );

        // Second call: everything live, fast path, nothing touched.
        let removed_again = strip_dead_nodes(&mut nodes, &mut terms);
        assert_eq!(removed_again, 0, "no dead nodes after the first strip");
        assert_eq!(nodes.len(), 4);
    }

    /// The construction-level regression guard: a series carrying a
    /// *zero-weight* deep bracket subtree and duplicate terms must compile
    /// to exactly the same DAG as the clean series — the zero-weight filter
    /// runs before interning (no dead subtree may leak into the node list)
    /// and duplicates reuse the interned root (no re-evaluation) — and the
    /// kept-term list grows by exactly the duplicates, preserving order.
    #[test]
    fn from_bch_series_zero_weight_and_duplicate_terms_leave_no_dead_nodes() {
        let bch_basis = LyndonBasis::<u8>::new(2, Sort::Lexicographical);
        let series: LieSeries<u8, R> =
            BchSeriesGenerator::new(bch_basis, 3).generate_lie_series();
        let baseline = CommutatorDag::<R>::from_bch_series(&series);

        // The last kept term must be an internal node so the duplicate
        // exercises the Node-source remap/intern-reuse path.
        let last_kept = baseline
            .structure
            .terms
            .last()
            .expect("m=3 series keeps terms")
            .clone();
        let last_root = match last_kept.0 {
            TermSource::Node(id) => id,
            other => panic!("expected an internal-node root, got {other:?}"),
        };
        let last_idx = (0..series.commutator_basis.len())
            .rev()
            .find(|&i| !series.coefficients[i].is_zero())
            .expect("non-zero term exists");

        // Pathological series: original terms + duplicate of the B atom +
        // duplicate of the last kept internal term + a deep zero-weight
        // subtree that would leave six dead nodes if it were ever interned.
        let mut modified = series.clone();
        let mut trees = series.commutator_basis.to_vec();
        let mut weights = series.coefficients.clone();
        trees.push(series.commutator_basis[1].clone()); // atom B duplicate
        weights.push(series.coefficients[1].clone());
        trees.push(series.commutator_basis[last_idx].clone()); // term duplicate
        weights.push(series.coefficients[last_idx].clone());
        let mut deep = atom(0);
        for _ in 0..6 {
            deep = bracket(deep, atom(1));
        }
        trees.push(deep); // zero weight: must be dropped BEFORE interning
        weights.push(R::zero());
        modified.commutator_basis = std::sync::Arc::new(trees);
        modified.coefficients = weights;

        // The same series WITHOUT the zero-weight subtree: folding it must
        // be indistinguishable from folding the pathological one — that is
        // the whole point of filtering before interning. (The duplicates
        // themselves are legitimately kept terms: they contribute twice,
        // so the fold is compared against this dup-only series, not the
        // clean baseline.)
        let mut dup_only = series.clone();
        let mut trees2 = series.commutator_basis.to_vec();
        let mut weights2 = series.coefficients.clone();
        trees2.push(series.commutator_basis[1].clone());
        weights2.push(series.coefficients[1].clone());
        trees2.push(series.commutator_basis[last_idx].clone());
        weights2.push(series.coefficients[last_idx].clone());
        dup_only.commutator_basis = std::sync::Arc::new(trees2);
        dup_only.coefficients = weights2;
        let dag_dup_only = CommutatorDag::<R>::from_bch_series(&dup_only);

        let dag = CommutatorDag::<R>::from_bch_series(&modified);

        // Bit-identical node structure; no dead nodes; terms grow by the
        // two duplicates in order, reusing the interned roots.
        assert_eq!(
            dag.structure.nodes, baseline.structure.nodes,
            "zero-weight subtree must not create nodes"
        );
        let live = reachable_roots(&dag);
        assert!(
            (2..dag.structure.nodes.len()).all(|i| live[i]),
            "every internal node must be reachable from a kept term root"
        );
        let mut expected_terms = baseline.structure.terms.clone();
        expected_terms.push((TermSource::Displacement, series.coefficients[1].clone()));
        expected_terms.push((TermSource::Node(last_root), series.coefficients[last_idx].clone()));
        assert_eq!(dag.structure.terms, expected_terms);

        // Evaluation through both DAGs must be bit-identical.
        let mut seed = 0xbeef_u64;
        let n = series.coefficients.len();
        let a: Vec<R> = (0..n)
            .map(|_| R::from_integer(lcg(&mut seed)))
            .collect();
        let b: Vec<R> = (0..n)
            .map(|k| {
                if k < 2 {
                    R::from_integer(lcg(&mut seed))
                } else {
                    R::from_integer(0)
                }
            })
            .collect();
        let a_nz: Vec<usize> = a
            .iter()
            .enumerate()
            .filter(|(_, c)| !c.is_zero())
            .map(|(i, _)| i)
            .collect();
        let b_nz: Vec<usize> = b
            .iter()
            .enumerate()
            .filter(|(_, c)| !c.is_zero())
            .map(|(i, _)| i)
            .collect();
        let (mut acc_base, mut acc_mod, mut acc_dup) = (
            vec![R::zero(); n],
            vec![R::zero(); n],
            vec![R::zero(); n],
        );
        let mut base = baseline.clone_shallow();
        let mut modified_dag = dag.clone_shallow();
        let mut dup_only_dag = dag_dup_only.clone_shallow();
        base.evaluate(&series, &a, &a_nz, &b, &b_nz);
        base.accumulate_terms(&mut acc_base, &b);
        modified_dag.evaluate(&series, &a, &a_nz, &b, &b_nz);
        modified_dag.accumulate_terms(&mut acc_mod, &b);
        dup_only_dag.evaluate(&series, &a, &a_nz, &b, &b_nz);
        dup_only_dag.accumulate_terms(&mut acc_dup, &b);
        assert_ne!(
            acc_base, acc_dup,
            "duplicates are kept terms: each copy must contribute"
        );
        assert_eq!(
            acc_mod, acc_dup,
            "zero-weight subtree must not perturb the fold — bit-identical to the dup-only series"
        );
    }
}
