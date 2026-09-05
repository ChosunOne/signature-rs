use super::structure::strip_dead_nodes;
use super::*;
use commutator_rs::CommutatorTerm;
use crate::LogSignatureBuilder;
use lie_rs::LieSeries;
use lyndon_rs::generators::ENotation;
use lyndon_rs::lyndon::{LyndonBasis, Sort};
use num_rational::Ratio;
use num_traits::{One, Zero};
use std::collections::HashMap;

type R = Ratio<i128>;

fn lcg(seed: &mut u64) -> i128 {
    *seed = seed
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    ((*seed >> 33) % 19) as i128 - 9
}

/// The pre-DAG recursive evaluator, kept as an independent oracle: same
/// tree walk, memoization and accumulation semantics as the old fold.
#[allow(clippy::too_many_arguments)]
fn reference_eval(
    term: &CommutatorTerm<R, u8>,
    series: &LieSeries<u8, R>,
    a: &[R],
    a_nonzero: &[usize],
    b: &[R],
    b_nonzero: &[usize],
    memo: &mut HashMap<CommutatorTerm<R, u8>, (Vec<R>, Vec<usize>)>,
) -> (Vec<R>, Vec<usize>) {
    match term {
        CommutatorTerm::Atom { atom, .. } => {
            if *atom == 0 {
                (a.to_vec(), a_nonzero.to_vec())
            } else {
                (b.to_vec(), b_nonzero.to_vec())
            }
        }
        CommutatorTerm::Expression { .. } => {
            if let Some(hit) = memo.get(term) {
                return hit.clone();
            }
            let (la, lnz) = reference_eval(
                term.left().unwrap(),
                series,
                a,
                a_nonzero,
                b,
                b_nonzero,
                memo,
            );
            let (rb, rnz) = reference_eval(
                term.right().unwrap(),
                series,
                a,
                a_nonzero,
                b,
                b_nonzero,
                memo,
            );
            let mut coefficients = vec![R::from_integer(0); a.len()];
            LieSeries::commutator_coefficients_with_nonzero(
                series,
                &la,
                &lnz,
                &rb,
                &rnz,
                &mut coefficients,
            );
            let nonzero = series.nonzero_coefficient_indices(&coefficients);
            memo.insert(term.clone(), (coefficients.clone(), nonzero.clone()));
            (coefficients, nonzero)
        }
    }
}

/// The DAG fold must reproduce the recursive evaluator exactly: same CSE,
/// same term weights, same accumulation order.
#[test]
fn dag_fold_matches_recursive_reference() {
    for (d, m) in [(2usize, 4usize), (2, 5), (3, 4)] {
        let mut seed = 0x5eed_u64.wrapping_mul(d as u64 * 31 + m as u64);
        let builder = LogSignatureBuilder::<u8>::new()
            .with_num_dimensions(d)
            .with_max_degree(m);

        // Accumulator: random dense coefficients.
        let mut log_sig = builder.build::<R>();
        let basis = LyndonBasis::<ENotation>::new(d, Sort::Lexicographical).generate_basis(m);
        let mut acc: Vec<R> = (0..basis.len())
            .map(|_| R::from_integer(lcg(&mut seed)))
            .collect();
        log_sig.series.coefficients.clone_from(&acc);

        // Reference state mirrors the accumulator and the fold inputs.
        let weights: Vec<R> = log_sig.bch_series.coefficients.clone();
        let trees: Vec<CommutatorTerm<R, u8>> = log_sig.bch_series.commutator_basis.to_vec();

        for step in 0..5 {
            // Perturb the accumulator between folds: zero a few
            // coefficients (support shrink) or raise a few zeros (support
            // growth), alternating with value-only changes. This drives
            // the DAG through both its collecting-rebuild and its
            // fixed-point-reuse paths, including the soundness corner
            // where a previously canceled scatter target goes non-zero.
            if step % 2 == 0 {
                for _ in 0..3 {
                    let k = (lcg(&mut seed) as usize) % acc.len();
                    acc[k] = R::from_integer(0);
                    log_sig.series.coefficients[k] = R::from_integer(0);
                }
            } else {
                for _ in 0..3 {
                    let k = d + (lcg(&mut seed) as usize) % (acc.len() - d);
                    let v = R::from_integer(lcg(&mut seed));
                    acc[k] = v.clone();
                    log_sig.series.coefficients[k] = v;
                }
            }

            // Production-shaped displacement: full-length coefficient
            // vector with only the degree-1 letters non-zero.
            let displacement: Vec<R> = (0..basis.len())
                .map(|k| {
                    if k < d {
                        R::from_integer(lcg(&mut seed))
                    } else {
                        R::from_integer(0)
                    }
                })
                .collect();

            // --- DAG fold (the code under test) ---
            let disp_sig = builder.build::<R>();
            let mut disp_sig = disp_sig;
            disp_sig.series.coefficients.clone_from(&displacement);
            log_sig.concatenate_assign(&disp_sig);

            // --- recursive reference fold ---
            // The atom-`A` input of every term is the *pre-fold*
            // accumulator snapshot (the old evaluator's
            // `original_coefficients`), not the mutating accumulator.
            let snapshot = acc.clone();
            let a_nonzero: Vec<usize> = snapshot
                .iter()
                .enumerate()
                .filter(|(_, c)| !c.is_zero())
                .map(|(i, _)| i)
                .collect();
            let b_nonzero: Vec<usize> = displacement
                .iter()
                .enumerate()
                .filter(|(_, c)| !c.is_zero())
                .map(|(i, _)| i)
                .collect();
            let mut memo = HashMap::new();
            for (i, tree) in trees.iter().enumerate() {
                if i == 0 || weights[i].is_zero() {
                    continue;
                }
                let (ct, _nz) = reference_eval(
                    tree,
                    &log_sig.series,
                    &snapshot,
                    &a_nonzero,
                    &displacement,
                    &b_nonzero,
                    &mut memo,
                );
                for (acc_coeff, comm_coeff) in acc.iter_mut().zip(&ct) {
                    *acc_coeff += comm_coeff * weights[i].clone();
                }
            }

            let diffs: Vec<_> = log_sig
                .series
                .coefficients
                .iter()
                .zip(&acc)
                .enumerate()
                .filter(|(_, (g, w))| g != w)
                .take(6)
                .map(|(k, (g, w))| format!("k={k} dag={g} ref={w}"))
                .collect();
            assert!(
                diffs.is_empty(),
                "d={d} m={m} step={step}: DAG fold diverged: {}",
                diffs.join("; ")
            );
        }
    }
}

/// Liveness of the one-dispatch fold under nested parallelism. When the
/// shared pool is saturated by outer tasks, a fold's slot tasks queue
/// behind them; a fold that blocked on a *queued* releaser's pack would
/// never finish (the fixed slot→pack assignment starved exactly there).
/// Pack claims go to whoever runs, so folds nested inside an
/// oversubscribed outer dispatch must complete — and match the serial
/// result bit for bit, since a claim relabels which slot sweeps a pack,
/// never the order of adds within a word.
#[test]
fn dag_fold_survives_nested_parallel_oversubscription() {
    use rayon::iter::{IntoParallelIterator, ParallelIterator};
    use std::sync::atomic::AtomicBool;

    let d = 3usize;
    let m = 6usize;
    let mut seed = 0x51ea_u64;
    let builder = LogSignatureBuilder::<u8>::new()
        .with_num_dimensions(d)
        .with_max_degree(m);
    let basis = LyndonBasis::<ENotation>::new(d, Sort::Lexicographical).generate_basis(m);

    // One accumulator and one displacement, folded once on the calling
    // thread to produce the expected post-fold coefficients.
    let acc: Vec<R> = (0..basis.len())
        .map(|_| R::from_integer(lcg(&mut seed)))
        .collect();
    let displacement: Vec<R> = (0..basis.len())
        .map(|k| {
            if k < d {
                R::from_integer(lcg(&mut seed))
            } else {
                R::from_integer(0)
            }
        })
        .collect();
    let mut template = builder.build::<R>();
    template.series.coefficients.clone_from(&acc);
    let mut disp_sig = builder.build::<R>();
    disp_sig.series.coefficients.clone_from(&displacement);
    let mut reference = template.clone();
    reference.concatenate_assign(&disp_sig);
    let expected = reference.series.coefficients.clone();

    // Watchdog: a liveness bug shows up as a hang, not a wrong value —
    // turn it into a hard failure (abort fails the test process).
    let done = std::sync::Arc::new(AtomicBool::new(false));
    let flag = done.clone();
    std::thread::spawn(move || {
        for _ in 0..600 {
            std::thread::sleep(std::time::Duration::from_millis(100));
            if flag.load(std::sync::atomic::Ordering::SeqCst) {
                return;
            }
        }
        eprintln!("nested fold stress did not complete in 60s: fold walk starved");
        std::process::abort();
    });

    // 2× oversubscription: 64 nested folds on the default pool.
    (0..64usize).into_par_iter().for_each(|_| {
        let mut sig = template.clone();
        sig.concatenate_assign(&disp_sig);
        assert_eq!(sig.series.coefficients, expected);
    });
    done.store(true, std::sync::atomic::Ordering::SeqCst);
}

// --- dead-node strip: census and adversarial coverage ---

use lie_rs::bch_series_generator::BchSeriesGenerator;
use lie_rs::LieSeriesGenerator;

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
        let series: LieSeries<u8, R> = BchSeriesGenerator::new(bch_basis, m).generate_lie_series();
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

/// The stripped hand DAG is executable: evaluating the compacted
/// 4-node plan (atoms + [A,B] + [[A,B],A]) through the production
/// rebuild path must match the recursive oracle bit for bit.
#[test]
fn strip_then_evaluate_matches_recursive_oracle() {
    let d = 3usize;
    let m = 3usize;
    let mut seed = 0xdead_u64;
    let builder = LogSignatureBuilder::<u8>::new()
        .with_num_dimensions(d)
        .with_max_degree(m);
    let basis = LyndonBasis::<ENotation>::new(d, Sort::Lexicographical).generate_basis(m);

    // The compacted DAG from the renumbering test: [A,B] and [[A,B],A].
    let mut dag = CommutatorDag::<R> {
        structure: std::sync::Arc::new(DagStructure {
            nodes: vec![
                DagNode::Atom(0),
                DagNode::Atom(1),
                DagNode::Binary {
                    left: 0,
                    right: 1,
                },
                DagNode::Binary {
                    left: 2,
                    right: 0,
                },
            ],
            levels: vec![vec![0, 1], vec![2], vec![3]],
            terms: vec![
                (TermSource::Node(2), R::from_integer(3)),
                (TermSource::Node(3), R::from_integer(-7)),
            ],
        }),
        buffers: Vec::new(),
        nonzeros: Arc::new(Vec::new()),
        dirty: Arc::new(Vec::new()),
        atom_a: Vec::new(),
        atom_b: Vec::new(),
        lists_built: false,
        gating_cache: GatingCache::default(),
        class_order: OnceLock::new(),
    };

    let a: Vec<R> = (0..basis.len())
        .map(|_| R::from_integer(lcg(&mut seed)))
        .collect();
    let b: Vec<R> = (0..basis.len())
        .map(|k| if k < d { R::from_integer(lcg(&mut seed)) } else { R::from_integer(0) })
        .collect();
    let a_nonzero: Vec<usize> = a
        .iter()
        .enumerate()
        .filter(|(_, c)| !c.is_zero())
        .map(|(i, _)| i)
        .collect();
    let b_nonzero: Vec<usize> = b
        .iter()
        .enumerate()
        .filter(|(_, c)| !c.is_zero())
        .map(|(i, _)| i)
        .collect();

    dag.evaluate(&builder.build::<R>().series, &a, &a_nonzero, &b, &b_nonzero);

    // Recursive oracle on the equivalent bracket trees. Node buffers
    // live in the DAG's class-contiguous order (the public-order
    // epilogue only runs in `accumulate_terms`), so the oracle's
    // public-order coefficients are gathered through the same
    // `class[inv[w]] = public[w]` map the fold's inputs use before
    // comparing.
    let t2 = bracket(atom(0), atom(1)); // node 2 = [A, B]
    let t3 = bracket(t2.clone(), atom(0)); // node 3 = [[A,B], A]
    let series = builder.build::<R>().series;
    let inv = dag
        .class_order
        .get()
        .expect("evaluate built the class order")
        .inv();
    let mut memo = HashMap::new();
    for (node, tree) in [(2u32, &t2), (3, &t3)] {
        let (ct, _nz) = reference_eval(
            tree,
            &series,
            &a,
            &a_nonzero,
            &b,
            &b_nonzero,
            &mut memo,
        );
        let mut cls = vec![R::zero(); ct.len()];
        for (w, value) in ct.iter().enumerate() {
            cls[inv[w] as usize] = value.clone();
        }
        assert_eq!(
            dag.buffers[node as usize - 2], cls,
            "node {node} buffer diverged from the recursive oracle"
        );
    }
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
    let series: LieSeries<u8, R> = BchSeriesGenerator::new(bch_basis, 3).generate_lie_series();
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
        .map(|k| if k < 2 { R::from_integer(lcg(&mut seed)) } else { R::from_integer(0) })
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
