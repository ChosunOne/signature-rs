use std::fmt::{Debug, Formatter};
use std::ops::{AddAssign, MulAssign, Neg};
use std::sync::Arc;
use std::{collections::HashMap, hash::Hash};

use commutator_rs::CommutatorTerm;
use lie_rs::LieSeries;
use lyndon_rs::Generator;
use num_traits::{One, Zero};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum DagNode {
    Atom(u8),
    Binary { left: u32, right: u32 },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum TermSource {
    Node(u32),
    Displacement,
}

struct DagStructure<U> {
    /// Topologically sorted: ids 0 and 1 are the atoms, every
    /// internal node's children have strictly smaller ids.
    nodes: Vec<DagNode>,
    /// `(source, BCH weight)` per accumulated term, in the
    /// original term order (float summation order preserved).
    terms: Vec<(TermSource, U)>,
}

pub(crate) struct CommutatorDag<U> {
    structure: Arc<DagStructure<U>>,
    buffers: Vec<Vec<U>>,
    nonzeros: Vec<Vec<usize>>,
}

/// Coefficient slice for a node `idx`: an atom reads the
/// fold's inputs, an internal node its own result buffer.
#[inline]
fn node_slice<'a, U>(idx: u32, before: &'a [Vec<U>], a: &'a [U], b: &'a [U]) -> &'a [U] {
    if idx < 2 {
        if idx == 0 { a } else { b }
    } else {
        &before[idx as usize - 2]
    }
}

/// Non-zero index list counterpart of [`node_slice`].
fn node_nonzeros<'a>(
    idx: u32,
    before: &'a [Vec<usize>],
    a_nonzero: &'a [usize],
    b_nonzero: &'a [usize],
) -> &'a [usize] {
    if idx < 2 {
        if idx == 0 { a_nonzero } else { b_nonzero }
    } else {
        &before[idx as usize - 2]
    }
}

impl<U: Clone + Zero> CommutatorDag<U> {
    pub(crate) fn from_bch_series(bch_series: &LieSeries<u8, U>) -> Self
    where
        U: Eq + Hash + One,
    {
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

        let mut nodes = vec![DagNode::Atom(0), DagNode::Atom(1)];
        let mut interned = HashMap::new();

        let mut terms = Vec::new();
        for (i, term) in bch_series.commutator_basis.iter().enumerate() {
            if i == 0 {
                continue;
            }
            let weight = &bch_series.coefficients[i];
            if weight.is_zero() {
                continue;
            }
            let root = intern(term, &mut nodes, &mut interned);
            let source = match root {
                0 => continue,
                1 => TermSource::Displacement,
                id => TermSource::Node(id),
            };
            terms.push((source, weight.clone()));
        }

        Self {
            structure: Arc::new(DagStructure { nodes, terms }),
            buffers: Vec::new(),
            nonzeros: Vec::new(),
        }
    }
}

impl<U: Clone + Default + One + Zero + Eq + MulAssign + Neg<Output = U> + Hash + AddAssign>
    CommutatorDag<U>
{
    pub(crate) fn evaluate<T: Clone + Ord + Generator + Hash + Eq>(
        &mut self,
        series: &LieSeries<T, U>,
        a: &[U],
        a_nonzero: &[usize],
        b: &[U],
        b_nonzero: &[usize],
    ) {
        let internal = self.structure.nodes.len() - 2;
        self.ensure_buffers(series.coefficients.len(), internal);

        for k in 2..self.structure.nodes.len() {
            let (left, right) = match self.structure.nodes[k] {
                DagNode::Binary { left, right } => (left, right),
                DagNode::Atom(_) => unreachable!("atoms are the first two nodes"),
            };
            let (before, rest) = self.buffers.split_at_mut(k - 2);
            let result = &mut rest[0];
            result.fill(U::default());

            let lbuf = node_slice(left, before, a, b);
            let rbuf = node_slice(right, before, a, b);
            let (nz_before, nz_rest) = self.nonzeros.split_at_mut(k - 2);
            let lnz = node_nonzeros(left, nz_before, a_nonzero, b_nonzero);
            let rnz = node_nonzeros(right, nz_before, a_nonzero, b_nonzero);

            LieSeries::commutator_coefficients_with_nonzero(series, lbuf, lnz, rbuf, rnz, result);

            nz_rest[0] = series.nonzero_coefficient_indices(result);
        }
    }

    pub(crate) fn terms(&self) -> impl Iterator<Item = (TermSource, &U)> {
        self.structure
            .terms
            .iter()
            .map(|(source, weight)| (*source, weight))
    }
    pub(crate) fn buffer(&self, node: u32) -> &[U] {
        &self.buffers[node as usize - 2]
    }

    fn ensure_buffers(&mut self, d: usize, count: usize) {
        if self.buffers.len() != count || self.buffers.first().is_none_or(|b| b.len() != d) {
            self.buffers = (0..count).map(|_| vec![U::default(); d]).collect();
            self.nonzeros = vec![Vec::new(); count];
        }
    }
}

impl<U> CommutatorDag<U> {
    pub(crate) fn clone_shallow(&self) -> Self {
        Self {
            structure: Arc::clone(&self.structure),
            buffers: Vec::new(),
            nonzeros: Vec::new(),
        }
    }
}

impl<U> Debug for CommutatorDag<U> {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("CommutatorDag")
            .field("nodes", &self.structure.nodes.len())
            .field("terms", &self.structure.terms.len())
            .field("buffers", &self.buffers.len())
            .finish_non_exhaustive()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{LogSignature, LogSignatureBuilder};
    use lyndon_rs::{
        LyndonWord,
        lyndon::{LyndonBasis, Sort},
    };
    use num_rational::Ratio;

    type R = Ratio<i128>;

    fn lcg(seed: &mut u64) -> i128 {
        *seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((*seed >> 33) % 19) as i128 - 9
    }

    /// The pre-DAG recursive evaluator, kept as an independent oracle: same
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
        match *term {
            CommutatorTerm::Atom { ref atom, .. } => {
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
                let mut coefficients: Vec<Ratio<i128>> = vec![R::from_integer(0); a.len()];
                LieSeries::commutator_coefficients_with_nonzero(
                    series,
                    &la,
                    &lnz,
                    &rb,
                    &rnz,
                    &mut coefficients,
                );
                let nonzero: Vec<usize> = series.nonzero_coefficient_indices(&coefficients);
                memo.insert(term.clone(), (coefficients.clone(), nonzero.clone()));
                (coefficients, nonzero)
            }
        }
    }

    #[test]
    fn dag_fold_matches_recursive_reference() {
        for (d, m) in [(2usize, 4usize), (2, 5), (3, 4)] {
            let mut seed: u64 = 0x5eed_u64.wrapping_mul(d as u64 * 31 + m as u64);
            let builder: LogSignatureBuilder<u8> = LogSignatureBuilder::<u8>::new()
                .with_num_dimensions(d)
                .with_max_degree(m);

            // Accumulator: random dense coefficients.
            let mut log_sig: LogSignature<u8, Ratio<i128>> = builder.build::<R>();
            let basis: Vec<LyndonWord<u8>> =
                LyndonBasis::<u8>::new(d, Sort::Lexicographical).generate_basis(m);
            let mut acc: Vec<R> = (0..basis.len())
                .map(|_| R::from_integer(lcg(&mut seed)))
                .collect();
            log_sig.series.coefficients.clone_from(&acc);

            // Reference state mirrors the accumulator and the fold inputs.
            let weights: Vec<R> = log_sig.bch_series.coefficients.clone();
            let trees: Vec<CommutatorTerm<R, u8>> = log_sig.bch_series.commutator_basis.to_vec();

            for step in 0..3 {
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

                // --- DAG fold (the code under test) ---   Not Committed Yet
                let disp_sig: LogSignature<u8, Ratio<i128>> = builder.build::<R>();
                let mut disp_sig: LogSignature<u8, Ratio<i128>> = disp_sig;
                disp_sig.series.coefficients.clone_from(&displacement);
                log_sig.concatenate_assign(&disp_sig);

                // --- recursive reference fold ---
                // The atom-`A` input of every term is the *pre-fold*
                // accumulator snapshot (the old evaluator's
                // `original_coefficients`), not the mutating accumulator.
                let snapshot: Vec<Ratio<i128>> = acc.clone();
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
                        *acc_coeff += comm_coeff * weights[i];
                    }
                }

                let diffs: Vec<_> = log_sig
                    .series
                    .coefficients
                    .iter()
                    .zip(&acc)
                    .enumerate()
                    .filter(|&(_, (g, w))| g != w)
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
}
