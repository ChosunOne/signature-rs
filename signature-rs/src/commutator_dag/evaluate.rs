//! Single-fold plan evaluation, batch eligibility and the public-order
//! term-accumulation epilogue.

use lie_rs::{ClassOrderedCommutation, KernelJob, LieSeries};
use lyndon_rs::generators::Generator;
use num_traits::{One, Zero};
use std::hash::Hash;
use std::ops::{AddAssign, MulAssign, Neg};
use std::sync::Arc;

use super::structure::{DagNode, TermSource};
use super::{CommutatorDag, node_nonzeros, node_slice};

impl<U> CommutatorDag<U>
where
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
{
    /// Evaluates the whole plan for one concatenation.
    ///
    /// `a` is the accumulated log-signature's coefficient slice *before* this
    /// fold's BCH additions, `b` the segment displacement; after the call the
    /// per-term results are available via [`Self::terms`]/[`Self::buffer`]
    /// for accumulation into the log signature.
    pub(crate) fn evaluate<T>(
        &mut self,
        series: &LieSeries<T, U>,
        a: &[U],
        a_nonzero: &[usize],
        b: &[U],
        b_nonzero: &[usize],
    ) where
        T: Clone + Ord + Generator + Hash + Eq,
        U: Send + Sync,
    {
        let internal = self.structure.nodes.len() - 2;
        self.ensure_buffers(series.coefficients.len(), internal);

        // The class-contiguous ordering: created on the first fold, shared
        // by every kernel call of every fold through this DAG. The fold's
        // whole working set — node buffers and support lists — lives in the
        // ordering; the single public-order epilogue runs in
        // `accumulate_terms`.
        let order = self.class_order.get_or_init(|| series.class_order());
        let inv = order.inv();
        // The fold's inputs arrive in public basis order: convert once per
        // fold. (The rebuild trigger below keeps comparing the caller's
        // public support lists, so it is unaffected by the relabeling.)
        let a_cls = series.class_coefficients(order, a);
        let b_cls = series.class_coefficients(order, b);
        let a_nz_cls: Vec<usize> = a_nonzero.iter().copied().map(|i| inv[i] as usize).collect();
        let b_nz_cls: Vec<usize> = b_nonzero.iter().copied().map(|i| inv[i] as usize).collect();

        // The node lists are a fixed point of the DAG's support propagation:
        // they change only when an atom support changes. Value-level changes
        // with an unchanged support leave them valid; a support change forces
        // one collecting pass in which the kernel re-derives every node's
        // scatter-target set top-down (children before parents).
        let rebuild = !self.lists_built || self.atom_a != a_nonzero || self.atom_b != b_nonzero;

        if rebuild {
            // Copy-on-write: a dag sharing an adopted snapshot's lists must
            // never mutate them through the `Arc`. Unwrap when exclusively
            // held (free), clone when shared, mutate the owned copy, and
            // install it as this dag's own `Arc`.
            let mut nz = Arc::unwrap_or_clone(std::mem::take(&mut self.nonzeros));
            let mut dt = Arc::unwrap_or_clone(std::mem::take(&mut self.dirty));
            // One serial topological sweep in which the kernel re-derives
            // every node's scatter-target set: topological order guarantees
            // a node's children's lists were rebuilt earlier in the pass.
            for k in 2..self.structure.nodes.len() {
                let (left, right) = match self.structure.nodes[k] {
                    DagNode::Binary { left, right } => (left, right),
                    DagNode::Atom(_) => unreachable!("atoms are the first two nodes"),
                };
                let (before, rest) = self.buffers.split_at_mut(k - 2);
                let result = &mut rest[0];
                result.fill(U::default());
                let lbuf = node_slice(left, before, &a_cls, &b_cls);
                let rbuf = node_slice(right, before, &a_cls, &b_cls);
                let (nz_before, nz_rest) = nz.split_at_mut(k - 2);
                let lnz = node_nonzeros(left, nz_before, &a_nz_cls, &b_nz_cls);
                let rnz = node_nonzeros(right, nz_before, &a_nz_cls, &b_nz_cls);
                let list = &mut nz_rest[0];
                list.clear();
                let dirty = &mut dt[k - 2];
                series.class_commutation_with_nonzero_collecting(
                    order, lbuf, lnz, rbuf, rnz, result, dirty, list,
                );
            }
            self.nonzeros = Arc::new(nz);
            self.dirty = Arc::new(dt);
            self.atom_a.clear();
            self.atom_a.extend_from_slice(a_nonzero);
            self.atom_b.clear();
            self.atom_b.extend_from_slice(b_nonzero);
            self.lists_built = true;
            return;
        }

        // Steady state: one fold-parallel dispatch. All levels' jobs are
        // built up front (zeroing every buffer first — the kernel
        // accumulates); the engine plans class-contiguous balanced packs per
        // level and walks them with dynamic pack claims, the level counters
        // ordering cross-level operand reads.
        let mut levels: Vec<Vec<KernelJob<U>>> = self
            .structure
            .levels
            .iter()
            .skip(1)
            .map(|level| {
                level
                    .iter()
                    .filter_map(|&k| match self.structure.nodes[k as usize] {
                        DagNode::Binary { .. } => {
                            // The kernel accumulates into the buffer: it must
                            // start from zero on every fold.
                            let buffer = &mut self.buffers[k as usize - 2];
                            buffer.fill(U::default());
                            Some(KernelJob {
                                a: std::ptr::null(),
                                a_len: 0,
                                a_nonzero: &[],
                                b: std::ptr::null(),
                                b_len: 0,
                                b_nonzero: &[],
                                result: buffer.as_mut_ptr(),
                                result_len: buffer.len(),
                                a_shift: lie_rs::IDENTITY_SHIFTS.as_ptr(),
                                b_shift: lie_rs::IDENTITY_SHIFTS.as_ptr(),
                                r_shift: lie_rs::IDENTITY_SHIFTS.as_ptr(),
                            })
                        }
                        DagNode::Atom(_) => None,
                    })
                    .collect()
            })
            .collect();
        // Fill the operand views after the zeroing pass (the mutable borrows
        // above have ended; operand buffers belong to strictly earlier
        // levels or are the fold's converted inputs, and each result buffer
        // is a distinct allocation, so the raw views never alias within a
        // level).
        for (li, level) in self.structure.levels.iter().enumerate().skip(1) {
            for (jj, &k) in level.iter().enumerate() {
                let (left, right) = match self.structure.nodes[k as usize] {
                    DagNode::Binary { left, right } => (left, right),
                    DagNode::Atom(_) => continue,
                };
                let (a_ptr, a_len, a_nz): (*const U, usize, &[usize]) = match left {
                    0 => (a_cls.as_ptr(), a_cls.len(), &a_nz_cls),
                    1 => (b_cls.as_ptr(), b_cls.len(), &b_nz_cls),
                    id => {
                        let v = &self.buffers[id as usize - 2];
                        (v.as_ptr(), v.len(), &self.nonzeros[id as usize - 2])
                    }
                };
                let (b_ptr, b_len, b_nz): (*const U, usize, &[usize]) = match right {
                    0 => (a_cls.as_ptr(), a_cls.len(), &a_nz_cls),
                    1 => (b_cls.as_ptr(), b_cls.len(), &b_nz_cls),
                    id => {
                        let v = &self.buffers[id as usize - 2];
                        (v.as_ptr(), v.len(), &self.nonzeros[id as usize - 2])
                    }
                };
                // SAFETY: the operand buffers are stable and unmodified for
                // this level's duration (they belong to earlier levels,
                // ordered by the engine's level counters); the result
                // pointers of one level never alias each other.
                let job = &mut levels[li - 1][jj];
                job.a = a_ptr;
                job.a_len = a_len;
                job.a_nonzero = a_nz;
                job.b = b_ptr;
                job.b_len = b_len;
                job.b_nonzero = b_nz;
            }
        }
        series.class_commutation_batch_fold(
            order,
            &mut levels,
            &mut self.gating_cache,
        );
    }

    /// Batch-readiness: the node lists are the built fixed point for these
    /// atom supports, and the accumulator's support already equals the
    /// reachable set — the union of every node's scatter-target list and
    /// the displacement support. Future accumulation can then only change
    /// values at already-nonzero positions (a later cancellation merely
    /// shrinks the support, and the superset lists stay sound), so a whole
    /// run of folds can share one plan and one set of gating masks.
    pub(crate) fn batch_eligible(&self, a_nonzero: &[usize], b_nonzero: &[usize]) -> bool {
        if !self.lists_built || self.atom_a != a_nonzero || self.atom_b != b_nonzero {
            return false;
        }
        let Some(order) = self.class_order.get() else {
            return false;
        };
        let inv = order.inv();
        let words = inv.len();
        let words64 = words.div_ceil(64);
        let mut reach = vec![0u64; words64];
        for list in self.nonzeros.iter() {
            for &i in list {
                reach[i / 64] |= 1 << (i % 64);
            }
        }
        for &w in b_nonzero {
            let i = inv[w] as usize;
            reach[i / 64] |= 1 << (i % 64);
        }
        let mut acc = vec![0u64; words64];
        for &w in a_nonzero {
            let i = inv[w] as usize;
            acc[i / 64] |= 1 << (i % 64);
        }
        acc == reach
    }

    /// Adds the evaluated terms, BCH-weighted, into `target` (a public
    /// basis-order coefficient slice, e.g. the accumulating log
    /// signature's). `rhs` is the displacement atom's public-order slice.
    ///
    /// The accumulation itself runs in the class-contiguous ordering — the
    /// term buffers are class-ordered, so the per-position summation order
    /// matches the public one exactly — and the single epilogue back to
    /// public order runs once, at the end.
    pub(crate) fn accumulate_terms(&mut self, target: &mut [U], rhs: &[U]) {
        let order = self
            .class_order
            .get()
            .expect("accumulate_terms called before evaluate");
        let inv = order.inv();

        // Gather the accumulator's current values and the displacement into
        // class positions. From here to the epilogue every write is
        // sequential in class order.
        let mut acc: Vec<U> = vec![U::default(); inv.len()];
        for (w, c) in target.iter().enumerate() {
            acc[inv[w] as usize] = c.clone();
        }
        let mut rhs_cls: Vec<U> = vec![U::default(); rhs.len().min(inv.len())];
        for (w, c) in rhs.iter().enumerate().take(inv.len()) {
            rhs_cls[inv[w] as usize] = c.clone();
        }

        // Per-position term accumulation in the DAG's term order — the same
        // summation order the public-loop version produced. `raw_mul`/
        // `raw_add_assign` skip the float wrappers' per-op NaN checks
        // (raw-float fast path, see lie-rs `raw_mul`'s NaN policy).
        for (source, weight) in self.structure.terms.iter() {
            let ct: &[U] = match source {
                TermSource::Node(node) => &self.buffers[*node as usize - 2],
                TermSource::Displacement => &rhs_cls,
            };
            for (acc_coeff, comm_coeff) in acc.iter_mut().zip(ct) {
                lie_rs::raw_add_assign(acc_coeff, &lie_rs::raw_mul(comm_coeff, weight));
            }
        }

        // The one epilogue: class-contiguous accumulator back to public
        // basis order (assignment — `acc` already holds the full sum).
        for (k, &src) in inv.iter().enumerate().take(target.len()) {
            target[k] = acc[src as usize].clone();
        }
    }
}
