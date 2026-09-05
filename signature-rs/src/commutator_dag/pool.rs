//! Cheap plan sharing for the dag pool: `clone_shallow` and the debug
//! rendering of the compiled plan.

use std::fmt;
use std::sync::Arc;

use super::CommutatorDag;

impl<U> CommutatorDag<U> {
    /// Shares the compiled plan; scratch buffers are not copied and are
    /// reallocated on the clone's next [`Self::evaluate`].
    pub(crate) fn clone_shallow(&self) -> Self {
        Self {
            structure: Arc::clone(&self.structure),
            // Scratch buffers stay per-value, but the compiled target
            // lists are inherited by `Arc` reference — a deep-copied
            // accumulator has identical coefficients, hence identical atom
            // supports, hence the same fixed point, and the lists are a
            // support-determined snapshot that the clone's first mutation
            // copies (the original is never mutated through the shared
            // `Arc`). The clone's first fold reuses the lists instead of
            // paying a collecting rebuild — now without the copy.
            buffers: Vec::new(),
            nonzeros: Arc::clone(&self.nonzeros),
            dirty: Arc::clone(&self.dirty),
            atom_a: self.atom_a.clone(),
            atom_b: self.atom_b.clone(),
            lists_built: self.lists_built,
            gating_cache: self.gating_cache.clone(),
            // Scratch: the clone re-allocates its own term buffers on its
            // first plan (see the field doc).
            class_order: self.class_order.clone(),
        }
    }
}

impl<U> fmt::Debug for CommutatorDag<U> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("CommutatorDag")
            .field("nodes", &self.structure.nodes.len())
            .field("terms", &self.structure.terms.len())
            .field("buffers", &self.buffers.len())
            .finish()
    }
}
