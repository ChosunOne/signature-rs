use lie_rs::LieSeriesGenerator;
use lie_rs::{
    ClassOrderedCommutation, ClassOrder, FeasibleDecompositions, LieSeries,
    bch_series_generator::BchSeriesGenerator,
};
use lyndon_rs::lyndon::LyndonWord;
use lyndon_rs::{
    generators::Generator,
    lyndon::{LyndonBasis, Sort},
};
use ndarray::{ArrayView, Axis, Dimension, RemoveAxis};

use num_traits::{FromPrimitive, One, Zero};
use std::any::{Any, TypeId};
use std::collections::HashMap;
use std::ops::{Mul, Sub};
use std::sync::{Arc, Mutex, OnceLock};
use std::{
    fmt::Debug,
    hash::Hash,
    ops::{AddAssign, Div, Index, IndexMut, MulAssign, Neg, SubAssign},
};

use crate::commutator_dag::{
    CohortLane, CommutatorDag, SharedDirtyLists, SharedNodeLists, cohort_capable,
    cohort_type_capable,
};

/// Process-wide cache of per-configuration series plans: the Lyndon basis
/// and the structure-constant (`FeasibleDecompositions`) table that
/// [`LogSignatureBuilder::build`] used to reconstruct on every call.
///
/// The plan is a pure function of the generator type `T` (letter values),
/// the coefficient type `U` (table arithmetic), the alphabet size, the
/// basis `Sort`, and the maximum degree — so entries are keyed by exactly
/// that tuple and shared across every builder (and every builder copy:
/// the builder is `Copy`) with the same configuration. Tables are
/// read-only after construction ([`LieSeries::with_feasible_decompositions`]
/// injects them unchanged), and their interior lazy fields
/// (`class_order`, full-support gating) are `OnceLock`-guarded pure
/// computations, so sharing them across threads and folds is sound.
///
/// A cache miss builds the plan **outside** the map lock; the
/// insert-if-absent then hands back whichever identical plan won the race,
/// so concurrent first builds never block each other and never observe a
/// different table.
mod series_plan_cache;

/// Process-wide cache of tournament LEAF steady plans: the fixed-point
/// reachable-support plan (`a_plan`) plus the plan dag's collected node
/// scatter lists, keyed by the displacement support they were derived
/// from.
///
/// Every field is a pure function of (feasible table, class order,
/// supports): collection is support-determined — the kernel's scatter set
/// is fixed by the supports' gating alone, the collect sweep's values are
/// throwaway — so a snapshot adopted by a later group is bit-identical to
/// what the group's own fixed-point derivation would produce, minus the
/// collecting kernel sweeps (two per derivation; with ~8-32 leaf groups
/// per tournament batch this was measured at 4.5% of e2e samples at 1
/// thread and 13.2% at 4 threads — it grows with pool width because the
/// sweep divides by threads while per-group planning does not).
///
/// The key's table id pins the whole basis configuration (the table was
/// built from exactly that basis: alphabet, max degree, sort, letter
/// type), and the id is a monotonic counter — never reused, unlike an
/// address. Snapshot contents are letter- and coefficient-type-free
/// (`usize`/`u64`/class order), so the map is concretely typed.
mod leaf_steady_cache;

mod engine;

mod builder;
mod signature;

#[cfg(test)]
mod test;

pub use builder::LogSignatureBuilder;
pub use signature::LogSignature;
