use super::*;

#[cfg(feature = "progress")]
use indicatif::{ProgressBar, ProgressStyle};

/// Cheap deterministic per-basis setup shared by every construction path:
/// the commutator-basis terms, their unit-hash membership set, and the
/// unit-hash -> basis-index map.
struct CommutatorPrologue<U, T> {
    commutator_basis: Vec<CommutatorTerm<U, T>>,
    commutator_basis_set: HashSet<u64>,
    commutator_basis_index_map: HashMap<u64, usize>,
}

impl<
    T: Clone + Ord + Generator + Hash + Eq,
    U: Clone + One + Zero + Eq + MulAssign + Neg<Output = U> + Hash + AddAssign + Ord,
> LieSeries<T, U>
{
    #[must_use]
    #[cfg_attr(
        feature = "tracing",
        tracing::instrument(skip_all, fields(basis_len = basis.len()))
    )]
    /// DEBUG/telemetry: per-position Lyndon degree and the degree-slice
    /// starts of the feasible-decomposition table (used by support-size
    /// probes in downstream crates).
    #[doc(hidden)]
    pub fn debug_degree_layout(&self) -> (Vec<u8>, Vec<u32>) {
        self.feasible_decompositions.debug_degree_layout()
    }

    /// Cheap deterministic per-basis setup shared by every construction
    /// path: the commutator-basis terms, their unit-hash membership set,
    /// and the unit-hash -> basis-index map.
    fn commutator_prologue(basis: &[LyndonWord<T>]) -> CommutatorPrologue<U, T> {
        let mut commutator_basis = Vec::<CommutatorTerm<U, T>>::with_capacity(basis.len());
        // Only the tracing build reports it; kept for the debug! shorthand.
        #[cfg_attr(not(feature = "tracing"), allow(unused_variables))]
        let max_degree = if basis.is_empty() {
            0
        } else {
            basis[basis.len() - 1].len()
        };
        for word in basis {
            commutator_basis.push(CommutatorTerm::from(word));
        }
        #[cfg(feature = "tracing")]
        tracing::debug!(
            max_degree,
            basis_len = commutator_basis.len(),
            "converted Lyndon basis to commutator terms"
        );
        let commutator_basis_set = commutator_basis
            .iter()
            .map(commutator_rs::CommutatorTerm::unit_hash)
            .collect::<HashSet<_>>();

        let mut commutator_basis_index_map = HashMap::new();
        for (i, a) in commutator_basis.iter().enumerate() {
            commutator_basis_index_map.insert(a.unit_hash(), i);
        }
        #[cfg(feature = "tracing")]
        tracing::debug!("built commutator basis maps");
        CommutatorPrologue {
            commutator_basis,
            commutator_basis_set,
            commutator_basis_index_map,
        }
    }

    /// Builds the feasible-decomposition (structure-constant) table for a
    /// basis: for every feasible pair `(i, j)` the commutator `[w_i, w_j]`
    /// is rewritten into the Lyndon basis. Deterministic — pairs are
    /// visited in index order and each decomposition's output terms are
    /// sorted — so the table is a pure function of the basis (plus the
    /// coefficient type's arithmetic). This is the expensive part of
    /// [`LieSeries::new`]; it is factored out so callers can build the
    /// table once and share it via
    /// [`LieSeries::with_feasible_decompositions`].
    #[must_use]
    pub fn build_feasible_decompositions(basis: &[LyndonWord<T>]) -> FeasibleDecompositions<U> {
        let prologue = Self::commutator_prologue(basis);
        Self::structure_constants(
            basis,
            &prologue.commutator_basis,
            &prologue.commutator_basis_set,
            &prologue.commutator_basis_index_map,
        )
    }

    /// The structure-constant pass over the pair table. `basis` must be the
    /// basis `commutator_basis` / `commutator_basis_set` /
    /// `commutator_basis_index_map` were derived from.
    fn structure_constants(
        basis: &[LyndonWord<T>],
        commutator_basis: &[CommutatorTerm<U, T>],
        commutator_basis_set: &HashSet<u64>,
        commutator_basis_index_map: &HashMap<u64, usize>,
    ) -> FeasibleDecompositions<U> {
        let max_degree = if basis.is_empty() {
            0
        } else {
            basis[basis.len() - 1].len()
        };
        // Per-word letter contents: the table is grouped by (target degree,
        // content class), which requires each word's letter multiset.
        let mut alphabet: Vec<T> = basis
            .iter()
            .filter(|w| w.len() == 1)
            .map(|w| w.letters[0].clone())
            .collect();
        alphabet.sort();
        alphabet.dedup();
        let basis_contents: Vec<Vec<u8>> = basis
            .iter()
            .map(|w| {
                let mut c = vec![0u8; alphabet.len()];
                for l in &w.letters {
                    let k = alphabet
                        .binary_search(l)
                        .expect("basis letter missing from alphabet");
                    c[k] += 1;
                }
                c
            })
            .collect();
        let mut feasible_builder = feasible_decompositions::Builder::new(&basis_contents);

        #[cfg(feature = "tracing")]
        let _structure_constants_span = tracing::debug_span!(
            "compute_structure_constants",
            pairs = basis.len() * basis.len(),
            max_degree
        )
        .entered();
        #[cfg(feature = "progress")]
        let pb = {
            let style = ProgressStyle::with_template(
                "[{elapsed_precise}] [{bar:35.cyan/blue}] {pos}/{len} structure-constant rows ({msg})",
            )
            .unwrap()
            .progress_chars("=>-");
            let pb = ProgressBar::new(basis.len() as u64).with_style(style);
            pb
        };
        for (i, a) in commutator_basis.iter().enumerate() {
            #[cfg(feature = "progress")]
            {
                pb.set_message(format!("degree {}", a.degree()));
                pb.inc(1);
            }
            #[cfg(feature = "tracing")]
            if i == 0 || commutator_basis[i - 1].degree() != a.degree() {
                tracing::debug!(
                    row = i,
                    degree = a.degree(),
                    "computing structure constants"
                );
            }
            for (j, b) in commutator_basis.iter().enumerate() {
                if j <= i || max_degree < a.degree() + b.degree() {
                    continue;
                }
                #[cfg(feature = "tracing")]
                tracing::trace!(i, j, "computing commutator term");
                let mut term = comm![a, b];
                #[cfg(feature = "tracing")]
                tracing::trace!(i, j, "sorting commutator term");
                term.lyndon_sort();

                // For some non-basis term T, and its decomposition to basis terms A, B, and C, ...
                // T -> [c_1 * A, c_2 * B, c_3 * C, ...]
                #[cfg(feature = "tracing")]
                tracing::trace!(i, j, "decomposing into Lyndon basis");
                let basis_terms = term
                    .lyndon_basis_decomposition(&commutator_basis_set)
                    .into_iter()
                    .filter(|x| !x.coefficient().is_zero())
                    .collect::<Vec<_>>();
                // Get the coefficients [c_1, c_2, c_3, ...]
                let basis_term_coefficients = basis_terms
                    .iter()
                    .map(|x| x.coefficient().clone())
                    .collect::<Vec<_>>();
                // Get the indices i of A, B, C, ...
                // [i_A, i_B, i_C, ...]
                let basis_term_indices = basis_terms
                    .into_iter()
                    .map(|x| commutator_basis_index_map[&x.unit_hash()])
                    .collect::<Vec<_>>();
                #[cfg(feature = "tracing")]
                tracing::trace!(
                    i,
                    j,
                    terms = basis_term_indices.len(),
                    "decomposed commutator into basis terms"
                );
                feasible_builder.push(i, j, &basis_term_indices, &basis_term_coefficients);
            }
        }

        #[cfg(feature = "tracing")]
        tracing::debug!(
            basis_len = basis.len(),
            "finished computing structure constants"
        );
        #[cfg(feature = "progress")]
        pb.finish_with_message("done");
        feasible_builder.finish()
    }

    /// Builds an empty series over `basis` with a pre-built
    /// feasible-decomposition table, skipping the structure-constant
    /// computation. The table must have been produced for exactly this
    /// basis (see [`LieSeries::build_feasible_decompositions`]); sharing
    /// one table across series over the same basis is sound because the
    /// table is read-only after construction.
    pub fn with_feasible_decompositions(
        basis: Arc<Vec<LyndonWord<T>>>,
        coefficients: Vec<U>,
        feasible_decompositions: Arc<FeasibleDecompositions<U>>,
    ) -> Self {
        let prologue = Self::commutator_prologue(&basis);
        let max_degree = if basis.is_empty() {
            0
        } else {
            basis[basis.len() - 1].len()
        };
        Self {
            basis,
            commutator_basis: Arc::new(prologue.commutator_basis),
            commutator_basis_index_map: Arc::new(prologue.commutator_basis_index_map),
            coefficients,
            feasible_decompositions,
            max_degree,
        }
    }

    /// Identity of the shared structure-constant table (diagnostics and
    /// tests: series over the same cached plan report the same id).
    #[doc(hidden)]
    #[must_use]
    pub fn table_id(&self) -> usize {
        Arc::as_ptr(&self.feasible_decompositions) as usize
    }

    pub fn new(basis: Vec<LyndonWord<T>>, coefficients: Vec<U>) -> Self {
        let prologue = Self::commutator_prologue(&basis);
        let feasible_decompositions = Self::structure_constants(
            &basis,
            &prologue.commutator_basis,
            &prologue.commutator_basis_set,
            &prologue.commutator_basis_index_map,
        );
        let max_degree = if basis.is_empty() {
            0
        } else {
            basis[basis.len() - 1].len()
        };
        Self {
            basis: Arc::new(basis),
            commutator_basis: Arc::new(prologue.commutator_basis),
            commutator_basis_index_map: Arc::new(prologue.commutator_basis_index_map),
            coefficients,
            feasible_decompositions: Arc::new(feasible_decompositions),
            max_degree,
        }
    }
}

impl<T, U> LieSeries<T, U> {
    #[must_use]
    pub fn with_coefficients(&self, coefficients: Vec<U>) -> Self {
        Self {
            basis: Arc::clone(&self.basis),
            commutator_basis: Arc::clone(&self.commutator_basis),
            feasible_decompositions: Arc::clone(&self.feasible_decompositions),
            commutator_basis_index_map: Arc::clone(&self.commutator_basis_index_map),
            coefficients,
            max_degree: self.max_degree,
        }
    }

    /// Decomposition of the bracket `[basis[i], basis[j]]` onto the basis,
    /// for pairs stored in the feasible table (`i != j` and degree-feasible).
    /// Returns parallel slices (basis indices, non-zero coefficients).
    #[must_use]
    pub fn decomposition(&self, i: usize, j: usize) -> Option<(&[u32], &[U], bool)> {
        self.feasible_decompositions.get(i, j)
    }
}
