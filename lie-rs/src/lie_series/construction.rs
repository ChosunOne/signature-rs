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

#[cfg(test)]
mod tests {
    use super::*;
    use lyndon_rs::lyndon::{LyndonBasis, Sort};
    use ordered_float::NotNan;

    /// The free Lie algebra is multigraded by letter content: the bracket of
    /// two content-homogeneous elements is content-homogeneous of the summed
    /// multiset, and the Lyndon basis refines the grading. So every
    /// decomposition word of `[basis[i], basis[j]]` must carry exactly the
    /// letter multiset of `basis[i]` + `basis[j]` — the invariant that lets
    /// the commutation kernel parallelize over target anagram classes with
    /// disjoint writes.
    #[test]
    fn decompositions_are_content_homogeneous() {
        for (d, m) in [(2usize, 8usize), (3, 8), (4, 8)] {
            let words = LyndonBasis::<u8>::new(d, Sort::Lexicographical).generate_basis(m);
            let mut letters: Vec<u8> = words
                .iter()
                .filter(|w| w.len() == 1)
                .map(|w| w.letters[0])
                .collect();
            letters.sort_unstable();
            let content = |w: &LyndonWord<u8>| -> Vec<u8> {
                let mut c = vec![0u8; letters.len()];
                for l in &w.letters {
                    let k = letters.iter().position(|a| a == l).unwrap();
                    c[k] += 1;
                }
                c
            };
            let contents: Vec<Vec<u8>> = words.iter().map(content).collect();

            let series = LieSeries::<u8, NotNan<f64>>::new(words, Vec::new());
            let d_len = series.commutator_basis.len();
            for i in 0..d_len {
                for j in (i + 1)..d_len {
                    if let Some((idx, _, _)) = series.decomposition(i, j) {
                        let mut target = contents[i].clone();
                        for k in 0..letters.len() {
                            target[k] += contents[j][k];
                        }
                        for &x in idx {
                            assert_eq!(
                                contents[x as usize], target,
                                "d={d} m={m}: pair ({i}, {j}) decomposition word {x}                                  violates content homogeneity"
                            );
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn decompositions_contain_no_zero_coefficients() {
        for (d, m) in [(2usize, 6usize), (3, 5), (4, 4)] {
            let words = LyndonBasis::<u8>::new(d, Sort::Lexicographical).generate_basis(m);
            let series = LieSeries::<u8, NotNan<f64>>::new(words, Vec::new());
            assert!(
                series
                    .feasible_decompositions
                    .iter_entries()
                    .flat_map(|(_, _, _, coefficients)| coefficients)
                    .all(|c| !c.is_zero()),
                "zero coefficient found in decompositions for d={d}, m={m}"
            );
        }
    }

    /// Independently recomputes the decomposition of every feasible canonical
    /// pair and pins the compact table's contents (indices, coefficients,
    /// ordering, feasibility, and canonicality) against it, plus the mirrored
    /// (negated) query path.
    #[test]
    fn feasible_decompositions_match_independent_recomputation() {
        for (d, m) in [(2usize, 6usize), (3, 5), (4, 4)] {
            let words = LyndonBasis::<u8>::new(d, Sort::Lexicographical).generate_basis(m);
            let series = LieSeries::<u8, NotNan<f64>>::new(words, Vec::new());
            let degree = |x: usize| series.commutator_basis[x].degree();
            let index_of = |term: &CommutatorTerm<NotNan<f64>, u8>| {
                series.commutator_basis_index_map[&term.unit_hash()]
            };
            let basis_set: HashSet<_> = series
                .commutator_basis
                .iter()
                .map(CommutatorTerm::unit_hash)
                .collect();

            let mut feasible = 0;
            for unit in series.feasible_decompositions.units() {
                let (p, _q, target) = (unit.p_mask, unit.target as usize, unit.target as usize);
                assert_eq!(
                    p.iter().map(|w| w.count_ones()).sum::<u32>() as usize >= 1,
                    true,
                    "empty unit"
                );
                assert_eq!(
                    unit.gamma.iter().map(|&x| x as usize).sum::<usize>(),
                    target,
                    "unit target degree mismatch"
                );
                let entries: Vec<_> = series
                    .feasible_decompositions
                    .unit_iter(unit)
                    .map(|(e, _, _)| (e.i as usize, e.j as usize))
                    .collect();
                assert!(
                    entries.windows(2).all(|w| w[0] < w[1]),
                    "unit {:?} entries not sorted by (i, j)",
                    unit.gamma
                );

                for (i, j) in entries {
                    assert!(i < j, "non-canonical pair stored");
                    assert!(
                        p[(degree(i) / 64) as usize] >> (degree(i) % 64) & 1 == 1,
                        "pair ({i}, {j}) in a unit of the wrong degrees"
                    );
                    assert_eq!(
                        degree(i) + degree(j),
                        target,
                        "pair ({i}, {j}) in a unit of the wrong target degree"
                    );
                    assert!(degree(i) + degree(j) <= m, "infeasible pair stored");

                    let mut term = comm![&series.commutator_basis[i], &series.commutator_basis[j]];
                    term.lyndon_sort();
                    let expected: Vec<_> = term
                        .lyndon_basis_decomposition(&basis_set)
                        .into_iter()
                        .filter(|x| !x.coefficient().is_zero())
                        .collect();

                    let expected_indices: Vec<u32> =
                        expected.iter().map(|x| index_of(x) as u32).collect();
                    let expected_coefficients: Vec<_> =
                        expected.iter().map(|x| x.coefficient().clone()).collect();

                    // The canonical query returns the stored data unflagged;
                    // comparing through `decomposition` exercises the
                    // degree-block lookup as well as the contents.
                    let (canonical_indices, canonical_coefficients, swapped) =
                        series.decomposition(i, j).expect("canonical pair query");
                    assert!(!swapped, "(i={i}, j={j}) canonical query flagged");
                    assert_eq!(
                        canonical_indices,
                        &expected_indices[..],
                        "(i={i}, j={j}) indices"
                    );
                    assert_eq!(
                        canonical_coefficients,
                        &expected_coefficients[..],
                        "(i={i}, j={j}) coeffs"
                    );
                    feasible += 1;

                    // The mirrored query returns the same (canonical) data,
                    // flagged so the caller negates it into
                    // [basis[j], basis[i]] orientation.
                    let (mirrored_indices, mirrored_coefficients, swapped) =
                        series.decomposition(j, i).expect("mirrored pair query");
                    assert!(swapped, "(i={i}, j={j}) mirrored query not flagged");
                    assert_eq!(mirrored_indices, &expected_indices[..]);
                    assert_eq!(
                        mirrored_coefficients,
                        &expected_coefficients[..],
                        "(i={i}, j={j}) mirrored data differs from canonical"
                    );
                }
            }
            assert_eq!(series.feasible_decompositions.len(), feasible);
        }
    }

    /// INVARIANT CHECK on the real tables: within one degree-target slice,
    /// two different units must never scatter onto the same basis word (the
    /// batch's packs cut at unit boundaries and rely on per-word
    /// single-writer accumulation).
    #[test]
    fn debug_real_table_unit_word_ownership() {
        for (d, m) in [
            (2usize, 4usize),
            (3usize, 4usize),
            (3usize, 5usize),
            (3usize, 8usize),
        ] {
            let basis =
                lyndon_rs::lyndon::LyndonBasis::<u8>::new(d, lyndon_rs::lyndon::Sort::Lexicographical)
                    .generate_basis(m);
            let len = basis.len();
            let a = LieSeries::<u8, NotNan<f64>>::new(
                basis.clone(),
                (0..len)
                    .map(|k| NotNan::new((k % 17) as f64 / 19.0 - 0.4).unwrap())
                    .collect(),
            );
            let _b = LieSeries::<u8, NotNan<f64>>::new(
                basis,
                (0..len)
                    .map(|k| NotNan::new((k % 13) as f64 / 17.0 - 0.3).unwrap())
                    .collect(),
            );
            let fd = &a.feasible_decompositions;
            let mut owner: std::collections::HashMap<(u8, u32), usize> =
                std::collections::HashMap::new();
            let mut violations = 0usize;
            for (uid, unit) in fd.units().iter().enumerate() {
                let span = fd.entry_span(unit);
                for (entry, next) in span[..span.len() - 1].iter().zip(span[1..].iter()) {
                    let from = entry.decomp_start as usize;
                    let to = next.decomp_start as usize;
                    for &rel in &fd.decomp_indices_rel()[from..to] {
                        let key = (unit.target, rel);
                        match owner.get(&key) {
                            Some(prev) if *prev != uid => {
                                violations += 1;
                                if violations <= 6 {
                                    println!(
                                        "VIOLATION d={d} m={m}: degree {} rel {} by units {} and {}",
                                        unit.target, rel, prev, uid
                                    );
                                }
                            }
                            _ => {
                                owner.insert(key, uid);
                            }
                        }
                    }
                }
            }
            println!(
                "d={d} m={m}: units={} rels={} violations={}",
                fd.units().len(),
                fd.decomp_indices_rel().len(),
                violations
            );
            assert_eq!(violations, 0);
        }
    }

    #[test]
    #[ignore = "probe"]
    fn probe_entry_contribution_ratio() {
        for (d, m) in [
            (2usize, 12usize),
            (3usize, 8usize),
            (4usize, 8usize),
            (2usize, 8usize),
            (3usize, 12usize),
        ] {
            let words: Vec<LyndonWord<u8>> =
                LyndonBasis::<u8>::new(d, Sort::Lexicographical).generate_basis(m);
            let len = words.len();
            let a: LieSeries<u8, NotNan<f64>> =
                LieSeries::new(words, vec![NotNan::new(1.0).unwrap(); len]);
            let table = &a.feasible_decompositions;
            let entries = table.entries();
            let e = entries.len() - 1;
            let c: usize = (0..e)
                .map(|ei| {
                    entries[ei + 1].decomp_start as usize - entries[ei].decomp_start as usize
                })
                .sum();
            let rows_gt1 = (0..e)
                .filter(|&ei| {
                    entries[ei + 1].decomp_start - entries[ei].decomp_start > 1
                })
                .count();
            println!(
                "d={d} m={m}: basis={len} entries={e} contribs={c} rho={:.3} rows_gt1={rows_gt1} ({:.1}%)",
                c as f64 / e as f64,
                100 * rows_gt1 / e.max(1)
            );
        }
    }

    /// Anagram-class scatter locality of the commutation kernel's write set:
    /// per anagram unit, how many distinct 64-byte lines its stores touch and
    /// how far they spread inside the public degree slice, versus the same
    /// stores routed into the class-contiguous (internal) layout.
    #[test]
    #[ignore = "layout stats: run explicitly with --ignored --nocapture"]
    fn dump_scatter_locality_stats() {
        for (d, m) in [(2usize, 12usize), (3, 8), (4, 10), (12, 2)] {
            let basis = lyndon_rs::LyndonBasis::<u8>::new(d, lyndon_rs::Sort::Lexicographical);
            let words = basis.generate_basis(m);
            let len = words.len();
            let a = LieSeries::new(words.clone(), vec![NotNan::new(1.0).unwrap(); len]);
            let table = &a.feasible_decompositions;
            let decomp = table.decomp_indices_rel();
            let entries = table.entries();
            let mut cur_lines = 0usize;
            let mut cls_lines = 0usize;
            let mut total_rmw = 0usize;
            let mut spread_max = 0usize;
            let mut report = String::new();
            for unit in table.units() {
                let (rs, re) = table.degree_range(unit.target);
                let slice = (re - rs) as usize;
                let mut words_seen = std::collections::BTreeSet::new();
                let mut rmw = 0usize;
                for ei in unit.start..unit.end {
                    let from = entries[ei as usize].decomp_start as usize;
                    let to = entries[ei as usize + 1].decomp_start as usize;
                    for &rel in &decomp[from..to] {
                        words_seen.insert(rel as usize);
                    }
                    rmw += to - from;
                }
                let distinct = words_seen.len();
                if distinct == 0 {
                    continue;
                }
                let lo = *words_seen.first().unwrap();
                let hi = *words_seen.last().unwrap();
                let spread = hi - lo + 1;
                let lines_touched = words_seen
                    .iter()
                    .map(|&w| (rs as usize + w) / 8)
                    .collect::<std::collections::HashSet<_>>()
                    .len();
                let lines_cls = (distinct + 7) / 8;
                cur_lines += lines_touched;
                cls_lines += lines_cls;
                total_rmw += rmw;
                spread_max = spread_max.max(spread);
                report.push_str(&format!(
                    "    unit t={:2} words={:5} rmw={:6} slice={:7} spread={:6} ({:5.1}% of slice) lines_cur={:5} lines_cls={:5}\n",
                    unit.target,
                    distinct,
                    rmw,
                    slice,
                    spread,
                    100.0 * spread as f64 / slice.max(1) as f64,
                    lines_touched,
                    lines_cls,
                ));
            }
            println!(
                "{d}x{m}: len={len} store_rmw={total_rmw} lines_cur_sum={cur_lines} lines_cls_sum={cls_lines} ({:.2}x reduction) max_unit_spread={spread_max}",
                cur_lines as f64 / cls_lines.max(1) as f64
            );
            print!("{report}");
        }
    }
}
