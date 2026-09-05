# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

Changes on top of the initial `0.1.0` workspace release, from the parallel
log-signature performance campaign.

### Added

- New public API in `lie-rs` around the feasible-decomposition
  structure-constant tables: `LieSeries::with_feasible_decompositions`
  (constructs a series around an existing table instead of rebuilding one),
  `LieSeries::build_feasible_decompositions`, `LieSeries::feasible_table_len`,
  the `FeasibleDecompositions` type re-exported from the crate root, and a new
  `FeasibleDecompositions::len` accessor.
- A 4-lane SIMD "cohort" fold engine for the concatenation kernel: when the
  coefficient type is a representation-transparent raw float and the CPU
  supports AVX2, up to four independent folds share one plan walk and are
  evaluated together in structure-of-arrays lane vectors. Results are
  bit-identical to the scalar engine (plain autovectored array code, no FMA);
  other coefficient types, such as exact rationals, keep the scalar kernel
  bit-for-bit.
- This changelog.

### Changed

- The parallel commutator kernel now divides work by write slot: all words
  and units writing a given kernel output slot form one class, with exact
  gating (each decomposition writes exactly one target slot) replacing the
  previous heuristic letter-content classes. The fold itself is a
  single-phase, per-unit sweep that replaces the former two-phase
  term/fan-in pipeline.
- `LogSignatureBuilder` caches its series plans (Lyndon basis plus
  structure-constant tables) process-wide per configuration, so every build
  after the first skips table construction: steady-state `build()` at
  3x8-class configurations is roughly an order of magnitude cheaper
  (approx. 48 ms to approx. 3.5 ms; approximate). The per-build BCH series
  and commutator DAG template are shared per configuration as well: the
  series is cloned with `Arc`-shared internals and a fresh coefficient
  vector, and the DAG is shallow-cloned, while fold scratch stays per
  instance.
- Leaf groups share their converged steady-state plans through a
  process-wide cache keyed by the feasible-table identity and the exact
  displacement support, so tournament groups with uniform displacements no
  longer re-derive bit-identical plans.
- Commutator DAG node lists are `Arc`-shared copy-on-write: adopting a
  cached steady plan is zero-copy, and a fold never mutates a shared
  snapshot.
- The `SIG_NO_COHORT=1` environment-variable kill switch now swaps the fold
  engine inside the tournament only (cohort vs scalar); the association tree
  is untouched, so results are bit-stable across switch states and thread
  counts. Cohort engagement is capability-gated only (raw-float coefficient
  type plus AVX2, minus the kill switch); the former wide-pool
  thread/work-size gate was removed after measurement showed the SIMD
  engine winning or tying in every tested regime.
- **Breaking:** `LogSignatureBuilder::build` and `build_from_path` now
  require `T: 'static` (and `U: 'static`), because the process-wide plan
  caches key on `TypeId`.

### Performance (approximate)

- End-to-end log-signature builds at the 3x8 configuration on 8 threads are
  roughly 55% faster cumulatively across this campaign versus the
  pre-campaign baseline.
- Light-fold configurations at high thread counts improved sharply: 8x3 at
  32 threads is roughly 29% faster (1.50 ms to 1.09 ms) now that the cohort
  engine engages there.
- Thread scaling at 3x8 improved from 3.9x to 5.8x+ (1 thread to 8
  threads); the per-build table reconstruction that capped scaling is now
  cached.
- The dense-concatenation regressions observed earlier in the campaign are
  eliminated versus the pre-campaign baseline: the final regression floor
  keeps every measured cell within +/-1.4% of its predecessor, and no
  measured regime favors the scalar fallback at width.

### Fixed

- The full-support gating fast path required only the *length* of the exact
  prefix, not the prefix itself, so a batch-recorded node list of matching
  length could fire all-entries-active gating over slots the operand layout
  does not cover. This produced an intermittent exact-rational divergence
  that appeared under CPU oversubscription; both gating prologues now
  require the exact prefix.
- Duplicate entries in a displacement support corrupted batch folds (each
  duplicate re-added the displacement, inflating degree-1 coefficients up to
  6x over the sequential oracle); supports are now deduplicated before
  folding.
- BCH commutator DAG construction strips unreachable ("dead") nodes before
  the level tables are built. Today's construction produces zero dead nodes
  (measured across the supported-degree census), so the strip is a near-free
  invariant guard for future construction changes.
