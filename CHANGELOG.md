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

### Performance (measured against the 0.1.1 release)

Methodology: the current criterion harness was ported onto unmodified 0.1.1
library sources, so both sides run an identical benchmark driver — only the
library under test differs. Cells run interleaved x2 with 10 samples per
binary; end-to-end (e2e) = a 1000-point path build. Medians:

| e2e cell | 0.1.1 release | this version | delta |
|---|---|---|---|
| 2x8, 1 thread | 197.9 ms | 10.6 ms | -94.7% |
| 2x8, 8 threads | 199.8 ms | 3.5 ms | -98.3% |
| 2x8, 32 threads | 199.9 ms | 4.2 ms | -97.9% |
| 8x3, 8 threads | 11.8 ms | 1.0 ms | -91% |
| 8x3, 32 threads | ~31.7-49.8 ms (noisy) | 12.1 ms | ~-62% |

Notes:

- The 0.1.1 build does not complete the 3x8 e2e cell at any thread count:
  it hangs during warmup (no output in a 25-minute window). That cell is
  reported as unmeasurable rather than inferred. For reference, this
  version's absolute numbers there are 1t ~262 ms, 8t ~47 ms, 32t ~51 ms.
- The 0.1.1 build shows ~198 ms at 2x8 regardless of thread count (1t, 8t,
  32t) — pre-campaign code has no fold parallelism; scaling in this version
  at 3x8-class workloads reaches ~5.8x at 8 threads (1t->8t, measured within
  the campaign as 3.93x -> 5.83x across the plan-cache work).
- Engine-choice A/Bs inside this version (cohort SIMD engine vs scalar
  fallback, not a version comparison): the cohort engine wins or ties every
  measured regime at 1-32 threads; the previous wide-pool engagement gate
  was removed on that evidence (8x3 32t 1.54 -> 1.09 ms on the
  counterfactual measurement).

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
