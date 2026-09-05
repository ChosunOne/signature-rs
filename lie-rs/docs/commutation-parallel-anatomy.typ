#set page(paper: "a4", margin: 1.25cm)
#set text(size: 8.6pt)
#set par(justify: true, leading: 0.52em)
#let sw = rgb("#3d6fa8")   // sweep
#let gt = rgb("#c8873a")   // gating
#let gl = rgb("#8a8f98")   // glue/dispatch
#let c1 = rgb("#3d6fa8")
#let c2 = rgb("#c0504d")
#let c3 = rgb("#6aa84f")
#let c4 = rgb("#b8860b")
#let cols = (c1, c2, c3, c4, c1, c2, c3, c4, c1, c2)
#let cap(sz, body) = text(fill: black.lighten(35%), size: sz, body)

== What is parallelized in the Lie-series commutation kernel

#block(width: 100%, inset: 5pt, fill: rgb("#fff3cd"), radius: 2pt, stroke: 0.5pt + rgb("#c8873a"), text(size: 8pt)[
  *Status: partially superseded (2026-09-05).* The algebra here — feasible pairs, content homogeneity, the anagram-unit partition — is unchanged and still the foundation; see `parallelism-in-the-fold.typ` for the canonical treatment. The execution mechanics changed since: bundles are now cut per unit under a single-phase direct-add sweep (no two-phase term/fan-in), balance by unit pair-weight, and the gating/ticket machinery is cached per exact operand supports. Numbers below are as of 2026-09-01.
])

#grid(columns: (auto, 1fr), column-gutter: 8pt,
  align(left)[*Your model is right for the partition.*],
  align(left)[The pairs to commute are grouped into #emph[anagram units], and units have no dependencies — sweeping them is embarrassingly parallel. What the "2t = ½ runtime" intuition misses is (a) the serial prologue and dispatch #emph[around] the sweep, and (b) the sweep's own parallel #emph[efficiency] — a memory-system property, not a dependency property.],
)

=== The work being divided: feasible pairs, grouped by anagram class

$[A, B] = Sigma_(i < j) (a_i b_j - a_j b_i) [w_i, w_j]$

Every basis word $w_i$ has a #emph[degree] $d_i$ (its length) and a #emph[content] $gamma_i$ (its letter multiset). Brackets are #emph[content-homogeneous]: $[w_i, w_j]$ decomposes only onto words of degree $d_i + d_j$ and content $gamma_i + gamma_j$. The kernel precomputes all #emph[feasible] canonical pairs ($d_i + d_j <= m$) into a flat table and groups them into #emph[anagram units] — all pairs sharing $(t, gamma) = (d_i + d_j, gamma_i + gamma_j)$:

- All of a unit's decompositions land on #emph[that unit's own words] — the degree-$t$ words with content $gamma$. Two different units never write the same result word: #emph[units are disjoint write regions]. No locks; no ordering across units (per-word accumulation order is preserved because a word is only ever written from inside its one unit).
- Entries within a unit are sorted by $p = d_i$ (so $q = t - p$ is forced) → contiguous #emph[p-runs]; degree-support gating keeps only the runs whose $(p, q)$ degrees are backed by nonzeros in $A$ and $B$.

#align(center)[
  #block(width: 92%)[
    #grid(columns: (1fr,) * 10, column-gutter: 1.5pt,
      ..range(10).map(i => rect(width: 100%, height: 15pt, radius: 1pt, fill: cols.at(i), stroke: 0.4pt + black, align(center + horizon, text(size: 5.8pt, fill: white, "u" + str(calc.rem(i, 4) + 1)))))
    )
    #cap(6.4pt)[schematic (4 of the 212 units): degree-$t$ slice of the result vector — one cell per basis word, colour = owning anagram unit (a unit's words are interleaved in memory with its neighbours')]
    #v(3pt)
    #grid(columns: (1fr,) * 10, column-gutter: 1.5pt,
      ..range(10).map(i => rect(width: 100%, height: 12pt, radius: 1pt, fill: cols.at(i).lighten(55%), stroke: 0.4pt + luma(140), align(center + horizon, text(size: 5.6pt, luma(70), "p" + str(calc.rem(i, 4) + 1)))))
    )
    #cap(6.4pt)[the unit's pairs p1…p4, all sharing the same $(t, gamma)$ — they scatter only onto their own colour above]
  ]
]

Real shape (measured, 2×12 grid): 747 basis words, *212 anagram units*, 1 694 feasible pairs; a unit owns ≈3.5 words and ≈8 pairs here. Units are packed into #emph[bundles] of target size `min(max(total_entries/(2·threads), 16), 2048)` — a unit is never split across bundles — and handed to rayon work-stealing, not a fixed assignment. Measured 2t assignment for 2×12: worker A took bundles of {444, 431} entries, worker B {499, 320} — balanced to within 7%.

=== But one kernel call ≠ just the sweep

#align(center)[
  #block(width: 100%)[
    #grid(columns: (8fr, 25fr, 67fr), column-gutter: 2pt, row-gutter: 2pt,
      rect(width: 100%, height: 26pt, fill: gl, radius: 1pt, stroke: 0.4pt + luma(120)),
      rect(width: 100%, height: 26pt, fill: gt.lighten(20%), radius: 1pt, stroke: 0.4pt + luma(120)),
      grid(columns: (1fr,), row-gutter: 2pt,
        rect(width: 100%, height: 12pt, fill: sw.lighten(15%), radius: 1pt, stroke: 0.4pt + luma(120), align(right + horizon, text(size: 5.6pt, fill: white, "worker A ⋯ its bundles ⋯"))),
        rect(width: 100%, height: 12pt, fill: sw, radius: 1pt, stroke: 0.4pt + luma(120), align(right + horizon, text(size: 5.6pt, fill: white, "worker B ⋯ its bundles ⋯"))),
      ),
    )
    #v(3pt)
    #grid(columns: (10pt, auto, 12pt, auto, 12pt, auto), column-gutter: 3pt, align: left,
      rect(width: 7pt, height: 7pt, fill: sw, stroke: 0.3pt + luma(120)), text(size: 6.4pt)[sweep 19.4 µs — the only parallel region (fork · work-stealing · join)],
      rect(width: 7pt, height: 7pt, fill: gt.lighten(20%), stroke: 0.3pt + luma(120)), text(size: 6.4pt)[gating 7.4 µs (walk + presence + masks)],
      rect(width: 7pt, height: 7pt, fill: gl, stroke: 0.3pt + luma(120)), text(size: 6.4pt)[nonzero scan · alloc · bundle build · join ≈ 2.4 µs],
    )
    #cap(6.2pt)[one 3×8 kernel call, 29.1 µs at 1 thread (span-measured); segment widths to scale. The sweep lane shows both workers: each runs its own bundles over disjoint anagram units.]
  ]
]

One call = #emph[prologue] (① nonzero scan → presence bitsets + degree masks, ② gating walk over the whole unit table → active p-run segments) + #emph[dispatch] (③ flatten segments into bundles, fork) + #emph[sweep] (④ per entry: orientation presence tests, term $a_i b_j - a_j b_i$, scatter the stored decomposition) + ⑤ join. The prologue is serial per call but #emph[memoized]: the `GatingCache` keys the active-segment list by the degree-mask pair, so the fold (cache reused across calls) pays ② once ever; and since today a #emph[full-support fast path] (nonzero lists covering `0..degree_start(max_degree)`) skips the bitset scan and ② entirely — presence is all-ones and the active list is precomputed per table.

=== Proportion of the work: anagram-class commutations vs everything else

#align(center)[
  #block(width: 100%)[
    #grid(columns: (92pt, 1fr, 44pt, 68pt), column-gutter: 4pt, row-gutter: 3.5pt,
      align(right + horizon, text(size: 6.4pt, "single 3×8 call")),
      grid(columns: (66fr, 25fr, 9fr), column-gutter: 1pt,
        rect(width: 100%, height: 11pt, fill: sw, stroke: none, align(center + horizon, text(size: 5.8pt, fill: white, "sweep 66%"))),
        rect(width: 100%, height: 11pt, fill: gt, stroke: none, align(center + horizon, text(size: 5.8pt, fill: white, "gating 25%"))),
        rect(width: 100%, height: 11pt, fill: gl, stroke: none, align(center + horizon, text(size: 5.8pt, fill: white, "8%")))),
      text(size: 6.4pt)[29.1 µs], cap(6.4pt)[span-measured],
      align(right + horizon, text(size: 6.4pt, "…after fast path")),
      grid(columns: (92fr, 2fr, 6fr), column-gutter: 1pt,
        rect(width: 100%, height: 11pt, fill: sw, stroke: none, align(center + horizon, text(size: 5.8pt, fill: white, "sweep ≈ 92%"))),
        rect(width: 100%, height: 11pt, fill: gt, stroke: none),
        rect(width: 100%, height: 11pt, fill: gl, stroke: none, align(center + horizon, text(size: 5.8pt, fill: white, "8%")))),
      text(size: 6.4pt)[21.1 µs], cap(6.4pt)[criterion 1t median],
      align(right + horizon, text(size: 6.4pt, "batch 4×10, 64 jobs")),
      grid(columns: (99fr, 0.5fr, 0.5fr), column-gutter: 1pt,
        rect(width: 100%, height: 11pt, fill: sw, stroke: none, align(center + horizon, text(size: 5.8pt, fill: white, "sweep ≈ 99%"))),
        rect(width: 100%, height: 11pt, fill: gt, stroke: none),
        rect(width: 100%, height: 11pt, fill: gl, stroke: none)),
      text(size: 6.4pt)[4.32 ms/job], cap(6.4pt)[batch probe],
    )
  ]
]

The commutations themselves dominate: 66% of a cold single call (25% was gating — fixed today), ≈92% after the fast path, ≈99% in the batch steady state the fold actually runs. The serial remainder is small #emph[proportionally] but its #emph[absolute] size decides small calls: all of 12×2 is 0.87 µs, so ≈2 µs of fork/join alone makes 2t slower (0.24×).

=== So why isn't 2t = ½ the runtime, even at 99% parallel work?

$T_2 approx underbrace(G + delta, "serial: prologue + dispatch") + overbrace(P / E_2, "sweep on 2 threads")$

#text(size: 7pt)[
  #grid(columns: (76pt, 40pt, 48pt, 56pt, 58pt), column-gutter: 4pt, row-gutter: 2.5pt,
    align(left)[*workload*], align(right)[*serial s*], align(right)[*E₂ (sweep)*], align(right)[*2t predicted*], align(right)[*2t measured*],
    align(left)[single 3×8], align(right)[8%], align(right)[≈1.5×], align(right)[1.42×], align(right)[*1.30×*],
    align(left)[batch 4×10], align(right)[1%], align(right)[≈1.8×], align(right)[1.78×], align(right)[*1.74–1.80×*],
  )
]

$E_2 < 2$ is a property of the sweep itself — measured with both workers 100% busy (perf at 9999 Hz; `/proc` tick deltas; ~3.1 GHz on both cores; no idle cores; disjoint memory):
+ *Memory system, not dependencies*: the sweep is a streaming gather–scatter over the pair table and decomposition arrays; at 2t a batch burns ≈15% more cycles for identical work. Bigger bundles don't help (cap 2048 → 2 M: no change), same-CCX pinning doesn't help, page faults ruled out.
+ *Fork/join + bundle build* ≈ 1–2 µs per call — invisible at 4.3 ms/job (4×10 batch), fatal at 0.9 µs/job (12×2).
+ *History*: the bundle builder used to split units across bundles at every boundary — two workers `+=`-ing the same words (lost updates: flaky wrong results #emph[and] worst-case cache-line sharing). Fixed today; the old oracle missed it because its grids fit inside a single bundle.

#v(3pt)
#block(width: 100%, fill: luma(245), inset: 5pt, radius: 2pt, text(size: 6.8pt)[
  *Take-away.* The anagram-unit partition is sound and fully parallel. Where 2t lands is decided by (i) how much of the call is sweep — 66% → 92% → 99% for cold single call → fast path → batch — and (ii) the sweep's memory-bound efficiency $E_2 ≈ 1.5×$ (small grids) to $1.8×$ (4×10). Measured outcomes: single calls 0.24–1.30× at 2t; batch 1.44–1.80× at 2t and up to 5.9× at 8t (`commutator_kernel_batch` bench, baseline `postrace`). Sources: criterion medians (today), tracing spans (pre-fastpath phase split), perf 9999 Hz profile, `/proc` per-thread accounting.
])
