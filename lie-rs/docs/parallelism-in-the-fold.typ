#set page(paper: "a4", margin: 1.1cm)
#set text(size: 9.8pt)
#set par(justify: true, leading: 0.58em, spacing: 0.75em)
#show heading.where(level: 1): set text(size: 13pt)
#show heading.where(level: 2): set block(above: 10pt, below: 5pt)
#show heading.where(level: 3): set block(above: 7pt, below: 3.5pt)
#let cap(sz, body) = text(fill: black.lighten(35%), size: sz, body)
#let sw = rgb("#3d6fa8")    // sweep / commutation work
#let gt = rgb("#c8873a")    // gating / planning
#let gv = rgb("#6aa84f")    // gather
#let ac = rgb("#8e6bb0")    // accumulate
#let rs = rgb("#c0504d")    // serial / synchronization / caution
#let idl = luma(228)        // idle / dead
#let boxc(w, h, fillc, label, txtc) = rect(width: w, height: h, fill: fillc, radius: 1.5pt, stroke: 0.4pt + luma(150), align(center + horizon, text(size: 5.9pt, fill: txtc, label)))

#align(center)[
  #text(size: 14.5pt, weight: "bold")[Parallelism in the commutation fold]
  #v(2pt)
  #text(size: 10pt)[From the algebra of series commutation to the execution schedule of a log-signature computation]
  #v(3pt)
  #cap(7.5pt)[signature-rs / lie-rs — 2026-09-02. Audience: readers comfortable with Lyndon bases, Lie series commutation, and log signatures; no familiarity with this codebase assumed. Concepts and methods only.]
]

== The three scales of the computation

A log signature of a path is built by folding the displacements of consecutive path points into a running Lie series, one Baker–Campbell–Hausdorff (BCH) concatenation at a time. The whole computation nests three scales, and the honest first answer to "where can threads go?" is different at each scale:

#grid(columns: (70pt, 1fr, 105pt), column-gutter: 7pt, row-gutter: 4pt,
  align(horizon, [*The path* \ (scale 1)]),
  [An associative reduction of the per-segment displacements by BCH concatenation: $S = (((d_1 ⊕ d_2) ⊕ d_3) ⊕ dots.c)$. Computed today as a *strictly sequential chain* — fold $k+1$'s input is fold $k$'s output.],
  align(horizon, cap(7pt)[_serial between folds; associativity is the only escape hatch (§8)_]),
  align(horizon, [*The fold* \ (scale 2)]),
  [One concatenation $A ⊕ B$ evaluates a fixed *directed acyclic graph* of Lie bracket nodes, compiled once from the BCH series truncated at the working degree. Nodes on one dependency level are independent; levels are dependency-ordered.],
  align(horizon, cap(7pt)[_parallel within each level; barriers between levels (§5)_]),
  align(horizon, [*The commutation* \ (scale 3)]),
  [Each node evaluates one series commutation $[X, Y]$ over the Lyndon basis: thousands of independent pair contributions, partitioned into *anagram units* that never share a target word.],
  align(horizon, cap(7pt)[_embarrassingly parallel inside a unit-partition (§2–4)_]),
)

#align(center)[
  #block(width: 98%, inset: 4pt)[
    #stack(dir: ltr, spacing: 3pt,
      cap(6.5pt)[path~~], boxc(56pt, 13pt, rs, "d₁ ⋈ d₂ ⋈ …", white),
      sym.arrow.r, boxc(62pt, 13pt, gv.darken(10%), "fold k: BCH DAG", white),
      sym.arrow.r,
      boxc(42pt, 13pt, sw, "level ℓ", white),
      boxc(44pt, 13pt, sw.lighten(20%), "level ℓ+1", white),
      boxc(44pt, 13pt, sw.lighten(40%), "level ℓ+2", black),
    )
    #v(4pt)
    #stack(dir: ltr, spacing: 1.5pt,
      cap(6.5pt)[units~~],
      boxc(30pt, 9pt, sw.lighten(35%), "unit", black), boxc(30pt, 9pt, sw.lighten(50%), "unit", black),
      boxc(30pt, 9pt, sw.lighten(35%), "unit", black), boxc(30pt, 9pt, sw.lighten(50%), "unit", black),
      boxc(30pt, 9pt, sw.lighten(35%), "unit", black), boxc(30pt, 9pt, sw.lighten(50%), "unit", black),
      cap(6pt)[~~… hundreds of anagram units per node])
    #v(3pt)
    #cap(6.5pt)[Nesting: a serial chain of folds; each fold a level-ordered DAG of commutations; each commutation a bag of anagram units. Parallelism exists at the two inner scales; the outermost scale is a data dependency, not a scheduling choice.]
  ]
]

== The algebra: what one commutation computes

=== Brackets over the Lyndon basis

A Lie series $A = sum_u a_u w_u$ is a formal sum over the Lyndon words — the canonical basis of the free Lie algebra on $d$ letters. The commutator of two series distributes over pairs of basis words:

$ [A, B] = sum_(u < v) (a_u b_v - a_v b_u) [u, v] $

The sum runs over *canonical pairs* $u < v$ of Lyndon words (antisymmetry halves the work immediately). Each bracket $[u, v] = u v - v u$ re-expands in the Lyndon basis; the expansion has two rigid structural symmetries, and they are the entire foundation of the parallel layout:

#text(size: 8.4pt)[
  #grid(columns: (110pt, 1fr), column-gutter: 6pt, row-gutter: 3pt,
    align(left)[*Degree additivity*], align(left)[$|u| + |v| = t$: every word in the expansion of $[u,v]$ has degree $t$. Truncating the series at degree $m$ discards all pairs with $|u| + |v| > m$ up front.],
    align(left)[*Content homogeneity*], align(left)[$[u,v]$ expands only over words with the same letter multiset as $u v$. Writing $gamma(w)$ for a word's content vector in $NN^d$: a pair with $gamma(u) + gamma(v) = gamma$ can only ever touch basis words of content $gamma$.],
  )
]

A *feasible decomposition* of a basis word $w$ is a canonical pair $(u, v)$ such that $w$ appears (with nonzero integer coefficient) in the expansion of $[u,v]$, together with that coefficient $c^w_(u,v) in QQ$. The coefficient of $w$ in $[A,B]$ is then one scalar fan-in over all its feasible decompositions:

#align(center)[
  #block(width: 92%, inset: 6pt, fill: luma(250), radius: 2pt)[
    $ [A,B]_w = sum_((u,v) arrow w) c^w_(u,v) dot underbrace((a_u b_v - a_v b_u), "term — 1–2 multiplies") $
    #v(2pt)
    #cap(6.6pt)[Each feasible pair contributes its bracket-coefficient times the antisymmetric product of the operand values. All $c^w_(u,v)$ are exact small rationals — precomputable once per basis, before any data arrives.]
  ]
]

#align(center)[
  #block(width: 96%, inset: 4pt)[
    #stack(dir: ltr, spacing: 2.5pt,
      boxc(46pt, 12pt, sw, "w", white),
      cap(6.5pt)[~~←~], boxc(46pt, 12pt, sw.lighten(25%), "(u,v)₁", white),
      boxc(14pt, 12pt, white, "·", luma(80)), boxc(46pt, 12pt, sw.lighten(25%), "(u,v)₂", white),
      boxc(14pt, 12pt, white, "·", luma(80)), boxc(46pt, 12pt, sw.lighten(25%), "(u,v)₃", white),
      cap(6.5pt)[~~(few: usually 1–3 pairs feed one word)],
    )
    #v(3pt)
    #stack(dir: ltr, spacing: 2.5pt,
      boxc(46pt, 12pt, white, " ", luma(60)),
      cap(6.5pt)[~~~~], boxc(46pt, 12pt, white, "c₁", luma(60)),
      boxc(14pt, 12pt, white, " ", white), boxc(46pt, 12pt, white, "c₂", luma(60)),
      boxc(14pt, 12pt, white, " ", white), boxc(46pt, 12pt, white, "c₃", luma(60)),
    )
    #cap(6.5pt)[One word $w$'s coefficient = a short inner product: its (pair, rational) decomposition row dotted with the antisymmetric pair-terms. Both the pair table and the rational rows are *structure constants*: they depend only on the basis, never on the data.]
  ]
]

=== The precomputed structure-constant table

Because every ingredient except the $a_u$ and $b_v$ is data-independent, the entire commutation is driven by a single precomputed table:

#text(size: 8.4pt)[
  #grid(columns: (96pt, 1fr), column-gutter: 6pt, row-gutter: 3pt,
    align(left)[*Entries*], align(left)[All feasible canonical pairs $(i, j)$, packed as tiny fixed-size records, *degree-blocked*: sorted by (target degree $t$, then $p = |u|$, $q = t - p$). Degree blocking is what lets a truncated fold skip whole $(p,q)$ regions of the table without testing them.],
    align(left)[*Decomposition rows*], align(left)[Per entry, the run of (target word, coefficient) pairs it feeds — the words whose expansion contains $[u,v]$. Most rows have length 1–3: almost every pair feeds its own distinguished word (the word $u v$ itself when it is Lyndon) and almost nothing else.],
    align(left)[*Anagram-unit bounds*], align(left)[Markers ($§3$) delimiting the groups of entries that can write to the *same* set of target words — the unit is the conflict boundary for parallel writes.],
  )
]

== Anagram classes and units: the partition that makes parallelism safe

=== Grouping by content

Group all degree-$t$ basis words by their content vector $gamma$. Each group — an *anagram class* — is the set of words obtainable from one another by permuting letters. Now consider a canonical pair $(u,v)$ with $gamma(u) + gamma(v) = gamma$, $|u|+|v| = t$: by content homogeneity, *everything it produces lands on the single class $(t, gamma)$*.

The kernel groups entries by this pair — a $(t, gamma)$-group is an *anagram unit* — and the payoff is the key structural fact of the whole engine:

#align(center)[
  #block(width: 84%, inset: 6pt, fill: rgb("#eef4ea"), radius: 2pt)[
    #align(center)[*Two distinct units never write to the same result word.*]
    #v(2pt)
    #align(left, text(size: 8.2pt)[Each result word belongs to exactly one anagram class, and each unit feeds exactly one class. Therefore a whole commutation decomposes into units that are *disjoint write regions*: any two of them can be swept concurrently with no locks, no atomics, and no order constraints between them.])
  ]
]

#align(center)[
  #block(width: 98%)[
    #grid(columns: (52pt,) + (1fr,) * 12, column-gutter: 1.5pt, row-gutter: 2pt,
      align(left + horizon, cap(6.5pt)[result slice \ (degree $t$)]),
      ..(("w",)*12).enumerate().map(((i, _)) => rect(width: 100%, height: 14pt, radius: 1pt, fill: (sw, sw.darken(12%), gv, gv.darken(12%), sw, ac, gv, sw.darken(20%), sw, gv.darken(20%), ac.lighten(30%), ac.lighten(10%)).at(i).lighten(20%), stroke: 0.4pt + luma(120), align(center + horizon, text(size: 5.4pt, fill: white, "w" + str(i+1)))))
    )
    #cap(6.4pt)[A degree-$t$ slice of the result vector, one cell per basis word, coloured by its anagram class: four classes (blue, green, purple and their shades). Words of one class are *contiguous in memory* (the basis is permuted so classes are consecutive — the class-contiguous ordering, §5).]
    #v(4pt)
    #grid(columns: (52pt,) + (1fr,) * 12, column-gutter: 1.5pt,
      align(left + horizon, cap(6.5pt)[units \ (pair groups)]),
      ..((sw, sw, sw, idl, gv, gv, idl, idl, sw, sw, ac, ac)).enumerate().map(((i, c)) => rect(width: 100%, height: 12pt, radius: 1pt, fill: c.lighten(if c == idl { 0% } else { 45% }), stroke: 0.4pt + luma(140), align(center + horizon, text(size: 5.4pt, fill: if c == idl { luma(120) } else { luma(60) }, if c == idl { "·" } else { "pairs" }))))
    )
    #cap(6.4pt)[The corresponding units in the entry table: each unit's canonical pairs (one sub-row per unit shown for the blue, green and purple classes; grey cells are units gated off because the current operand supports have no words of those degrees). Every unit scatters only onto words of its own colour above — never anywhere else.]
  ]
]

Measured shape of this partition (2 letters, degree ≤ 12): 747 basis words, 212 units, ≈1700 feasible pairs in one commutation; ≈3.5 words and ≈8 pairs per unit. At high degree the classes grow quickly (the class at content $(6,6)$ holds dozens of words), so units are numerous enough to fill many threads — but sized wildly differently, which matters for scheduling ($§5$).

=== Why a unit may never be split

Inside a unit, several canonical pairs feed the *same* result word (the word's decomposition row, §2). The result word accumulates its contributions sequentially — addition of floating-point values is not associative, so *the mathematical definition includes a fixed per-word accumulation order* (the order of the decomposition row). Two workers writing `+=` into the same word would race (lost updates) and would break reproducibility even when they don't (interleaving changes the rounding). Hence:

#align(center)[
  #block(width: 88%, inset: 5pt, fill: rgb("#f7ecec"), radius: 2pt)[
    *The anagram unit is the atomic scheduling unit.* Packs of parallel work are cut *only* at unit boundaries — a whole unit's entries always stay together, in table order, on one worker. Within a unit the sequential row order is honored by construction; across units no ordering exists at all. This is what makes results *bit-identical to the serial sweep at any thread count* — a property the codebase maintains as an invariant (and guards with cross-thread-count identity tests). The alternative — splitting units and compensating with atomics or privatized partial sums — costs either synchronization on the hot path or extra passes over memory, and an early version of the packer that accidentally cut inside units produced exactly the flaky wrong-results race this boundary exists to prevent.
  ]
]

== One commutation, parallelized end to end

A single evaluation of $[X,Y]$ has four phases; two are parallel, two are serial-but-amortized. The architecture's recurring theme is: *anything whose value depends only on the operand supports (which words are nonzero), not on their values, can be precomputed and cached — supports change rarely, values change every fold.*

#text(size: 8.2pt)[
  #grid(columns: (66pt, 1fr, 72pt), column-gutter: 6pt, row-gutter: 3.5pt,
    align(left)[*① Gating*], align(left)[From the two operands' supports, build presence bit-masks, then walk the entry table once to find *active segments*: the degree-runs whose $(p,q)$ could produce anything. Output: per-unit active entry lists.], align(horizon, cap(6.6pt)[serial per call, *cached* by support fingerprint]),
    align(left)[*② Tickets*], align(left)[For each active segment, resolve *once* the per-entry presence tests (is $u$ present in $X$? $v$ in $Y$? both reversed?) and pack the survivors as compact *tickets*: the entry's table index plus two orientation bits.], align(horizon, cap(6.6pt)[part of the cached gating — not recomputed per fold]),
    align(left)[*③ Bundle*], align(left)[Cut the active units into a handful of *packs* of roughly equal total work, cutting only at unit boundaries.], align(horizon, cap(6.6pt)[cheap planning, parallel-safe by §3]),
    align(left)[*④ Sweep*], align(left)[Workers claim packs dynamically. Per ticket: compute the antisymmetric term (1–2 multiplies) and scatter-multiply its decomposition row into the result.], align(horizon, cap(6.6pt)[the *only* heavy parallel region]),
  )
]

#align(center)[
  #block(width: 100%)[
    #grid(columns: (30pt, 46fr, 20fr, 34fr), column-gutter: 2pt, align: horizon,
      [],
      boxc(100%, 15pt, gt, "gating + tickets (cached)", white),
      boxc(100%, 15pt, luma(150), "bundle", white),
      stack(dir: ttb, spacing: 2pt,
        boxc(100%, 11.5pt, sw, "worker A: pack 1", white),
        boxc(100%, 11.5pt, sw.darken(15%), "worker B: pack 2", white),
      ),
    )
    #v(2.5pt)
    #grid(columns: (1fr,) * 12, column-gutter: 1.2pt,
      ..range(12).map(i => {
        let dead = (1, 3, 8).contains(i)
        rect(width: 100%, height: 10pt, radius: 1pt, fill: if dead { idl } else { sw.lighten(30%) }, stroke: 0.4pt + luma(140), align(center + horizon, text(size: 5.2pt, fill: if dead { luma(110) } else { white }, if dead { "✗" } else { "t" }))) })
    )
    #cap(6.3pt)[Ticket compaction in one segment: surviving entries (blue) become a dense ticket stream; entries that fail presence (grey, ✗) vanish from the sweep entirely. Measured at 2×12: 31.7% of all entry visits were dead — tickets eliminate both the dead visits and the repeated bit-tests from every fold of the batch, for the price of one pass per support-change.]
  ]
]

Two numbers characterize why these phases are arranged this way:

- *The sweep visits hundreds of thousands of entries per call* (measured ≈257,000 tickets per fold at 2×12); the gating/ticket phases are O(table) once per support change. Across a batch of 1000 folds, the per-fold marginal cost of ① and ② is zero.
- *The serial fraction around a single call* (fork/join, task flattening, bundle build) is a fixed few microseconds — negligible against a millisecond sweep at 2×12, but *larger than the entire kernel* at 12×2 (66 entries). This asymmetry returns in §7 as the work-adaptive scheduler.

=== The sweep itself

Per ticket the arithmetic is minimal: one or two multiplies for the term, then one fused multiply-add per decomposition-row element. What surrounds it is memory: the entry record, the operand words it indexes, the decomposed coefficients, and the scattered result words. Measured on real hardware, the per-entry time is *uniform* (~11 ns/entry at 1 thread, ≈44 cycles) regardless of the pair's degree or row length — which is why balancing packs by entry count balances time (a fancier work metric was implemented, measured slower, and reverted) — and why the sweep parallelizes cleanly until *cache-coherence and streaming bandwidth*, not dependencies, set the limit (per-entry cost rises ≈13% merely because two cores stream the same tables; ≈9× at 32 threads in the worst place — the accumulate traffic of §6).

== The fold: a DAG of commutations with level barriers

=== From the BCH series to a DAG

Truncated at degree $m$, the BCH series $A ⊕ B = A + B + ½[A,B] + 1/12[A,[A,B]] - 1/24[B,[A,B]] + dots.c$ has a fixed bracket structure — it depends only on $m$, not on data. The fold compiler turns this bracket forest into a DAG:

#text(size: 8.4pt)[
  #grid(columns: (96pt, 1fr), column-gutter: 6pt, row-gutter: 3pt,
    align(left)[*Nodes*], align(left)[Bracket subterms, *interned* so shared subtrees (the same nested bracket can appear in many BCH terms) are computed once. Each node is one commutation kernel call — and, by content homogeneity, each node writes *exactly one anagram class* of target words.],
    align(left)[*Atoms*], align(left)[The two operands: the running accumulator ($A$) and the folded-in displacement ($B$).],
    align(left)[*Levels*], align(left)[A node's level is the longest path from the atoms. Nodes of level $ell$ read only results of levels $< ell$, so each level is internally independent — many kernel calls with disjoint result buffers — while level $ell+1$ must wait for *all* of level $ell$.],
  )
]

#align(center)[
  #block(width: 88%, inset: 4pt)[
    #import "@preview/cetz:0.4.2": canvas, draw
    #canvas(length: 1pt, {
      import draw: *
      let node(x, y, w, label, fillc, txtc) = {
        rect((x - w / 2, y), (x + w / 2, y + 11), fill: fillc, radius: 1.5pt, stroke: 0.4pt + luma(140))
        content((x, y + 5.5), anchor: "center", text(size: 5.8pt, fill: txtc)[#label])
      }
      let edge(p, q) = line(p, q, stroke: 0.55pt + luma(110))
      let dedge(p, q) = line(p, q, stroke: (paint: luma(150), dash: "dashed", thickness: 0.5pt))
      // edges first (drawn under the nodes)
      edge((85, 11), (111, 26))            // A -> [A,B]
      edge((145, 11), (119, 26))           // B -> [A,B]
      edge((85, 11), (58, 52))             // A -> [A,[A,B]]
      edge((108, 37), (60, 52))            // [A,B] -> [A,[A,B]]
      edge((145, 11), (177, 52))           // B -> [B,[A,B]]
      edge((122, 37), (175, 52))           // [A,B] -> [B,[A,B]]
      dedge((60, 63), (103, 82))           // into the middle band
      dedge((175, 63), (127, 82))
      dedge((108, 92), (100, 106))         // out of the band to the top
      dedge((122, 92), (134, 106))
      // funnel hint: light dashed guides, narrow at atoms and top, wide mid-levels
      dedge((70, -4), (18, 87)); dedge((18, 87), (82, 120))
      dedge((160, -4), (212, 87)); dedge((212, 87), (148, 120))
      // nodes
      node(85, 0, 26, "A", luma(150), white)
      node(145, 0, 26, "B", luma(150), white)
      node(115, 26, 34, "[A,B]", sw.lighten(30%), white)
      node(60, 52, 52, "[A,[A,B]]", sw.lighten(15%), white)
      node(175, 52, 52, "[B,[A,B]]", sw.lighten(15%), white)
      node(100, 106, 30, "[A,·]", sw.darken(10%), white)
      node(134, 106, 30, "[B,·]", sw.darken(10%), white)
      content((115, 87), anchor: "center", text(size: 11pt, fill: luma(80))[⋮])
      content((158, 87), anchor: "west", cap(6pt)[\ ~150 nodes/level])
      // level markers
      content((22, 5.5), anchor: "east", cap(6pt)[atoms])
      content((22, 31.5), anchor: "east", cap(6pt)[level 1])
      content((22, 57.5), anchor: "east", cap(6pt)[level 2])
      content((22, 87), anchor: "east", cap(6pt)[⋮ mid])
      content((22, 111.5), anchor: "east", cap(6pt)[top])
    })
    #v(3pt)
    #cap(6.3pt)[The fold as a DAG — the $[A,B]$-iterate spine of the truncated BCH structure, shown schematically. Solid edges: consumption of an operand (a node's two operands are results of earlier levels and/or the atoms — note how $[A,B]$ is *shared* by every node above it: interning computes it once). The dashed guides suggest the *funnel* shape of the real DAG at working scale (2×12: 11 levels, ~600 nodes, ~150 nodes at the widest level): thin at the atoms and the deepest levels, wide in the middle. The levels are the longest path from the atoms; width per level is what a parallel schedule can harvest, and the barrier after each full level is what it must pay.]
  ]
]

=== One parallel dispatch for the whole fold: the pack-slot walk

Naively, a fold with $L$ levels pays $L$ separate fork/join dispatches — measured at ≈35% of a small fold's time. The folded engine instead plans *all* levels up front (legal because every node's support is fixed for the fold's duration) and dispatches *once*. Dependency order is enforced not by dispatch boundaries but by a small counter protocol:

#text(size: 8.2pt)[
  #grid(columns: (84pt, 1fr), column-gutter: 6pt, row-gutter: 3.5pt,
    align(left)[*Packs per level*], align(left)[Each level's (node, unit) tasks are cut into work-balanced packs — again, only at unit boundaries (§3). Because each node is one anagram class and tasks are node-ordered, every pack is *class-contiguous for free*, which keeps a pack's scatter traffic inside a few cache-resident class blocks.],
    align(left)[*Slots*], align(left)[A small number of long-lived worker tasks — never more than workers — each *claims packs off an atomic cursor* per level, executes them, then bumps the level's completion counter.],
    align(left)[*Gates*], align(left)[At each level boundary *every* slot waits for the previous level's counter to reach its pack count (acquire/release ordering). A slot that arrives early is not allowed to sit idle-*and*-blocking: before waiting, it *drains any unclaimed packs of the level itself*. That one rule is the whole liveness proof: a level's remaining work is always either claimed-and-running or claimable by the waiter, so the walk cannot deadlock no matter how the pool is oversubscribed.],
  )
]

#align(center)[
  #block(width: 92%)[
    #grid(columns: (30pt, 58pt, 58pt, 58pt, 58pt), column-gutter: 2pt, row-gutter: 2pt,
      [], align(center, cap(6.2pt)[level 0 (thin)]), align(center, cap(6.2pt)[level 1 (wide)]), align(center, cap(6.2pt)[level 2 (wide)]), align(center, cap(6.2pt)[level 3 (thin)]),
      align(right + horizon, cap(6.2pt)[slot 0]), boxc(58pt, 13pt, sw, "pack α", white), boxc(58pt, 13pt, sw.lighten(15%), "pack β", white), boxc(58pt, 13pt, sw.lighten(25%), "pack γ", white), boxc(58pt, 13pt, sw.darken(15%), "pack δ", white),
      align(right + horizon, cap(6.2pt)[slot 1]), boxc(58pt, 13pt, rs.lighten(0%), "wait→drains β", white), boxc(58pt, 13pt, sw.lighten(35%), "pack ε", black), boxc(58pt, 13pt, sw.lighten(45%), "pack ζ", black), boxc(58pt, 13pt, idl, "— done", luma(120)),
      align(right + horizon, cap(6.2pt)[slot 2]), boxc(58pt, 13pt, rs, "wait (parked)", white), boxc(58pt, 13pt, sw.lighten(50%), "pack η", black), boxc(58pt, 13pt, idl, "parks", luma(120)), boxc(58pt, 13pt, idl, "—", luma(120)),
    )
    #v(2pt)
    #cap(6.3pt)[The walk across one fold (3 slots, 4 levels; shaded = sweep work, red = waiting). Barriers are implicit in the counters; empty levels or low-work phases are where slots park (and free their workers for co-scheduled work) after draining the cursor. Per-word accumulation order is exactly the serial schedule's: packs never split a unit.]
  ]
]

== The batch: amortizing everything except the serial dependence

A path of 1000 segments is 1000 folds. Running them through the per-fold machinery would re-plan, re-gate, re-gather, and re-allocate 1000 times. The batched driver instead recognizes what is invariant *across* folds and pays it once:

#text(size: 8.2pt)[
  #grid(columns: (96pt, 1fr), column-gutter: 6pt, row-gutter: 3.5pt,
    align(left)[*Support fixed point*], align(left)[After the first few folds the accumulator is dense; every node's "which words does it write" list becomes a fixed point. From then on: no per-fold gating, no list rebuilds — the planner's whole output (gates, tickets, packs) is reused verbatim for the remaining folds. Eligibility for batching is exactly this fixed point plus a shared displacement support.],
    align(left)[*Class-space accumulator*], align(left)[The accumulator stays in the class-contiguous basis order across the whole batch (gathered once, written back once at the end), removing two full-basis permutations per fold.],
    align(left)[*Compact node buffers*], align(left)[Each node's result lives in a private scratch buffer holding only the degree slices it can touch (degree-blocked, shift-addressed), not the full basis — hundreds of nodes × full-basis buffers would be megabytes and cache-hostile; degree slices cut that ~5×.],
    align(left)[*One stage chain, one walk*], align(left)[The whole batch is a single chain of stages — per fold: *gather* the displacement (block-parallel), the $L$ sweep levels, then *accumulate* the terms into the accumulator (block-parallel, simultaneously re-zeroing the node buffers for the next fold) — driven by exactly one pack-slot walk, sharing all cursors and counters.],
  )
]

#align(center)[
  #block(width: 100%)[
    #stack(dir: ltr, spacing: 1.5pt,
      cap(6.2pt)[lead~], boxc(28pt, 13pt, gv.lighten(30%), "zero", black),
      cap(6.2pt)[~fold 1:~], boxc(38pt, 13pt, gv, "gather₁", white),
      boxc(11pt, 13pt, white, "→", luma(80)),
      ..range(5).map(i => boxc(17pt, 13pt, sw.lighten(i * 8%), "L" + str(i+1), white)),
      boxc(11pt, 13pt, white, "→", luma(80)), boxc(48pt, 13pt, ac, "accumulate₁", white),
      cap(6.2pt)[~fold 2:~], boxc(38pt, 13pt, gv, "gather₂", white),
      cap(6.2pt)[~…],
    )
    #v(2.5pt)
    #cap(6.3pt)[The batch as one stage chain (5 of $L$ sweep levels shown). Parallel work exists *inside* every stage; between stages, barrier counters; between the last accumulate of fold $k$ and the first sweep of fold $k+1$, a true data dependency — the accumulator. The chain is long and thin by design: the planning and allocation it replaces cost far more than the barrier hops when a fold's sweep work is large, and far less (§7) when it is small.]
  ]
]

The accumulate stage deserves one sentence of its own, because it is the batch's second-largest cost: each fold's terms must be added (with BCH weights) into the accumulator, and every node buffer re-zeroed, per fold. It runs block-parallel over class-position ranges — but measured element volume is large (≈575,000 slot-touches per fold at 2×12, ≈85% of the sweep's own multiply count) because node buffers store whole degree slices while only support positions are live. Trimming the walk to live slots and fusing read-and-zero into one pass made this $-63\%$ smaller; what remains is *coherence-limited*: per-element cost grows ~9× from 4 to 32 threads as the shared accumulator region crosses many cores — the one remaining place where more threads actively hurt inside a stage.

== Scheduling to the shape of the work

The one structural number that decides whether parallelism helps is *planned sweep work per fold* — visible exactly, for free, inside the plan. The engine therefore sizes its own parallelism:

#align(center)[
  #block(width: 97%, breakable: false)[
    #grid(columns: (96pt, 92pt, 96pt, 1fr), column-gutter: 5pt, row-gutter: 3pt,
      [], align(center, cap(6.6pt)[*planned entries per fold*]), align(center, cap(6.6pt)[*chosen slots (32-thread pool)*]), align(center, cap(6.6pt)[*e2e at 32t vs 1t*]),
      [12×2 (high-dim, low-degree)], [66], [1 (serial walk)], cap(6.6pt)[2.55 → 1.25 ms ✓ \ (was 16.4 ms!)],
      [8×3], [700], [1], cap(6.6pt)[28.6 → 7.7 ms ✓],
      [3×8], [≈74,000], [19], cap(6.6pt)[338 → 161 ms ✓],
      [2×12 (low-dim, high-degree)], [≈257,000], [32], cap(6.6pt)[1393 → 549 ms ✓],
    )
  ]
]

#align(center)[
  #block(width: 88%)[
    #grid(columns: (60pt, 62fr, 26fr, 8fr, 4fr), column-gutter: 2.5pt, row-gutter: 2.5pt,
      [12×2 \@32t], boxc(100%, 11pt, rs.lighten(30%), text(size: 5.6pt, fill: white)[barrier/sync 92%], white), boxc(100%, 11pt, sw, "sweep", white), boxc(100%, 11pt, ac, "", white), boxc(100%, 11pt, white, "", white),
      [2×12 \@32t], boxc(100%, 11pt, rs.lighten(30%), text(5.6pt, white)[sync], white), boxc(100%, 11pt, sw, text(5.6pt, white)[sweep], white), boxc(100%, 11pt, sw.lighten(20%), text(5.6pt, white)[…], white), boxc(100%, 11pt, ac, text(5.6pt, white)[acc], white),
    )
    #cap(6.3pt)[Same engine, two shapes, before work-adaptivity (illustrative shares, measured): at 12×2 the fold's 66 entries vanish inside the stage machinery (16t = 6.6× *slower* than 1t); at 2×12 the sweep dominates and 16–32 threads pay. The scheduler's rule — *slots ≈ planned work ÷ a fixed quantum, capped by the pool, floored at 1* — makes each shape run in the regime that suits it, at the cost of nothing when the work is large.]
  ]
]

This single rule closed the high-dimension/low-degree pathology: the parallel infrastructure is still there, but a 66-entry fold no longer pays for it.

== What threads cannot reach inside a fold — and the associative axis they can

It is worth stating plainly where parallelism *does not* exist, because it is a consequence of the mathematics rather than an implementation gap:

#block(width: 100%, inset: 6pt, fill: luma(248), radius: 2pt)[
  *The fold chain is serial along the path.* Fold $k+1$ evaluates its DAG on the output of fold $k$; level-0 jobs of fold $k+1$ read the accumulator written by the accumulate stage of fold $k$. No reordering of stages, no matter how clever, removes a *data* dependency. Everything shown so far — units, levels, packs, slots, batches — lives strictly inside one concatenation, replicated identically 1000 times down the chain.
]

The one honest escape hatch already exists in the algebra: BCH concatenation is *associative*. The path's log signature is an associative reduction of the segment displacements, so the reduction tree may have any shape — including a balanced tournament:

#align(center)[
  #block(width: 88%)[
    #grid(columns: 1fr, row-gutter: 2.5pt,
      grid(columns: (1fr,)*8, column-gutter: 1.5pt, ..range(8).map(i => boxc(100%, 11pt, luma(170), "d" + str(i+1), white))),
      grid(columns: (1fr,)*8, column-gutter: 1.5pt, ..range(4).map(i => boxc(100%, 12pt, sw.lighten(35%), "⋈", black)) + range(4).map(i => boxc(100%, 12pt, white, "", white))),
      grid(columns: (1fr,)*8, column-gutter: 1.5pt, ..(range(2).map(i => boxc(100%, 12pt, sw.lighten(15%), "⋈", white)) + range(6).map(i => boxc(100%, 12pt, white, "", white)))),
      grid(columns: (1fr,)*8, column-gutter: 1.5pt, ..((boxc(100%, 12pt, sw.darken(5%), "⋈", white),) + range(7).map(i => boxc(100%, 12pt, white, "", white)))),
    )
    #cap(6.3pt)[A tournament reduction over 8 segments: $log_2 n$ sequential rounds, every pair in a round independent. Rounds start *sparse* (degree-1 displacements: cheap, support-gated folds) and only the last $~2 log_2 n$ concatenations along the critical path pay full dense-fold cost. A 1000-segment path: from 1000 serial dense folds to ~10 of dependency depth, up to 500 concurrent in round 1.]
  ]
]

The engine now performs exactly this reduction, adaptively and with no new API surface: when a batch holds at least 16 displacements on a multi-thread pool, the batch driver reduces a *balanced tournament over contiguous chunks* instead of the serial chain. Leaves are not cold per-fold dispatches — each leaf folds its chunk through the same sequential batch engine (warm-up folds, then one steady batched dispatch). Cold per-fold leaves forfeited that amortization and regressed 2×12 by up to 1.27× in measurement; with batched leaves the per-leaf warm-up duplication is the only real tax, and the merge tail — the only remaining serial spine — is ~$2 log_2 (n\/8)$ dense folds deep. The tree shape is a pure function of the displacement count (at most 32 leaves), so results are *bit-stable at any pool size*; exact (rational) coefficients remain bit-identical to the sequential driver, while f64 reassociation shifts rounding only in the last ulps (measured ≈$10^(-13)$ abs on $\~10^2$ magnitudes — accepted, and pinned by tolerance tests against the sequential oracle plus exactness oracles for rationals). The §7 discipline applies at the tree level too: leaves and merges run through the same work-adaptive machinery, so a thin regime's early rounds do not pay for scheduling. Measured end-to-end (full grid in the regime map below): 12×2 1.84→0.50 ms \@8t ($6.7×$ at 10000 segments), 8×3 7.94→2.08 ms \@32t, 3×8 154→92 ms \@32t, and 2×12 — the regime that *needs* the intra-fold machinery most — 525→390 ms \@16t ($2.0×$ at 10000 segments): ~128 batched leaves plus 127 merges replace 1000 dense folds. Single-threaded runs are unaffected: the driver keeps the sequential chain there, so the reformulation is upside-only on every measured grid.

The tournament's folds then compose with a second multiplier, *cohort execution*. A round's folds all walk the same data-independent plan, so the engine executes four of them as one SIMD walk over lane-interleaved buffers:

#align(center)[
  #block(inset: 4pt)[
    #import "@preview/cetz:0.4.2": canvas, draw
    #canvas(length: 1pt, {
      import draw: *
      let lanes = (sw.lighten(45%), sw.lighten(20%), sw.darken(5%), sw.darken(30%))
      let slot(x, y, w, h, lab) = {
        for l in range(4) {
          rect((x, y + l * h / 4), (x + w, y + (l + 1) * h / 4), fill: lanes.at(3 - l), stroke: 0.35pt + luma(120))
          content((x + w / 2, y + (l + 0.5) * h / 4), anchor: "center", text(size: 5pt, fill: white)[f#(4 - l)])
        }
        content((x + w / 2, y + h + 2), anchor: "south", cap(6pt)[#lab])
      }
      // shared plan strip with visible walk steps
      rect((16, 84), (258, 96), fill: luma(238), radius: 1.5pt, stroke: 0.4pt + luma(140))
      for (i, cx) in ((0, 55), (1, 129), (2, 203)) {
        rect((cx - 33, 86.5), (cx + 33, 93.5), fill: white, radius: 1pt, stroke: 0.35pt + luma(150))
        content((cx, 90), anchor: "center", text(size: 5.4pt)[ $(u_#(i + 1), v_#(i + 1), w_#(i + 1)) arrow e_#(i + 1)$ ])
      }
      content((16, 80), anchor: "west", cap(5.8pt)[ONE walk of the shared plan — decomps, tickets, scatter indices, identical for the whole 4-fold group:])
      // the slot towers
      slot(36, 26, 22, 32, [$A_4[u_k]$])
      slot(110, 26, 22, 32, [$B_4[v_k]$])
      slot(204, 26, 22, 32, [$"out"_4[e]$])
      content((84, 46), anchor: "center", text(size: 12pt)[$times$])
      content((168, 46), anchor: "center", text(size: 12pt)[$+$])
      content((84, 32), anchor: "center", text(size: 5pt)[mul — then add\ (vertical per lane,\ no FMA)])
      content((168, 32), anchor: "center", text(size: 5pt)[accumulates,\ per lane])
      // dashed guides from plan strip down to slots
      line((47, 82), (47, 62), stroke: (paint: luma(150), dash: "dashed", thickness: 0.5pt))
      line((121, 82), (121, 62), stroke: (paint: luma(150), dash: "dashed", thickness: 0.5pt))
      line((215, 82), (215, 62), stroke: (paint: luma(150), dash: "dashed", thickness: 0.5pt))
      // the trick annotations
      content((47, 16), anchor: "center", cap(5.6pt)[ONE contiguous\ 4-wide load, all lanes])
      content((121, 16), anchor: "center", cap(5.6pt)[ONE contiguous\ 4-wide load, all lanes])
      content((215, 16), anchor: "center", cap(5.6pt)[ONE contiguous\ 4-wide store —\ 4 folds written at once])
      content((6, 42), anchor: "west", cap(5.6pt)[lane = one\ tournament\ fold\ (f₁…f₄)])
    })
    #v(2pt)
    #cap(6.3pt)[Cohort execution of one 4-fold tournament group. Coefficient buffers are interleaved per index (`[idx*4+lane]`): every index access of the shared plan walk touches all four lanes in ONE contiguous 4-wide load/store — no gather/scatter instructions anywhere; arithmetic is vertical per lane (mul, then add — deliberately *no* FMA, proven worthless in this load/scatter-bound phase). Indices, gating tickets, and control flow are shared; only coefficient values differ per lane. Each lane therefore replicates the scalar fold's operation order exactly, so a cohort run is *0-ulp bit-identical* to the scalar tournament over the same reduction tree. The kernel is capability-gated like the existing raw-float fast path: f64/f32 (including `NotNan`-wrapped) via TypeId + runtime AVX2 check; every other coefficient type keeps the scalar kernel bit-for-bit. Engagement is measured-out (≤16 threads, or ≥256-term folds, or enough groups; env off switch). A/B: −32% 2×12\@4t, −37% 3×8\@4t, −29% 8×3\@4t (the old 32-thread time on 4 threads), −20% 2×12\@8t; neutral elsewhere. Absolute bests unchanged — cohorts buy that 4–8 threads now reach ~1.5× of the 16–32-thread times.]
  ]
]

#pagebreak()

== Regime map: what was measured

#text(size: 7.7pt)[
  #grid(columns: (50pt, 82pt, 76pt, 72pt, 1fr), column-gutter: 4.5pt, row-gutter: 2.2pt,
    [*grid*], [*stock, best config*], [*today, best config*], [*gain*], [*why*],
    [2×12], [1393 ms \@32t], [390 ms \@16t], [3.6×], [Every intra-fold layer engaged — plus the tournament: 128 batched leaves + 127 merges replace 1000 dense folds ($2.0×$ alone at 10000 segments). #h(0.5em)_(low-d, high-m)_],
    [3×8], [331 ms \@8t], [92 ms \@32t], [3.6×], [Same; tournament gain grows with pool ($1.7×$ \@32t).],
    [12×2], [2.55 ms \@1t (16t: 16.7 ms!)], [0.50 ms \@8t], [5.1×, pathology gone — and threads come back], [66 entries/fold: adaptivity runs the fold serial — the tournament puts map lanes on the fold AXIS ($6.7×$ vs stock at 10000 segments). #h(0.5em)_(high-d, low-m)_],
    [8×3], [8.02 ms \@2t (32t: 28.6 ms)], [2.08 ms \@32t], [3.9×], [700 entries/fold: formerly on the serial/parallel boundary — same mechanism as 12×2.],
  )
]
#cap(6.6pt)[End-to-end log signatures of 1000-segment random paths, criterion medians, 32-core machine. "Stock" = the engine before the parallelism work documented in the companion notes. Single-threaded runs keep the sequential driver (bit-identical there); parallel paths reassociate floating-point folds — exact types unchanged, ulps pinned by oracle tests.]

#v(4pt)
#block(width: 100%, fill: luma(245), inset: 6pt, radius: 2pt, text(size: 8pt)[
  *Summary.* Parallelism in this engine is *structural*, not opportunistic: the algebra partitions work into anagram units that never share a write; the fold DAG partitions those units into dependency levels; the batch amortizes all support-invariant planning; the scheduler spends threads only where measured work justifies the barriers. The fold chain along the path is a data dependency of the BCH recursion itself — and the tree reduction of §8, now performed adaptively by the batch driver, is the only mechanism mathematics offers for crossing it.
])
