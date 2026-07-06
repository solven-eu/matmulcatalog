# Project conventions for AI assistants

## Scheme registration — be systematic

This repo aims to be a **research-grade catalog of fast matrix
multiplication algorithms**. The historical and provenance value comes
from coverage, not selectivity.

**Whenever a new scheme is presented or discovered** (in a paper, a
chat, a search run), register it:

1. **If explicit factor matrices are available** → add a JSON file under
   `src/main/resources/schemes/sectionN/` following the existing naming
   convention `{source}-{provenance}_{n}x{m}x{p}_r{rank}_a{additions}[_field].json`.
   Examples: `alphatensor-Z_2x3x3_r15_a58.json`, `perminov-ZT_5x5x5_r93_a843.json`.
   Make sure the JSON encodes the field (e.g. `"field": "Z"` or filename suffix
   `_F2_` / `_C_` / `_Z_`). (A `_ZT_` token is Perminov provenance for a *ternary*
   `Z` scheme, not a separate field — see the ZT note under "Field discipline".)
   Run `Verifier.isExactNonCubic` if cubic / non-cubic
   to validate the scheme actually computes matmul.
2. **If only the rank claim is published (no explicit matrices)** → add an
   entry to `docs/cited-bounds.json` via `GenerateCitedBounds.java` with
   the source citation. The bound becomes visible in the catalog without
   pretending we have the scheme.
3. **If the scheme is derivable from a formula** (Waksman, Rosowski,
   Pan TA, etc.) → wire the formula into the matching `*Bound.java`
   helper (e.g. `WaksmanBound`, `PanTrilinearAggregation`) and have
   `GenerateDerivedBounds.java` emit it to `docs/derived-from-cited-bounds.json`.
4. **If the scheme is non-bilinear or otherwise doesn't fit `SchemeIO`**
   (e.g. Rosowski 2019's Algorithm 1) → describe it in markdown under
   `references/` and reference it from `REFERENCES.md`. Wire the bound
   if possible.

Always cite the source. Always state the field and whether the scheme
exploits commutativity.

**Don't be selective about "improvements only"**: registering a
historical result with worse rank than current SOTA is still
valuable — it preserves the chronology of progress and prevents
re-discovery being misattributed.

## When a bound is missing — check FMM and Perminov first

Before declaring a `⟨n,m,p⟩` rank "unknown" or trying to derive a
fresh bound, check both community catalogs (non-commutative,
both publish factor matrices):

- **[FMM (Université de Lille)](https://fmm.univ-lille.fr/)** —
  per-format pages at `https://fmm.univ-lille.fr/{n}x{m}x{p}.html`
  (sorted n≤m≤p convention). Maple format
  (`*.mpl.bz2`), needs bunzip2 + parse.
- **Perminov FastMatrixMultiplication** —
  [status.json](https://github.com/dronperminov/FastMatrixMultiplication/blob/master/schemes/status.json)
  is the synthetic index with per-format `ranks` + `schemes` (per
  field Q / Z / ZT). Each entry has a `source` pointing to the
  originating sub-catalog (classic / alpha_tensor / alpha_evolve /
  a_60_addition / fmm_add_reduction / jakobmoosbauer_* /
  meta_flip_graph / tensor) — propagate that as
  `attribution_for_rank` per the "Distinguish discoveries from
  re-discoveries" rule above.

Whichever has the best (= lowest) rank wins. Perminov re-encodes
many FMM-known schemes plus his own derivations, so it usually
suffices to pull from Perminov alone — but cross-check FMM for
recent additions and for schemes Perminov doesn't carry.

## Distinguish discoveries from re-discoveries

A scheme imported under a given source (e.g. `alphatensor-Z_*`) is
NOT automatically a discovery by that source. The AlphaTensor 2022
paper (Fawzi et al., Nature) reports per-format ranks but only some
are colour-highlighted as improvements over prior SOTA — the rest
matched known bounds (Smirnov 2013, Strassen recursion, Pan, etc.).

**When importing schemes in bulk**, capture the discovery status:

- Add a `"discovery"` field to the JSON: `true` | `false` | `"TBD"`
  (default `"TBD"` if not verified against the source paper).
- If `discovery: false`, add `"attribution_for_rank"` pointing to the
  earliest source that established the bound (paper, OEIS, etc.).
- Do NOT cite the importing source (e.g. AlphaTensor) for ranks it
  merely recomputed — that misattributes the historical record.

This applies retroactively: existing `alphatensor-*` and similar bulk
imports should be audited and re-tagged. Track this in
`research/DISCOVERIES_PENDING_ANALYSIS.md` until done.

## Tracking provisional findings

When a search run produces a result that hasn't been digested into the
permanent catalog yet, add it to `research/DISCOVERIES_PENDING_ANALYSIS.md`
with: claim, repro command, breakdown, status. Don't let findings
evaporate when a process dies.

## Shape notation — always tag the algebra

A bare `⟨n,m,p⟩` is **ambiguous** — different algebras give wildly
different ranks (canonical example: `⟨4,4,4⟩` is 47 over F₂, 48 over
C/Q/R since the DPS-2025 rationalisation of AlphaEvolve, 49 over Z —
pre-DPS, R itself stood at 49). Our shape notation always includes the algebra context
in prose and code-comments. (This direction was informally validated
in a conversation with Sedoglavic — he confirmed it's a reasonable
way to go; the convention is ours, not his decision.)

**Preferred forms:**

| Form | Meaning | Example |
| --- | --- | --- |
| `K⟨n,m,p⟩` | Non-commutative matmul over field `K` | `R⟨7,7,7⟩=250`, `F2⟨4,4,4⟩=47` |
| `K⟨n,m,p⟩ᶜ` (superscript-c) | Commutative matmul over `K` | `R⟨3,3,3⟩ᶜ=21` (Rosowski) |
| `K⟨n,m,p⟩:m` | "Multiplications" count (small-matrix focus) | Sedoglavic convention |
| `K⟨n,m,p⟩:r` | Rank (asymptotic-complexity focus) | Equal to `:m` for fixed-size schemes |
| `K⟨n,m,p⟩:a` | Addition count | `R⟨7,7,7⟩:a=2417` |

**Why `:r` vs `:m`**: `r` (rank) is the right term when discussing
asymptotic complexity (`ω = log_n(r)`). `m` (multiplications) is the
right term when discussing concrete small-matrix algorithms. Sedoglavic
uses `m`; we adopt `:m` in new code and gradually migrate filenames
(`_r{count}_` → `_m{count}_`). For now both are accepted in
`FieldAwareLookup`.

**Filenames** already encode the field as a suffix token (`_Z_`,
`_Q_`, `_R_`, `_C_`, `_F2_`). Commutative-only schemes carry
`"commutative": true` in JSON metadata. The shape notation above is for
prose and code-comments, not for filename layout.

**ZT is NOT a field** (load-bearing — it was historically, wrongly,
conflated with F₂/Z₂; it has nothing to do with characteristic 2). `ZT`
is the **sub-class of `Z`** (integer) schemes whose every U/V/W
coefficient is in `{-1,0,1}` — **"ternary integer"** (Perminov's term).
Beware the word "ternary": `F3`/`Z3` (GF(3)) is also ternary, but
**ternary _modular_** (three residues `{0,1,2}` in characteristic 3) —
a different algebra that merely shares the "three values" count. `ZT` is
ternary _integer_ (`{-1,0,1}` over characteristic-0 ℤ). It is carried
as a per-scheme boolean — `"ZT": true|false` in `catalog.json`, emitted
by `GenerateCatalogManifest` **only when `Z` is among the scheme's
`fields[]`**, and computed from the actual coefficients
(`SchemeIO.isTernary`), NOT from the `_ZT_` filename token. So a ternary
scheme from any source (Strassen, AlphaTensor-Z, …) is `ZT:true`
regardless of its filename. The `_ZT_` filename token survives only as
Perminov provenance and classifies to `Field.Z`. In the SPA, "ZT" is a
field-selector value that filters `Z` schemes by `scheme.ZT === true`.

## Field discipline (load-bearing)

Every rank claim **must** name its field. The ⟨4,4,4⟩ landscape is the
canonical why: 47/F₂ (AlphaTensor 2022), 48/C (AlphaEvolve 2025) = 48/Q/R
(Dumas–Pernet–Sedoglavic 2025 rationalisation, arXiv:2506.13242), 49/Z
(Strassen²). Field-separation claims are date-stamped: between AlphaEvolve
and DPS, R stood at 49 while C had 48.

**Covered fields** (the algebras a scheme's `fields[]` can name). The
canonical user-facing version is the "Fields covered" legend in the SPA
(`docs/index.html`); keep the two in sync.

| Tag | Algebra | Coefficients | Note |
| --- | --- | --- | --- |
| `F2` (=`Z2`) | GF(2) | {0,1} mod 2 | Characteristic 2; own universe, no char-0 inclusion. |
| `F3` (=`Z3`) | GF(3) | {0,1,2} mod 3 | **Ternary _modular_** (characteristic 3). |
| `Z` | ℤ | any integer | Characteristic 0; `Z ⊂ Q ⊂ R ⊂ C`. |
| ↳ `ZT` | *sub-class of Z, not a field* | {−1,0,+1} ⊂ ℤ | **Ternary _integer_**; per-scheme boolean (see ZT note above). |
| `Q` | ℚ | any rational | Characteristic 0; implies R, C. Reduces mod prime `p` iff every denominator is coprime to `p` (`num·den⁻¹ mod p`): `1/2` works in F₃ (2⁻¹≡2) but not F₂; `1/3` works in F₂ (3≡1) but not F₃. Verified per scheme by `Verifier.residualNonCubicFp`; `NarrowFields` stamps the resulting F₂/F₃. |
| `R` | ℝ | any real | Characteristic 0; implies C. |
| `C` | ℂ | any complex | Accepts any characteristic-0 scheme (Z⇒Q⇒R⇒C inclusion); plus C-only schemes (e.g. AlphaEvolve's complex-coefficient originals — note its ⟨4,4,4⟩=48 rank is no longer C-only since the DPS-2025 rationalisation made 48 Q/R-valid). |

- **`Z`, `Q`, `R` are distinct, first-class fields** — there is **no
  "R/Q/Z cluster"**. Do not merge them anywhere: catalog `fields[]`,
  comparison tables, attribution, and bounds all keep them separate.
  `ZT` (ternary-integer) is likewise **first-class** (the per-scheme
  `{−1,0,1}` sub-class of `Z`). The only inclusion we use is for
  **field computation** — a `Z` scheme is also valid over `Q`, then `R`,
  then `C` (`Z ⇒ Q ⇒ R ⇒ C`), as exploited by field-widening / fallback.
  And the only place `Z`/`Q`/`R` are presented as a single grouped
  choice is the **SPA field selector**, where the user may pick any one.
  (All three are characteristic-0 and non-commutative-friendly, so they
  lift to recursive matmul — but that shared property is not a licence to
  fold them together.)
- **Sweeps default to `--field=Q`** (user requirement). All char-0 NC
  search/sweep/gap-closing work runs over `Q` — the FMM/Perminov digest
  scope; a `Q` lookup admits Z+Q ingredients. R and Q are often the same
  in practice (the catalog holds ~1 R-only scheme), but `Q` is required
  by default. Pass `R` only when R-only ingredients are deliberately
  wanted — and treat an R-only *derived* stub (a composition of Q/Z
  atoms stamped `["R","C"]`, e.g. `17x19x20-r3780`) as field drift to
  repair with `NarrowFields`, not a scheme to build on.
- **C**: complex extension; R-valid schemes work; plus AE 48 for ⟨4,4,4⟩.
- **F₂**: characteristic 2. AlphaTensor results. Integer (Z-exact) schemes
  also reduce here (mod 2) — and mod 3 to F₃ — as a theorem; the stampers
  (`StampFields`, `NarrowFields`) both grant F₂/F₃ from Z, and the SPA
  surfaces it purely via `fields[]` membership. (Historically `StampFields`
  omitted F₂/F₃ from Z, leaving ~8k integer schemes invisible under the
  F₂/F₃ selectors until re-stamped.)
- **F₃ vs ZT** — both are "ternary" but unrelated: F₃ is ternary
  _modular_ (GF(3)), ZT is ternary _integer_ ({−1,0,1} over ℤ).
- **Commutative**: a separate axis. Commutative-only schemes (Waksman,
  Rosowski, Makarov 1986) do **NOT** lift to recursive matmul over NC
  rings, but ARE valid for scalar matmul. Always tag with
  `"commutative": true` in JSON and bounds entries.

When comparing to historical tables (DIS09 Table 3 vs Table 4), pick
the right SOTA pipeline — `TestDIS09FullScan` for NC, `TestDIS09Table4Commutative`
for commutative. Cross-contamination silently produces wrong "wins".

## Scheme JSON: content-driven, canonically formatted, filename-as-label

The catalog was migrated (2026-06) so **filenames are pure labels** and **every
scheme property is read from JSON content, never parsed from the filename**.

- **Filename convention** (cosmetic only): `{n}x{m}x{p}-r{rank}-{note}-{hash7}.json`
  (e.g. `3x3x3-r23-laderman_1976-1a2b3c4.json`). `note` = source/author for imports,
  `derived` / `derived_strassen_recombine` / … for our output; `hash7` = first 7 of
  the content hash. Renaming a scheme file is safe — nothing reads the name.

- **Read metadata from CONTENT, never the filename**: shape ← `n`, rank ← `m`,
  field ← `fields[]`, source ← `source`, additions ← `additions`, commutative ←
  `commutative`. The lookup (`FieldAwareLookup.buildIndex`) and manifest
  (`GenerateCatalogManifest`) are content-driven — do NOT re-introduce a `NAME`-regex
  or `classifyFilenameField`-style heuristic. Lineage refs are `shape@hash7` (pinned)
  or bare `shape` (resolves to catalog-best), **never** filename stems. New schemes
  must carry `fields`/`source`/`additions` in JSON (the materialiser stamps them;
  `StampFields`/`StampSource`/`StampAdditions` backfill).

- **Write scheme JSON with `MatrixJsonFormatter.format(node)`** — what
  `SchemeIO.write` / `updateFields` emit. NEVER `toPrettyString` / Jackson default
  pretty-printing: it drifts the format (spaces around colons, inline matrix rows)
  and produces noisy diffs. `ReformatSchemes` re-canonicalises the whole tree if drift
  creeps back in.

## `docs` package layout — put new drivers in the right sub-package

`eu.solven.matmul.docs` holds `main()`-only drivers (no library code is
imported from non-`docs` production paths). They are split by role —
**place every new driver in the matching sub-package**, never loose in
`docs/`:

- **`docs.generate`** — regenerators of committed artifacts
  (`GenerateCatalogManifest`, `GenerateCitedBounds`, `GenerateDerivedBounds`,
  `GeneratePaperTables`, `EnrichSchemeMetrics`, …). Re-run to refresh
  `docs/*.json` / paper tables.
- **`docs.verify`** — validators / audits / cross-catalog comparisons
  (`VerifyAllSchemes`, `AuditLineageRefs`, `DetectCyclicStubs`,
  `FmmCrossCheck`, `SyncReferenceCatalogs`, `SanityCatalogMigration`, …).
- **`docs.migrate`** — one-shot catalog mutations / stampers / fixers
  (`RenameSchemes`, `RehashRationalsComplex`, `ReformatSchemes`,
  `Stamp*`, `NarrowFields`, `FieldWideningSweep`, `Fix*`, `RetireDerivableImports`,
  `Materialise{ZT,Kron888_336,Concat366_80}`, …). Already-run, but kept +
  re-runnable when drift recurs.
- **`docs.explore`** — throwaway probes / reports / experiments
  (`AllocationDeepDive`, `BaseFingerprintLister`, `Probe*`, `Scan*`,
  `Craft*`, `ProjectFmmGaps`, `RankBasesByOmega`, …).
- The search **engines** (`SchemeSweep`, `SerendipitousSweep`,
  `Sweep2x2x2AllSupports`) stay at the top level of `docs`.

White-box tests live in the same (sub-)package as the class under test
(e.g. `TestFieldWideningSweep` → `docs.migrate`).

## Tests use AssertJ

New / edited tests use `assertThat(...)`, not JUnit `Assert.*`.

## Every high-level / silent bug gets a regression test

When you fix a bug that was **high-level** (a whole search/sweep path silently
degraded, not a local off-by-one) and especially one that failed **silently**
(no crash — just worse results, an empty pool, a skipped scheme), you **must**
add a fast, parameter-specific test that would have caught it. The point is to
turn a one-time discovery into a permanent guard.

Canonical examples (all from real silent regressions in this repo):
- `extendedPool` returned empty after the known/derived/curated split
  (`listFiles("section*")` on the root found nothing) → search quietly weakened.
  Guard: `TestSweepSpotsSota.extended_pool_is_not_empty`.
- The bud-ordering bug hid ⟨8,9,9⟩=430 (U-first greedy masked the size-3 V-bud).
  Guard: `TestSerendipitousBudProduct.bud_ordering_recovers_8x9x9_430` +
  the SOTA-spotting rows in `TestSweepSpotsSota`.

Rules:
1. The guard must be **fast** — probe a handful of specific shapes
   (`RecursiveMaterialiser.materialise(n,m,p)` or the relevant engine entry),
   not the full `SchemeSweep` / `VerifyAllSchemes`. Build the lookup once.
2. Assert **SOTA-or-better** (`≤`), not equality — a real improvement must never
   break the test, only a regression. Add shapes that each exercise a *distinct*
   mechanism (disk hit / Kron / concat / recombination / serendipitous / …).
3. New SOTA-class shapes that were *hard to reach* (e.g. ⟨17,17,17⟩, the
   serendipitous family) belong in `TestSweepSpotsSota` so they can't silently
   regress later.

## Long-running procedures: progress logs + background execution

Anything that takes more than ~30s wall-clock should:

1. **Emit periodic progress lines** — every N units of work, log
   `[progress] X processed (… counters …) Yms elapsed`. With an ETA
   if the total is known: `[progress] X/N processed, ~Zs remaining`.
   The `MaterialiseRecursiveSweep` pattern is the canonical example.

   **Adapt the ETA to workload non-uniformity.** A naive linear
   extrapolation (`remaining_count × recent_per_item_time`) is wrong
   when later items are systematically heavier — e.g. the sweep
   iterates shape ⟨n,m,p⟩ with n≤m≤p, so the last batches process
   shapes whose matrices are 10–30× larger than the first batches.
   Either weight by predicted work (matrix entries, rank·dim, etc.),
   or fit the per-batch trend line (each batch takes ~constant more
   than the previous) instead of a flat throughput.

2. **Be launched in the background** when invoked by the assistant
   (`run_in_background: true` on the Bash tool), so:
   - the user can pivot to other work while CPU is busy;
   - the assistant can poll the output file periodically and keep
     the user posted (one short status update per poll, not silence);
   - the assistant doesn't burn its own context window waiting in a
     foreground command.

3. **Always enable heap-dump-on-OOM**. Any background driver MUST be
   launched with
   `-XX:+HeapDumpOnOutOfMemoryError -XX:HeapDumpPath=target/oom-dumps/`
   in `MAVEN_OPTS` (in addition to whatever `-Xmx` the workload needs).
   Without this, an OOM kill leaves us blind — exit code 144 with no
   evidence of *what* exhausted the heap. With it, a `.hprof` snapshot
   lands under `target/oom-dumps/` and can be opened in any heap-profiler
   (Eclipse MAT, VisualVM) to identify the retained-object culprit.
   Canonical invocation:
   ```bash
   MAVEN_OPTS="-Xmx4g -XX:+HeapDumpOnOutOfMemoryError -XX:HeapDumpPath=target/oom-dumps/" \
       mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.SomeBigDriver
   ```
   Create the dump directory ahead of time (`mkdir -p target/oom-dumps`)
   — `HeapDumpPath` is just a hint, not an mkdir.

   **Caveat: `HeapDumpOnOutOfMemoryError` only fires for JVM-thrown
   `OutOfMemoryError`.** If the JVM gets `SIGKILL`-ed by the OS
   (macOS / Linux killing the process because RSS / VM size exceeds a
   system pressure threshold), there's no Java-side OOM event and no
   dump is written. Symptom: `exit code 144`, empty `target/oom-dumps/`,
   process gone with no Java-stack trace. Cause is usually off-heap
   allocations (direct byte buffers from Jackson / NIO, large
   `double[][]` arrays from `CompactScheme.expand` on big schemes,
   retained `JsonNode` trees). Bumping `-Xmx` doesn't fix this — the
   heap was never the issue. The right move is to stream-process and
   drop references aggressively, or sample with `jcmd <pid> GC.heap_info`
   periodically to confirm where memory is going.

The combination means: the user always knows whether work is making
progress, how far along it is, and roughly when it will finish — and
if it dies, we have evidence rather than guesses.

### Process hygiene — kill what you started

Background JVMs / sweeps / monitors started by the assistant are the
assistant's responsibility to track and clean up. Don't let stale
processes accumulate.

Concrete discipline:
- **Track every background process** you launch. When you `run_in_background`
  a long task or arm a `Monitor`, remember its PID/task-id and what
  config it's running.
- **When the user redirects** (interrupt, re-scope, change shape, "let's
  focus elsewhere"), explicitly kill the previous run BEFORE starting
  the next one. Don't leave a 14-hour sweep CPU-pegged while you start
  a new 30-minute sweep.
- **When the user asks "interrupt all work"** or similar, `pkill` the
  target pattern AND verify with `ps aux | grep …` that nothing
  matching is still alive. Don't just trust the SIGTERM.
- **At the end of a thread** (after reporting results), check that no
  long-running JVM remains armed. If something is, name it explicitly
  to the user so they can decide.
- **Distinguish ours from the user's**: IntelliJ debug JVMs, Eclipse
  lemminx, IDE Maven helpers are NOT yours to kill. `ps aux` shows them
  too; identify by command line before pulling any trigger.

Failure mode to avoid: silent-CPU-burn while context moves on. A
1-hour sweep that finishes 30 minutes after the user has lost interest
in its result is wasted wall-clock AND fan noise AND power. Kill it
when the question changes.

## Math first, brute force only after you've understood the structure

When working on a rank or scheme-construction problem (HK same-method
pairs, residual-subspace decomposition, atom catalogs, …), start with
**math-driven analysis** — not enumeration. Concretely:

1. Write down the algebraic identity you're trying to satisfy and what
   the residual is (what's left after the obvious sub-products land).
2. Decompose the residual into invariant subspaces (rank, basis, image,
   kernel). Compute the *bilinear* rank — not just the joint matrix
   rank, since atoms can cross-mix residual subspaces.
3. From the rank, derive a lower bound on the number of new products
   needed AND on the shape of any valid atom (which row/col forms,
   which coefficients).
4. Only THEN write a sympy/Java enumerator — and constrain it by the
   shape derived above, not by "all atoms in some natural catalog".

The HK `(2,2,bridge-3)` case is the canonical why: my first claim
"joint matrix rank 4 → impossibility" was wrong because atoms with
cross-mixing row/col forms can contribute to multiple residual
subspaces simultaneously. The right invariant is *tensor rank*, which
the brute-force enumeration on its own didn't surface. ChatGPT
proposed the residual-decomposition framing in parallel and exposed
the mistake.

Brute-force enumeration is still useful — but as the **second step**
to confirm a constructive existence within a constrained shape, not
as the first step to discover the shape.

## Optimality discipline — local vs global, always labelled

Every computed minimum/maximum/"best" value **must** be labelled with
*what kind of optimum it is*. Heuristics and local minima are fine and
often necessary — but it must be **crystal clear, at every layer (code,
field names, JSON, logs, prose, the SPA, the paper)**, whether a number
is a proven global optimum or merely an upper/lower bound from a
restricted search.

**The three honesty tiers** (use these words):

| Tier | Meaning | How to surface it |
| --- | --- | --- |
| **proven-optimal** | Global optimum, exhaustively verified or proven (e.g. exact ILP/SAT on a small instance; de Groote uniqueness; an exhaustive enumeration that finished). | May be stated as "= optimum". |
| **bound** | Best found by a heuristic / bounded / greedy / anytime search — a valid **upper bound** (for ranks/additions) or **lower bound**, but *not* proven minimal. | Call it `min_*`/`best_*` and tag it a bound; never say "= optimum". |
| **optimal-within-scope** | Global optimum over a *restricted rule set / atom set / search space*, not over all constructions. | State the scope explicitly: "optimal over {rules}", not "optimal". |

**Concrete cases in this repo (all are NOT global optima):**

- **Additive complexity** (`LinearCircuitMinimizer` / `min_additions` /
  `slp`): a **greedy, cancellation-free** heuristic → an **upper bound**,
  not the minimal LSP (the problem is NP-hard). It matching Winograd's
  15 is luck on a small case, not a guarantee. Label `min_additions` as
  a bound; only a verified ILP/SAT result for small dims is
  proven-optimal.
- **Rank closure** (`SchemeSweep` / `RecursiveMaterialiser`): the
  fixpoint is **optimal-within-scope** — global only over the current
  rule set (Kron / concat / recombination) + atom set, with the
  recombination allocation B&B itself **bounded** (anytime node cap). It
  is *not* the matmul tensor rank. See
  `references/CLOSURE_OPTIMALITY_AND_PERF.md`.
- **Allocation / mask / orbit searches**: anytime / budgeted → bounds
  unless the run reports `exhaustive = true`.

**Rules:**

1. Never write "optimal", "minimal", or "the minimum" for a heuristic
   result. Use "best found", "upper bound", "≤", or name the scope.
2. Carry the tier in the data, not just in prose: prefer an explicit
   field (e.g. `optimality: "bound" | "proven" | "optimal-within:<scope>"`)
   over an unqualified number, so the SPA/paper can render it honestly.
3. When a search has an exhaustiveness flag (e.g. `Result.exhaustive()`),
   propagate it to the label — a budgeted run that *happened* to finish
   the whole space is proven-optimal; the same code hitting its cap is a
   bound.
4. Lower bounds are first-class: when we can bracket an optimum
   (form-count / rank LB vs heuristic UB), report **both** so the gap is
   visible and the claim is certified, not asserted.

This is the algorithmic sibling of the **field discipline** and the
**discoveries-vs-rediscoveries** rule: same spirit — say exactly what
we know, and exactly how strongly we know it.

## PDFs and paper sources: ask the user, don't auto-download

When a task needs a paper PDF (Pan 1980, DIS09, Schwartz-Zwecher 2025, etc.),
the **default is to ask the user to drop the PDF manually** (typically
`/Users/{user}/Downloads/...`) and tell you the path. The user often has
the file already, or can grab it from a paywalled source faster than the
agent can search.

**Exception**: when the user explicitly says "I can't find this paper,
please look for it" — then it's OK for the agent to try web/Wayback/
arXiv/Google Scholar fetches, including:
- WebFetch on canonical URLs (arXiv.org, paper home pages)
- Wayback Machine (`https://web.archive.org/web/.../...`) for dead links
- Google Scholar search for alternative hosts
- The author's personal page or institutional repository

For supplemental data ZIPs (factor matrix archives), Maple/Macaulay2
scripts, etc., the same rule applies — ask first, search only when the
user says they can't find it.

Don't waste cycles speculatively scraping when the user can paste a
local path. Don't assume a paper isn't behind a paywall without trying.
Don't commit downloaded PDFs to the repo unless the user explicitly OKs
it (most are copyrighted).

## Don't bypass safety in git operations

Never use `--no-verify` or skip GPG signing unless the user explicitly
asks. If a hook fails, investigate.
