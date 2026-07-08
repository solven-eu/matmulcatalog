---
name: fmm-gap
description: React to an FMM-Lille catalog gap — refresh the FMM cross-check, pick a random shape where our catalog is worse than FMM (or take an explicit ⟨n,m,p⟩), analyze the FMM scheme, and try to close the gap via projection, targeted sweep, or serendipitous bud product. Use when asked to "react to FMM", close FMM gaps, or work a shape from the FMM diff.
---

# fmm-gap — close a gap against the FMM-Lille catalog

Recurrent methodology: given a diff against FMM
(https://fmm.univ-lille.fr/), pick a shape where we lag, understand
*why* FMM is better there, and react with the matching engine. All
CLAUDE.md disciplines apply (field discipline, discovery attribution,
optimality tiers, math-first, long-run process hygiene).

## 0. Inputs

- **No argument (default):** run the diff, pick a **random** shape from
  the `WORSE` list (random so repeated invocations spread coverage).
- **Explicit shape** (`NxMxP` or `⟨n,m,p⟩`): skip selection, react to
  that shape.
- Optional modifiers the user may give: `--worst` (take the biggest gap
  instead of random), a strategy hint (`projection` / `sweep` /
  `serendipitous`), or a budget ("30 minutes", "overnight").

## 1. Refresh the diff and pick the target

0. **Compile first** — every command below is `mvn exec:java`, which
   runs whatever sits in `target/classes` WITHOUT recompiling. Stale
   classes silently change driver behaviour (real case: `ProjectFmmGaps`
   ran with a `FieldAwareLookup(R)` from an old build while the source
   said `Q`). One compile up front covers the whole session:

   ```bash
   mvn -q -ntp compile
   ```

1. Refresh the FMM rank digest (fault-tolerant; on network failure it
   keeps the committed digest and exits non-zero — then just reuse the
   committed one and say so):

   ```bash
   mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.verify.SyncReferenceCatalogs -Dexec.args="--fmm"
   ```

2. Regenerate the cross-check (reads the digest + raw on-disk schemes,
   writes `references/fmm-cross-check.md`):

   ```bash
   mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.verify.FmmCrossCheck
   ```

3. Pick a random `WORSE` row (columns: shape | ours | FMM | gap | source):

   ```bash
   awk -F'|' '/^## WORSE/{f=1} /^## /&&!/WORSE/{f=0} f && /^\| ⟨/ {gsub(/[⟨⟩ ]/,"",$2); print $2, $3+0, $4+0, $5+0}' \
       references/fmm-cross-check.md \
     | awk 'BEGIN{srand()} {a[NR]=$0} END{print a[int(rand()*NR)+1]}'
   ```

   Output: `n,m,p ours fmm gap`. Also glance at `## MISSING` — a
   MISSING shape is a valid target too (treat as an import gap).

4. **State the target** to the user before computing: shape, our rank,
   FMM's rank, gap, and FMM's source column if present. Field context is
   Q (the cross-check compares non-commutative characteristic-0 only).

## 2. Diagnose the gap BEFORE running anything heavy

**Always fetch the per-shape FMM page first** — it states the
construction outright, which usually settles the diagnosis in one read:

```bash
curl -sS --compressed "https://fmm.univ-lille.fr/{n}x{m}x{p}.html"
```

The "Algorithm definition" section names the mechanism and its exact
recipe, e.g. `⟨20×24×25:6466⟩ is serendipitous tensor product
(⟨5×12×5:204⟩ − 42) ⊗ ⟨4×2×5:32⟩ + ⟨16×2×5:126⟩ + 2⟨4×6×5:90⟩ +
16⟨4×4×5:61⟩`, or a block decomposition `⟨…⟩ = ⟨12×12×13:1184⟩ + …`.
That tells you (a) which engine to react with, (b) which ingredient
scheme to compare against ours (their base may be bud-richer), and
(c) whether the arithmetic re-prices with today's catalog ranks. If
the site is down (it happens — TLS handshake failures), try the
Wayback Machine (`https://web.archive.org/web/2026id_/<url>`), note
the snapshot date (an old snapshot describes an old rank), and record
a "re-check when the site is back" action in the research log.

Then classify the gap (math/structure first, CLAUDE.md "Math first"):

- **Closure gap** (typical for max-dim > 12): FMM's rank is reachable
  from schemes we already hold — our on-disk best is just a stale
  projection/recombination closure. The cross-check report's
  *Interpretation* paragraph tells you the current split. → react with
  **projection scatter / closure sweep** (§4a, §4b); no import needed.
- **Base gap** (small shapes, or MISSING): FMM holds a genuinely
  better base scheme. We do **NOT** copy it into the catalog — if a
  reaction needs its structure, import it as a **reaction base** (§3)
  and let projection/serendipity propagate the win as our own derived
  schemes.
- **Family gap**: shape belongs to a known parametric family — e.g.
  ⟨N−1,N,N⟩ next to a trilinear-aggregation cube → structured
  projection probe (§4a); ⟨2,3,k⟩ / Hopcroft–Kerr family → usually
  formula-derivable, check `*Bound.java` helpers before searching.

Useful context: `references/fmm-lineage-review.md`,
`references/PROJECTION_MARGIN_TRADEOFF.md`,
`references/BUD_STRUCTURE_THEORY.md`, and whether nearby divisor shapes
of the target are already SOTA (Kron/concat may close it cheaply).

## 3. Import a reaction base — NEVER a direct catalog import

**We never import an FMM scheme directly into the catalog trees**
(`known/`, `curated/`, `constructed/`, `derived/`). FMM's rank stays
visible through the digest + cross-check; bulk upstream import is the
CI Perminov pipeline's job, not this skill's. The only reason to pull
FMM factor matrices here is to serve a reaction, and the JSON goes
into the matching **base folder**:

- **Projection base** (margin-rich; rank-worse than catalog is fine if
  projection-valuable) → `src/main/resources/schemes/margin-bases/sectionN/`
  — policy in `references/PROJECTION_MARGIN_TRADEOFF.md`.
- **Serendipity base** (bud-rich) →
  `src/main/resources/schemes/bud-bases/sectionN/`
  — policy in `references/BUD_STRUCTURE_THEORY.md`.

Retention follows `references/PURGE_REFCOUNT_POLICY.md`: keep the base
only if the reaction actually uses it; if it produced nothing, delete
it before landing.

1. Download + parse + write JSON (canonical batch importer; caches in
   `/tmp/fmm-maple-cache`, falls back to the Perminov mirror):

   ```bash
   python3 tools/import_fmm_maple.py NxMxP
   ```

   Then **move** the emitted JSON out of the default location into the
   base folder chosen above — it must not land in a catalog tree.

2. **Verify** — exact `Verifier.isExactNonCubic` for small shapes,
   `Verifier.residualSampled` for big ones (pattern:
   `TestFmmMapleImport`). Never register an unverified scheme.

3. **Stamp metadata** (all dry-run by default):

   ```bash
   mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.migrate.StampFields    -Dexec.args=--execute
   mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.migrate.StampSource    -Dexec.args=--execute
   mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.migrate.StampAdditions -Dexec.args=--execute
   mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.migrate.NarrowFields   -Dexec.args=--apply
   mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.migrate.BackfillMissingFields -Dexec.args=--execute
   mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.migrate.StampImportHashes     -Dexec.args=--execute
   ```

   `BackfillMissingFields` is what actually stamps `fields[]` on a fresh
   import (StampFields is lineage-driven and skips imports; NarrowFields
   only narrows an existing `fields[]`). `StampImportHashes` stamps the
   content `hash` old-convention import filenames don't carry — without
   it `durableLeafRef` REFUSES to pin the base when a reaction wins
   (the ⟨20,28,28⟩=8434 persist failure, 2026-07-08).

4. **Attribution**: FMM is a catalog, not (usually) the discoverer. Set
   `discovery: false | "TBD"` and `attribution_for_rank` from FMM's
   source column / bibliography (`references/fmm-lille-biblio.json`).
   Scheme JSON is written via `MatrixJsonFormatter.format`, never
   Jackson pretty-print.

## 4. React — pick the engine(s) matching the diagnosis

Run heavy engines **in the background** with progress polling and
heap-dump-on-OOM (CLAUDE.md long-running-procedures section):

```bash
mkdir -p target/oom-dumps
MAVEN_OPTS="-Xmx4g -XX:+HeapDumpOnOutOfMemoryError -XX:HeapDumpPath=target/oom-dumps/"
```

### a) Projection

- **Coordinate-drop scatter over the FMM-gap list** (parses the WORSE
  table itself when `--shapes` is omitted):

  ```bash
  mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.explore.ProjectFmmGaps \
      -Dexec.args="--shapes=NxMxP --passes=2"
  ```

  Persists strict improvements via `RecursiveMaterialiser.projectScatter`.

- **Structured projection** for ⟨N−1,N,N⟩-type shapes next to a TA cube
  (arg2 = `TrilinearAggregations` enum member, e.g. `DIS`, `LITA`):

  ```bash
  mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.explore.ProbeStructuredProjection \
      -Dexec.args="N DIS"
  ```

### b) Targeted sweep / closure

`SchemeSweep` (`--field` mandatory; see class javadoc for the full flag
set). Typical reaction to a single shape:

```bash
mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.SchemeSweep \
    -Dexec.args="--mode=materialize --field=Q --shape=NxMxP --config=thorough --strategies=recomb,serendip,proj"
```

Scope alternatives: `--shape-file=PATH` (feed a curated gap list),
`--band=LO-HI`, `--cubic=LO-HI`. Budget with `--maxNodes` /
`--stagnation`; `--only-if-missing` / `--best-derived` gate writes.
Results are **bounds / optimal-within-scope**, never "optimal".

**Pool starvation caveat**: when §2's recipe names an outer base we
ALREADY pool (dim ≤ 5), do not conclude the engine "can't find it" —
the thorough pool's ~3k entries starve the per-base allocation budget,
so deep/uneven allocations get pruned. Re-run with
`--baseFilter=<file-stem>` (substring on pool label) to concentrate the
B&B on that base family. Precedent: ⟨7,14,24⟩=1514 = ⟨2,4,4⟩:26 with the
two-axis-uneven alloc `[3,4 | 3,4,3,4 | 6,6,6,6]` — invisible under the
full pool, found in seconds under a 6-entry filtered pool. Likewise use
`--base=NxMxP` for bases the `maxBaseDim=5` cap excludes entirely
(⟨3,4,7⟩, ⟨3,4,6⟩, ⟨2,5,6⟩, ⟨2,5,7⟩ — the 2026-07-07 harvest family).

### c) Serendipitous bud product (after importing the base)

The point of importing an FMM base is often its **bud structure**
(groups of rank-1 terms sharing a u/v/w factor — Smith 2002 §9.3):

1. Inspect the imported scheme's buds: `SerendipitousBudProduct.findBuds`
   / `.summarise` (e.g. `"4×U⟨1,1,2⟩ + 12×⟨1,1,1⟩"`). Remember bud
   decomposition is greedy — try `ALL_ORDERINGS` if the default U,V,W
   order looks poor (the ⟨8,9,9⟩=430 bug was exactly a masked V-bud).
2. If bud-rich, sweep products against small second schemes:

   ```bash
   mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.SerendipitousSweep \
       -Dexec.args="--baseCap=8 --secondCap=6 --targetCap=16 --field=Q"
   ```

   Set the caps so base×second covers the target shape (and its
   neighbours — an imported base often improves *other* shapes too).
3. Persist wins: `eu.solven.matmul.docs.migrate.MaterialiseSerendipitousWins`
   (takes `--shapes=NxMxP,…` to land a sweep's win list directly).
4. **Bud-representative caveat**: two same-rank schemes for the same shape
   can carry different bud profiles (direction classes are
   axis-permutation-invariant but NOT flip/GL-invariant). If FMM's recipe
   implies a bud partition that neither our base nor their published base
   shows (σ-arithmetic tells you: a size-k axis-class saves
   `k·R(inner) − R(enlarged)`), the composed representative is a third
   scheme in the rank orbit — react with a flip-graph walk from the base
   under the direct objective `FlipObjectives.serendipitous(Q, n2,m2,p2)`
   (pattern: `docs.explore.ProbeFlipWalkSerendip202425`).
   **BUT check the σ-arithmetic is even expressible first**: if FMM's saving
   requires a k-term class fused into a block NARROWER than k·(inner dim) —
   span compression, e.g. ⟨16,20,28⟩'s `(⟨4,5,7⟩−17) ⊗ ⟨4,4,4⟩ + 8·⟨4,4,8⟩`
   where size-2 buds save zero — then `findBuds`/`costOf` (span-BLIND) can
   neither price nor build it, and the flip walk CANNOT succeed either (its
   objective is the same span-blind cost; 2026-07-08 negative result,
   `ProbeFlipWalkSerendip162028`). That class needs the span-compressed-bud
   engine feature (see the DISCOVERIES log) — log and move on.

### d) GL-orbit frontier probe (closure gaps the others miss)

The recombination pool's `AXIS_FLIP` orbit only does row-reversals; a
full `GL×GL×GL` change-of-basis member of the ⟨2,2,2⟩ base can reach
strictly better allocations (precedent: ⟨17,19,20⟩=3780, −17 vs the
AXIS_FLIP closure — see `research/DISCOVERIES_PENDING_ANALYSIS.md`
2026-06-25). Cheap single-shape probe, worth firing whenever §4a/§4b
come back empty on a closure gap:

```bash
mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.verify.VerifyGLCandidate \
    -Dexec.args="N M P"
```

Band-sweep variant: `docs.explore.GLLargeSweep`. Persist a verified win
with `docs.migrate.RegisterGLWin N M P`.

## 5. Land the results

Whether or not the gap closed:

1. **Verify** every new scheme (`Verifier.passesRandomMatmulSpotCheck`
   at minimum; exact check when feasible) and confirm metadata
   (`fields[]`, `source`, `additions`, `discovery`/`attribution_for_rank`).
2. Optional metrics stamp: `docs.generate.EnrichSchemeMetrics`
   `-Dexec.args="--apply=true"`.
3. Regenerate committed artifacts in the canonical (mmc.sh) order —
   or just `scripts/mmc.sh index`:
   `GenerateDerivedBounds` → `GenerateCitedBounds` → `GenerateCatalogManifest`.
4. Re-run `FmmCrossCheck` and report the before/after for the target
   shape (and any collateral shapes that moved).
5. **Guard**: if the shape was hard to reach (new mechanism, new
   family), add a CSV row to
   `src/test/java/eu/solven/matmul/search/TestSweepSpotsSota.java`
   asserting `≤ SOTA` (never `=`).
6. **Log** the attempt in `research/DISCOVERIES_PENDING_ANALYSIS.md`:
   claim (with field tag, e.g. `Q⟨n,m,p⟩:m`), repro command, breakdown,
   status — **also when nothing improved** (record what was tried and
   with what budget, so the next run doesn't redo it).
7. Kill any background JVMs you started; name anything left running.

## Failure modes to avoid

- Importing an FMM scheme directly into a catalog tree (`known/`,
  `curated/`, `constructed/`, `derived/`) — reaction bases only, in
  `margin-bases/` / `bud-bases/` (§3).
- Reacting to a closure gap by importing a base (wasted) or to a base
  gap by sweeping without the base (hopeless) — diagnose first (§2).
- Citing FMM as discoverer of an imported rank.
- Comparing across fields/commutativity: the cross-check is
  Q-non-commutative; don't "close" it with an F2 or commutative scheme.
- Calling a sweep result "optimal" — it's a bound, or
  optimal-within-scope at best.
- Leaving a sweep CPU-pegged after the question changed.
