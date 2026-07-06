---
name: fmm-react
description: React to an FMM-Lille catalog gap — refresh the FMM cross-check, pick a random shape where our catalog is worse than FMM (or take an explicit ⟨n,m,p⟩), analyze the FMM scheme, and try to close the gap via projection, targeted sweep, or serendipitous bud product. Use when asked to "react to FMM", close FMM gaps, or work a shape from the FMM diff.
---

# fmm-react — close a gap against the FMM-Lille catalog

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

Math/structure first (CLAUDE.md "Math first"). Classify the gap:

- **Closure gap** (typical for max-dim > 12): FMM's rank is reachable
  from schemes we already hold — our on-disk best is just a stale
  projection/recombination closure. The cross-check report's
  *Interpretation* paragraph tells you the current split. → react with
  **projection scatter / closure sweep** (§4a, §4b); no import needed.
- **Import gap** (small shapes, or MISSING): FMM holds a genuinely
  better base scheme. → **import it** (§3), then react to its
  *structure* (§4c) so the win propagates to composite shapes.
- **Family gap**: shape belongs to a known parametric family — e.g.
  ⟨N−1,N,N⟩ next to a trilinear-aggregation cube → structured
  projection probe (§4a); ⟨2,3,k⟩ / Hopcroft–Kerr family → usually
  formula-derivable, check `*Bound.java` helpers before searching.

Useful context: `references/fmm-lineage-review.md`,
`references/PROJECTION_MARGIN_TRADEOFF.md`,
`references/BUD_STRUCTURE_THEORY.md`, and whether nearby divisor shapes
of the target are already SOTA (Kron/concat may close it cheaply).

## 3. Import the FMM scheme (import gaps only)

1. Download + parse + write JSON (canonical batch importer; caches in
   `/tmp/fmm-maple-cache`, falls back to the Perminov mirror):

   ```bash
   python3 tools/import_fmm_maple.py NxMxP
   ```

   Writes dronperminov-format JSON under
   `src/main/resources/schemes/section{max_dim}/`.

2. **Verify** — exact `Verifier.isExactNonCubic` for small shapes,
   `Verifier.residualSampled` for big ones (pattern:
   `TestFmmMapleImport`). Never register an unverified scheme.

3. **Stamp metadata** (all dry-run by default):

   ```bash
   mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.migrate.StampFields    -Dexec.args=--execute
   mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.migrate.StampSource    -Dexec.args=--execute
   mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.migrate.StampAdditions -Dexec.args=--execute
   mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.migrate.NarrowFields   -Dexec.args=--apply
   ```

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
3. Persist wins: `eu.solven.matmul.docs.migrate.MaterialiseSerendipitousWins`.

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

- Reacting to a closure gap by importing (wasted) or to an import gap
  by sweeping without the base (hopeless) — diagnose first (§2).
- Citing FMM as discoverer of an imported rank.
- Comparing across fields/commutativity: the cross-check is
  Q-non-commutative; don't "close" it with an F2 or commutative scheme.
- Calling a sweep result "optimal" — it's a bound, or
  optimal-within-scope at best.
- Leaving a sweep CPU-pegged after the question changed.
