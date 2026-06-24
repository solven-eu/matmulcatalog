package eu.solven.matmul.docs;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.atomic.AtomicLongArray;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Shape;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.catalog.CatalogLimits;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.Lineage;
import eu.solven.matmul.recombination.Recombination;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.recombination.BlockSplitSearch;
import eu.solven.matmul.search.ConcatSplitSearch;
import eu.solven.matmul.search.CitedBound;
import eu.solven.matmul.search.PoolConfig;
import eu.solven.matmul.search.RecursiveClosureSota;
import eu.solven.matmul.search.RecursiveMaterialiser;
import eu.solven.matmul.util.ProgressMonitor;
import lombok.extern.slf4j.Slf4j;

/**
 * Unified shape-sweep driver. Replaces {@code PoolConfigSweep},
 * {@code MaterialiseRecursiveSweep} and {@code MaterializeClosureLoop}.
 * Selects behaviour via {@code --mode}:
 *
 * <ul>
 *   <li>{@link Mode#EVALUATE} (default) — multi-config A/B harness. Calls
 *       {@link BlockSplitSearch#findBestStrategy} once per (shape, config),
 *       no writes; renders a markdown comparison report.</li>
 *   <li>{@link Mode#MATERIALIZE} — single config, parallel workers,
 *       {@link RecursiveMaterialiser}. Writes any newly-composed scheme
 *       to {@code src/main/resources/schemes/sectionN/}, recursively
 *       filling in sub-shapes. One pass.</li>
 *   <li>{@link Mode#CLOSURE} — single config, sequential, flat
 *       materialiser via {@link Recombination#constructWithAllocation}.
 *       Iterates rounds until no new file is written.</li>
 * </ul>
 *
 * <p>Full CLI surface. Flags are {@code --key=value}; most accept a couple of
 * spelling aliases (snake / kebab). {@code --field} is the ONLY mandatory flag.</p>
 *
 * <pre>
 * ── core ───────────────────────────────────────────────────────────────────
 *   --mode=evaluate|materialize|closure   default: evaluate
 *   --field=Q|Z|R|C|F2|F3                 MANDATORY — a sweep is field-specific;
 *                                         there is NO default (an unspecified field
 *                                         leaks cross-field phantoms, e.g. an F₂
 *                                         base into a Q search).
 *   --config=NAME[,NAME…]                 pool preset(s): simple | auditAxisFlip |
 *                                         axisFlipOnly | rectangular | includeDerived
 *                                         | thorough. EVALUATE takes many (A/B);
 *                                         materialize/closure take exactly ONE.
 *
 * ── which shapes (scope) ───────────────────────────────────────────────────
 *   --shape=NxMxP[,NxMxP…]                explicit target shape(s)
 *   --shape-file=PATH                     one "NxMxP" per line (leading token; extra
 *                                         text tolerated). The way to dig a curated
 *                                         gap list — e.g. the FMM-gap shapes.
 *   --cubic=LO-HI | N | off               cubic band ⟨d,d,d⟩ (default per mode)
 *   --noncubic=LO-HI | N | off            non-cubic band, all ⟨n,m,p⟩ in range
 *   --band=LO-HI | N                      keep shapes whose max(n,m,p) ∈ [LO,HI]
 *   --minMaxDim=N (band-floor)            keep shapes with max(n,m,p) ≥ N
 *
 * ── search space (pool & allocation enumeration) ───────────────────────────
 *   --base=NxMxP[,NxMxP…]                 use these catalog scheme(s) as the outer
 *                                         recombination base pool (each over its
 *                                         axis-perm orbit) instead of the default
 *                                         rootPool. "Sweep with base ⟨2,4,4⟩."
 *   --maxBaseDim=N                        cap max(bn,bm,bp) of pooled bases — the main
 *                                         pool-cost lever, but COARSE: includeDerived
 *                                         admits bases at dim 5, so dropping to 4
 *                                         yields a Strassen-only 5-base pool (loses
 *                                         ⟨2,4,4⟩ too), while 5 admits both the cheap
 *                                         ⟨2,4,4⟩ AND the expensive ⟨2,5,5⟩ (5-way
 *                                         split → minutes/shape). No single value
 *                                         separates them — for a curated cheap pool
 *                                         WITH ⟨2,4,4⟩, use --base=… instead.
 *   --maxImbalance=N                      cap allocation imbalance (default unbounded)
 *   --balancedOnly=true|false            only equal splits (default false)
 *   --maxCombinations=N                   cap #allocations enumerated per base
 *   --maxPadding=N                        cap zero-padding added to fit a base
 *   --orbitMode=…                         internal-orbit enumeration mode
 *   --commutative[=true]                  also admit commutative-only bases
 *   --baseFilter=… / --baseShape=…        restrict the pool by label / shape
 *
 * ── search budget (anytime caps; bound runtime on large bands) ─────────────
 *   --maxNodes=N    --stagnation=N        flat AllocationOptimizer anytime budget
 *   --assignmentMaxNodes=N                assignment-optimizer node cap
 *   --optimizer=flat|assignment           pick the allocation optimizer
 *
 * ── materialize/closure write-gate & strategies ───────────────────────────
 *   --only-if-missing[=true]              fill-only: write a shape ONLY when missing
 *                                         (never touch an existing entry).
 *                                         DEFAULT is the opposite (improve): rewrite
 *                                         an existing entry when a Verifier-passing
 *                                         composition is strictly better — closing
 *                                         re-derivable gaps is the tool's main job.
 *   --improve=true|false                  deprecated inverse alias of the above.
 *   --best-derived[=true]                 persist the best DERIVED scheme per shape on
 *                                         a TIE-or-better (reproducible derivation).
 *   --derive-all[=true]                   register the derivation even alongside an
 *                                         existing one (implies --best-derived).
 *   --buildable=required|optional         require sub-blocks be buildable (default
 *                                         required — no unbuildable stub churn).
 *   --strategies=recomb,serendip,proj     restrict to a subset of build strategies.
 *   --dry-run=true                        search/report only; write nothing (alias: --skip-materialise).
 *   --projection-only=true                projection strategy only.
 *
 * ── I/O ────────────────────────────────────────────────────────────────────
 *   --out=PATH                            report path (evaluate; default under target/)
 *   --schemes-root=PATH                   write root (materialize/closure). Default
 *                                         src/main/resources/schemes; point at
 *                                         target/staging-schemes to inspect before
 *                                         promoting. Writes land in derived/sectionN/.
 *   --threads=N                           parallelism (materialize only)
 * </pre>
 *
 * <p><b>Worked recipes.</b></p>
 * <pre>
 * # A/B-compare pool presets on the cubic band, no writes:
 *   --mode=evaluate --field=R --config=simple,includeDerived --cubic=2-16
 *
 * # Close a curated FMM-gap list with a CURATED cheap base pool (incl ⟨2,4,4⟩,
 * # excl the 5-way-split ⟨2,5,5⟩ that costs minutes/shape):
 *   --mode=materialize --field=R \
 *       --base=2x2x2,2x2x3,2x3x3,3x3x3,2x2x4,2x4x4,4x4x4 \
 *       --shape-file=target/fmm-gaps-non5.txt \
 *       --schemes-root=src/main/resources/schemes --threads=6
 *
 * # Dig one shape with an explicit base, staged write for inspection:
 *   --mode=materialize --field=R --shape=5x23x32 --base=2x4x4 \
 *       --schemes-root=target/staging-schemes
 * </pre>
 *
 * <p>After a materialize/closure run, regenerate the committed artifacts:
 * {@code GenerateCatalogManifest}, {@code GenerateDerivedBounds},
 * {@code GenerateFmmGapReport}.</p>
 */
@Slf4j
public final class SchemeSweep {

	private static final int MAX_CUBIC_DIM = CatalogLimits.MAX_CUBIC_DIM;
	private static final int MAX_NONCUBIC_DIM = CatalogLimits.MAX_NONCUBIC_DIM;

	/**
	 * Default evaluate-mode report path: under {@code target/} (gitignored).
	 * {@code mvn clean} wipes this — pass {@code --out=PATH} to point
	 * somewhere stable for runs you want to keep.
	 */
	private static final Path DEFAULT_REPORT_OUT = Path.of("target/scheme-sweep/coverage-summary.md");

	/** Default catalog write root for materialize / closure modes. */
	private static final Path DEFAULT_SCHEMES_ROOT = Path.of("src/main/resources/schemes");
	/** Materialize write-gate. DEFAULT = improve: revisit existing entries and
	 *  rewrite when a verified composition is strictly better — closing gaps is the
	 *  tool's primary job. Filling missing-only is the special case, opted into with
	 *  {@code --only-if-missing} (which sets this false). {@code --improve} is the
	 *  deprecated inverse alias. */
	private static volatile boolean MATERIALIZE_IMPROVE = true;
	/** Derive-best mode ({@code --best-derived}): materialise the best DERIVED scheme for
	 *  every shape and persist it on a TIE-or-better with the catalog incumbent — even when
	 *  a hand-crafted/imported equal already exists. The reproducible derivation is a
	 *  first-class per-shape artifact (the FMM-replacement policy). Unbounded search. */
	private static volatile boolean MATERIALIZE_DERIVE_BEST = false;
	/** Derive-all mode ({@code --derive-all}): like best-derived, but persist the
	 *  derivation for EVERY shape even when an imported atom beats it — so the catalog
	 *  carries a reproducible derived witness everywhere, enabling derived-vs-atom
	 *  comparison (e.g. ⟨4,4,4⟩=49 next to the imported 47). Implies best-derived. */
	private static volatile boolean MATERIALIZE_DERIVE_ALL = false;
	/**
	 * Whether a sub-shape's rank must be backed by a buildable construction
	 * ({@code --buildable=required|optional}, default {@code required}).
	 * <ul>
	 *   <li>{@code required} (default): value only buildable ranks — {@code findRank}
	 *       over explicit schemes and replayable stubs, i.e. ranks backed by an
	 *       ACTUAL construction. Constructive: every strategy the search picks can be
	 *       built. Right for writing schemes (materialize) — no unbuildable churn.</li>
	 *   <li>{@code optional}: buildability not required — also value NON-explicit
	 *       bounds (closed-form Pan-TA / Hopcroft-Kerr formulas, cited bounds) via
	 *       {@link CitedBound}. These are upper bounds with no constructive
	 *       scheme behind them, so the result is a RANK claim, not necessarily
	 *       buildable. Use for bound discovery / reporting, or once UB-marked stubs
	 *       are allowed; a writing materialize still rejects leaves it cannot build.</li>
	 * </ul>
	 */
	private static volatile boolean MATERIALIZE_BUILDABLE_REQUIRED = true;

	/** Dependency order for closure: smallest max-axis first, so larger-shape
	 *  splits can reuse discoveries at smaller shapes in the same round. */
	private static final java.util.Comparator<Shape> BY_MAX_AXIS_ASC =
			java.util.Comparator.comparingInt(Shape::maxDim);

	/** Apply {@code --baseFilter} to the outer-base pool (substring match on the
	 *  base label); identity when no filter set. Throws if the filter matches none. */
	private static List<BlockSplitSearch.NamedBase> applyBaseFilter(
			List<BlockSplitSearch.NamedBase> pool, String baseFilter) {
		if (baseFilter == null || baseFilter.isBlank()) return pool;
		List<BlockSplitSearch.NamedBase> filtered = new ArrayList<>();
		for (BlockSplitSearch.NamedBase nb : pool) {
			if (nb.label().contains(baseFilter)) filtered.add(nb);
		}
		System.out.printf("[base-filter] '%s' kept %d of %d pool entries: %s%n",
				baseFilter, filtered.size(), pool.size(),
				filtered.stream().map(BlockSplitSearch.NamedBase::label).toList());
		if (filtered.isEmpty()) {
			throw new IllegalArgumentException(
					"--baseFilter='" + baseFilter + "' matched 0 pool entries");
		}
		return filtered;
	}

	/**
	 * The outer-base pool for a run: the default {@code rootPool}-derived pool
	 * (config-driven), unless {@code --base=NxMxP[,…]} was given — then the pool
	 * is exactly the user's catalog scheme(s), each oriented over its axis-perm
	 * orbit (orientation matters for recombination tiling). {@code --baseFilter}
	 * is applied on top either way.
	 */
	private static List<BlockSplitSearch.NamedBase> buildPoolFor(
			PoolConfig config, RunSpec spec, FieldAwareLookup lookup) {
		List<BlockSplitSearch.NamedBase> pool =
				(spec.baseShapes == null || spec.baseShapes.isEmpty())
						? BlockSplitSearch.buildPool(config, spec.field)
						: userBasePool(spec.baseShapes, lookup);
		pool = applyBaseShapeFilter(pool, spec.baseShape);
		return applyBaseFilter(pool, spec.baseFilter);
	}

	/** Filter the pool by outer-base SHAPE so the ladder's tiers are disjoint:
	 *  {@code cubic} keeps only ⟨n,n,n⟩ bases (= the cubic tier), {@code rectangular}
	 *  keeps only non-cubic bases (the NEW work a rectangular pass adds over a cubic
	 *  pass), {@code all}/null = no filter. Lets T2 run STRICTLY rectangular bases
	 *  rather than re-walking the cubic ones T1 already covered (research/SWEEP_LADDER.md). */
	private static List<BlockSplitSearch.NamedBase> applyBaseShapeFilter(
			List<BlockSplitSearch.NamedBase> pool, String baseShape) {
		if (baseShape == null || baseShape.isBlank() || baseShape.equalsIgnoreCase("all")) return pool;
		boolean wantCubic = baseShape.equalsIgnoreCase("cubic");
		boolean wantRect = baseShape.equalsIgnoreCase("rectangular") || baseShape.equalsIgnoreCase("rect");
		if (!wantCubic && !wantRect) {
			throw new IllegalArgumentException("--baseShape must be cubic|rectangular|all, got '" + baseShape + "'");
		}
		List<BlockSplitSearch.NamedBase> filtered = new ArrayList<>();
		for (BlockSplitSearch.NamedBase nb : pool) {
			var b = nb.base();
			boolean cubic = b.n == b.m && b.m == b.p;
			if (cubic == wantCubic) filtered.add(nb);
		}
		System.out.printf("[base-shape] '%s' kept %d of %d pool entries%n",
				baseShape, filtered.size(), pool.size());
		if (filtered.isEmpty()) {
			throw new IllegalArgumentException("--baseShape='" + baseShape + "' matched 0 pool entries");
		}
		return filtered;
	}

	/** Build a pool from explicit {@code --base} shapes: look up each shape's best
	 *  catalog scheme and add it under every distinct axis orientation. */
	private static List<BlockSplitSearch.NamedBase> userBasePool(
			List<int[]> shapes, FieldAwareLookup lookup) {
		List<BlockSplitSearch.NamedBase> pool = new ArrayList<>();
		java.util.Set<String> seen = new java.util.LinkedHashSet<>();
		for (int[] s : shapes) {
			int n = s[0], m = s[1], p = s[2];
			// A width-1-axis base (e.g. ⟨1,2,2⟩, the symmetric-peel carrier) is the
			// elementary naïve scheme — it has no catalog file, so synthesise it
			// directly rather than dropping it (which would empty the pool).
			Optional<NonCubicBilinearAlgorithm> found = (n == 1 || m == 1 || p == 1)
					? Optional.of(NonCubicBilinearAlgorithm.naive(n, m, p))
					: lookup.find(n, m, p);
			if (found.isEmpty()) {
				log.warn("--base ⟨{},{},{}⟩ not found in catalog for this field — skipped", n, m, p);
				continue;
			}
			NonCubicBilinearAlgorithm base = found.get();
			for (int[] o : new int[][] { {n, m, p}, {n, p, m}, {m, n, p}, {m, p, n}, {p, n, m}, {p, m, n} }) {
				String key = o[0] + "x" + o[1] + "x" + o[2];
				if (!seen.add(key)) continue;
				base.orientAs(o[0], o[1], o[2]).ifPresent(or ->
						pool.add(new BlockSplitSearch.NamedBase("base<" + key + ">=" + or.r, or)));
			}
		}
		if (pool.isEmpty()) {
			throw new IllegalArgumentException("--base produced an empty pool (no shapes resolved in the catalog)");
		}
		log.info("--base pool: {} oriented base(s) from {} shape(s): {}",
				pool.size(), shapes.size(), pool.stream().map(BlockSplitSearch.NamedBase::label).toList());
		return pool;
	}

	private static final List<NamedConfig> ALL_CONFIGS = List.of(
			new NamedConfig("simple", PoolConfig.simple()),
			new NamedConfig("auditAxisFlip", PoolConfig.auditAxisFlip()),
			new NamedConfig("axisFlipOnly", PoolConfig.axisFlipOnly()),
			new NamedConfig("rectangular", PoolConfig.rectangular()),
			new NamedConfig("includeDerived", PoolConfig.includeDerived()),
			new NamedConfig("thorough", PoolConfig.thorough())
	// auditPermutation omitted: at maxBaseDim=5 it collapses to AXIS_FLIP via the cap.
	);

	private static final List<NamedConfig> DEFAULT_CONFIGS = List.of(
			new NamedConfig("simple", PoolConfig.simple()));

	private static final int OFF = -1;

	public enum Mode { EVALUATE, MATERIALIZE, CLOSURE }

	private SchemeSweep() {}

	private record NamedConfig(String name, PoolConfig config) {}

	private record RunSpec(
			Mode mode,
			String field,
			List<String> configNames,
			int cubicMin, int cubicMax,
			int noncubicMin, int noncubicMax,
			List<int[]> extraShapes,
			Path out,
			Path schemesRoot,
			int threads,
			Boolean balancedOnlyOverride,
			Integer maxImbalanceOverride,
			Integer maxBaseDimOverride,
			eu.solven.matmul.SymmetryTransforms.InternalOrbitMode orbitModeOverride,
			Integer maxCombinationsOverride,
			Integer maxPaddingOverride,
			boolean skipMaterialise,
			boolean projectionOnly,
			String baseFilter,
			String baseShape,
			List<int[]> baseShapes,
			int minMaxDim,
			int bandMin, int bandMax,
			java.util.Set<String> strategies,
			boolean allowCommutative) {}

	public static void main(String[] args) throws Exception {
		RunSpec spec = parseArgs(args);
		List<NamedConfig> configs = resolveConfigs(spec.configNames);
		// Per-run CLI overrides for the search-shape knobs. If the user passes
		// --balancedOnly or --maxImbalance, override every selected config's
		// preset values. Lets a user A/B "simple but unbalanced" without
		// inventing a new factory.
		if (spec.balancedOnlyOverride != null || spec.maxImbalanceOverride != null
				|| spec.maxBaseDimOverride != null || spec.orbitModeOverride != null
				|| spec.maxCombinationsOverride != null || spec.maxPaddingOverride != null) {
			List<NamedConfig> overridden = new ArrayList<>();
			for (NamedConfig nc : configs) {
				PoolConfig pc = nc.config;
				if (spec.balancedOnlyOverride != null) pc = pc.withBalancedOnly(spec.balancedOnlyOverride);
				if (spec.maxImbalanceOverride != null) pc = pc.withMaxImbalance(spec.maxImbalanceOverride);
				if (spec.maxBaseDimOverride != null) pc = pc.withMaxBaseDim(spec.maxBaseDimOverride);
				if (spec.orbitModeOverride != null) pc = pc.withOrbitMode(spec.orbitModeOverride);
				if (spec.maxCombinationsOverride != null) pc = pc.withMaxCombinations(spec.maxCombinationsOverride);
				if (spec.maxPaddingOverride != null) pc = pc.withMaxPadding(spec.maxPaddingOverride);
				overridden.add(new NamedConfig(nc.name, pc));
			}
			configs = overridden;
		}
		List<Shape> shapes = buildShapeSet(spec);
		if (shapes.isEmpty()) {
			throw new IllegalArgumentException("No shapes selected — check --shape / --cubic / --noncubic flags.");
		}
		log.info("Mode={} field={} — {} shapes × {} config(s) ({})",
				spec.mode, spec.field, shapes.size(), configs.size(),
				configs.stream().map(NamedConfig::name).toList());

		switch (spec.mode) {
			case EVALUATE -> runEvaluate(spec, configs, shapes);
			case MATERIALIZE -> runMaterialize(spec, requireSingle(configs, spec.mode), shapes);
			case CLOSURE -> runClosure(spec, requireSingle(configs, spec.mode), shapes);
		}
	}

	private static NamedConfig requireSingle(List<NamedConfig> configs, Mode mode) {
		if (configs.size() != 1) {
			throw new IllegalArgumentException("Mode " + mode
					+ " requires exactly one --config (got " + configs.size() + ").");
		}
		return configs.get(0);
	}

	// ---------------------------------------------------------------------
	// EVALUATE mode (was PoolConfigSweep)
	// ---------------------------------------------------------------------

	/**
	 * ETA weight for one ⟨n,m,p⟩ shape (#117). Per-shape work is dominated by
	 * the recombination search + verification, whose cost grows roughly as
	 * {@code O((n·m·p)²·r)} — super-linear in volume, not linear. A flat-volume
	 * weight (the old {@code n·m·p}) under-estimates the tail, where the last
	 * batches process matrices 10–30× larger than the first; a max-dim³ weight
	 * is worse still for UNBALANCED shapes (it scores ⟨2,2,16⟩ and ⟨16,16,16⟩
	 * equally despite a 64× volume gap). We weight by {@code (n·m·p)²} so the
	 * progress ETA tracks the actual super-quadratic cost across both balanced
	 * and unbalanced sweeps. (The extra {@code ·r} factor is omitted: r itself
	 * scales with volume, and keeping the weight rank-free avoids a lookup per
	 * remaining item.)
	 */
	private static long etaWeight(Shape s) {
		long vol = (long) s.n() * s.m() * s.p();
		return vol * vol;
	}

	private static void runEvaluate(RunSpec spec, List<NamedConfig> configs, List<Shape> shapes) throws IOException {
		FieldAwareLookup lookup = new FieldAwareLookup(spec.field);
		CitedBound sota = new CitedBound(lookup);

		log.info("Building pools for each config …");
		Map<String, List<BlockSplitSearch.NamedBase>> poolsByConfig = new LinkedHashMap<>();
		for (NamedConfig nc : configs) {
			long t0 = System.nanoTime();
			List<BlockSplitSearch.NamedBase> pool = buildPoolFor(nc.config, spec, lookup);
			poolsByConfig.put(nc.name, pool);
			log.info("  {}: {} entries  ({} ms)",
					nc.name, pool.size(), (System.nanoTime() - t0) / 1_000_000L);
		}

		Map<Shape, Map<String, Long>> results = new LinkedHashMap<>();
		for (Shape s : shapes) results.put(s, new LinkedHashMap<>());
		// Best PURE composition rank per shape (before folding in the catalog),
		// so we can report better/tie/worse vs the on-disk catalog SOTA.
		Map<Shape, Long> bestComp = new LinkedHashMap<>();

		ProgressMonitor monitor = ProgressMonitor.builder()
				.total(shapes.size())
				.weight(i -> etaWeight(shapes.get((int) i)))
				.label("SchemeSweep:evaluate")
				.logEveryMillis(5000)
				.logger(log)
				.build();
		monitor.start();
		for (int i = 0; i < shapes.size(); i++) {
			Shape s = shapes.get(i);
			for (NamedConfig nc : configs) {
				List<BlockSplitSearch.NamedBase> pool = poolsByConfig.get(nc.name);
				Optional<BlockSplitSearch.NonCubicStrategy> best =
						BlockSplitSearch.findBestStrategy(s.n(), s.m(), s.p(), pool, sota,
								nc.config.balancedOnly(), nc.config.maxImbalance(),
								nc.config.maxCombinations(), nc.config.maxPadding());
				long compRank = best.map(BlockSplitSearch.NonCubicStrategy::rank).orElse(-1L);
				if (compRank > 0) bestComp.merge(s, compRank, Math::min);
				long rank = compRank;
				long sotaDirect = sota.getRank(s.n(), s.m(), s.p());
				if (sotaDirect > 0 && sotaDirect < (long) s.n() * s.m() * s.p()
						&& (rank < 0 || sotaDirect < rank)) {
					rank = sotaDirect;
				}
				results.get(s).put(nc.name, rank);
			}
			monitor.tick(i);
		}
		monitor.done();

		String report = renderReport(shapes, results, configs)
				+ catalogComparison(shapes, bestComp, lookup);
		if (spec.out.getParent() != null) Files.createDirectories(spec.out.getParent());
		Files.writeString(spec.out, report);  // truncates + overwrites every run
		log.info("Wrote {} (overwritten, not appended)", spec.out);
		System.out.println();
		System.out.println("================ SchemeSweep evaluate report ================");
		System.out.println(report);
		System.out.println("================ end of report ================");
	}

	/**
	 * Classify each shape's best PURE composition (no catalog fold) against the
	 * on-disk catalog SOTA: better / tie / worse. "better" means the composition
	 * strictly beats the catalog (a discovery to investigate, or a bug); "tie"
	 * means the engine reproduces the catalog optimum by composition; "worse"
	 * means the catalog has an imported atom no composition reaches.
	 */
	private static String catalogComparison(List<Shape> shapes,
			Map<Shape, Long> bestComp, FieldAwareLookup lookup) {
		int better = 0, tie = 0, worse = 0, noCatalog = 0, noComp = 0;
		List<String> betterRows = new ArrayList<>();
		// worse rows kept as (gap, formatted-row) so we can sort by gap desc.
		List<long[]> worseGaps = new ArrayList<>();
		List<String> worseRows = new ArrayList<>();
		for (Shape s : shapes) {
			Long comp = bestComp.get(s);
			if (comp == null) { noComp++; continue; }
			int cat = lookup.findRank(s.n(), s.m(), s.p());
			boolean haveCat = cat < Recombination.SotaResolver.UNKNOWN_RANK;
			if (!haveCat) { noCatalog++; continue; }
			if (comp < cat) {
				better++;
				betterRows.add(String.format("| ⟨%d,%d,%d⟩ | %d | %d | %d |",
						s.n(), s.m(), s.p(), comp, cat, cat - comp));
			} else if (comp == cat) {
				tie++;
			} else {
				worse++;
				worseGaps.add(new long[] { comp - cat, s.n(), s.m(), s.p(), comp, cat });
			}
		}
		// Biggest composition-vs-catalog gap first (the most "imported" atoms).
		worseGaps.sort((a, b) -> Long.compare(b[0], a[0]));
		for (long[] g : worseGaps) {
			worseRows.add(String.format("| ⟨%d,%d,%d⟩ | %d | %d | +%d |",
					g[1], g[2], g[3], g[4], g[5], g[0]));
		}
		StringBuilder sb = new StringBuilder();
		sb.append("\n\n## Composition vs catalog (better / tie / worse)\n\n");
		sb.append("Best **pure composition** rank per shape (Kronecker / concat / "
				+ "single-base recombination) vs the on-disk catalog SOTA. "
				+ "_better_ = composition strictly beats catalog (investigate); "
				+ "_tie_ = engine reproduces the catalog optimum; "
				+ "_worse_ = catalog has an atom no composition reaches.\n\n");
		sb.append(String.format("- **better: %d**%n- tie: %d%n- worse: %d%n"
				+ "- no-catalog: %d%n- no-composition: %d%n",
				better, tie, worse, noCatalog, noComp));
		if (!betterRows.isEmpty()) {
			sb.append("\n### BETTER than catalog (investigate)\n\n");
			sb.append("| shape | composition | catalog | Δ |\n| --- | ---: | ---: | ---: |\n");
			betterRows.forEach(r -> sb.append(r).append('\n'));
		}
		if (!worseRows.isEmpty()) {
			sb.append("\n### WORSE than catalog (catalog atom no single-base composition reaches)\n\n");
			sb.append("Sorted by gap (composition − catalog) descending; the gap is the "
					+ "rank a single-base, non-recursive split leaves on the table versus the "
					+ "imported atom (HK / Smirnov / Pan TA / FMM / AlphaTensor-family / …).\n\n");
			sb.append("| shape | composition | catalog | gap |\n| --- | ---: | ---: | ---: |\n");
			worseRows.forEach(r -> sb.append(r).append('\n'));
		}
		return sb.toString();
	}

	// ---------------------------------------------------------------------
	// MATERIALIZE mode (was MaterialiseRecursiveSweep)
	// ---------------------------------------------------------------------

	private static void runMaterialize(RunSpec spec, NamedConfig config, List<Shape> shapes) {
		FieldAwareLookup lookup = new FieldAwareLookup(spec.field);
		List<BlockSplitSearch.NamedBase> pool = buildPoolFor(config.config, spec, lookup);

		int nThreads = spec.threads > 0 ? spec.threads
				: Math.max(1, Runtime.getRuntime().availableProcessors() - 1);
		// DEFAULT = false (maximum unbalance): silently pruning unbalanced splits is a
		// footgun (feedback_dont_silently_prune_search_space) — e.g. ⟨25,26,32⟩ missed the
		// (6,26) p-concat (10954) and stored 11343; ⟨5,32,32⟩ missed the ⟨2,4,4⟩ recombination
		// allocA=[3,2] route (3320) and stored 3446. Balanced-only is now an OPT-IN
		// (--balancedOnly=true), consistent with PoolConfig's own back-compat default.
		// A user who wants the cheaper balanced search can still ask for it.
		boolean composeBalancedOnly = spec.balancedOnlyOverride != null
				? spec.balancedOnlyOverride : false;
		ConcurrentLinkedQueue<RecursiveMaterialiser> mats = new ConcurrentLinkedQueue<>();
		for (int t = 0; t < nThreads; t++) {
			// FLAT + CONSTRUCTIVE resolver — materialize is a SINGLE PASS valuing
			// sub-shapes at their CURRENT catalog rank via findRank (live: reflects
			// schemes written earlier this sweep; excludes commutative; points only at
			// real/replayable entries). NOT RecursiveClosureSota (recursive → multi-hour),
			// and NOT CitedBound — its Pan-TA/Hopcroft-Kerr/Waksman FORMULA bounds
			// have no constructive scheme behind them, so the search picks unbuildable
			// strategies that then fail the spot-check / replay (the 1819-failure corruption,
			// 2026-06-13). Recursion/fixed-point belongs to Mode.CLOSURE, which iterates
			// these flat findRank passes — exactly the resolver it already uses (line ~646).
			// findRank now returns the NAÏVE rank a·b·c (always buildable) for absent shapes —
			// never the MAX/100 sentinel that used to poison AllocationOptimizer's root lower
			// bound and drop good bases (the ⟨12,13,13⟩=1274 via ⟨2,4,4⟩ loss). So a raw findRank
			// is a safe, constructive sota here.
			Recombination.SotaResolver sota = MATERIALIZE_BUILDABLE_REQUIRED
					? (a, b, c) -> lookup.findRank(a, b, c)         // --buildable=required (default): constructive
					: new CitedBound(lookup);                // --buildable=optional: incl. non-explicit bounds
			RecursiveMaterialiser m = new RecursiveMaterialiser(
					lookup, pool, sota, spec.schemesRoot, true, composeBalancedOnly,
					MATERIALIZE_IMPROVE, MATERIALIZE_DERIVE_BEST);
			m.registerDerivedAnyway(MATERIALIZE_DERIVE_ALL);
			m.setStrategies(spec.strategies);
			m.setAllowCommutative(spec.allowCommutative());
			// Opt-in allocation-enumeration caps (default unbounded). The lever that
			// makes a ⟨2,5,5⟩-class base tractable under the rich pool — see
			// RecursiveMaterialiser.setMaxCombinations. Silent default stays unbounded.
			if (spec.maxImbalanceOverride != null) m.setMaxImbalance(spec.maxImbalanceOverride);
			if (spec.maxCombinationsOverride != null) m.setMaxCombinations(spec.maxCombinationsOverride);
			mats.add(m);
		}
		ExecutorService exec = Executors.newFixedThreadPool(nThreads);
		log.info("Parallel sweep on {} threads, writing to {}", nThreads, spec.schemesRoot);

		// results layout: [wins, ties, skipped, noImprovement, totalProcessed,
		// startNanos, errors]. "noImprovement" = compose found nothing strictly
		// better than the committed scheme (the common, expected outcome) — NOT a
		// failure; genuine exceptions go in "errors".
		AtomicLongArray results = new AtomicLongArray(7);
		results.set(5, System.nanoTime());
		System.out.printf("%-12s  %-8s  %-8s  %s%n", "shape", "predict", "catalog", "outcome");
		System.out.println("-".repeat(60));

		List<Future<?>> futures = new ArrayList<>();
		for (Shape sk : shapes) {
			final int n = sk.n(), m = sk.m(), p = sk.p();
			futures.add(exec.submit(() -> {
				RecursiveMaterialiser worker = mats.poll();
				if (worker == null) worker = mats.peek();
				try {
					materializeOne(worker, lookup, n, m, p, results);
				} finally {
					if (worker != null) mats.offer(worker);
				}
			}));
		}
		// FAIL LOUD, but RESILIENT and LEAK-FREE. Two bugs lived here before:
		//   (1) throwing on the FIRST failed future aborted the join loop, so a single
		//       per-shape predict/build divergence (assertRebuildNotWorse) HID every OTHER
		//       divergence — you had to rerun N times to find N bugs. We now COLLECT all
		//       failures and fail ONCE at the end, so a single sweep surfaces every phantom
		//       over-claim (e.g. the 5 band-17 divergences ⟨4,16,17⟩…⟨13,13,17⟩).
		//   (2) that same early throw skipped exec.shutdown(), leaving the non-daemon fixed
		//       thread pool alive → the JVM never exited (mvn joinNonDaemonThreads hung
		//       forever at 0% CPU, work already written). The shutdown now runs in finally.
		List<Throwable> failures = new ArrayList<>();
		try {
			for (Future<?> f : futures) {
				try {
					f.get();
				} catch (Exception e) {
					log.error("sweep worker future failed", e);
					failures.add(e);
				}
			}
		} finally {
			exec.shutdownNow();
		}

		System.out.println("-".repeat(60));
		System.out.printf("Summary: %d wins, %d ties, %d skipped (direct), %d no-improvement, %d errors%n",
				results.get(0), results.get(1), results.get(2), results.get(3), results.get(6));

		if (!failures.isEmpty()) {
			throw new IllegalStateException(failures.size()
					+ " shape(s) failed during materialize (see ERROR logs above) — failing loud "
					+ "AFTER completing the full sweep so all divergences surface in one run",
					failures.get(0));
		}
	}

	private static void materializeOne(RecursiveMaterialiser mat, FieldAwareLookup lookup,
			int n, int m, int p, AtomicLongArray results) {
		String shape = String.format("⟨%d,%d,%d⟩", n, m, p);
		int catalogRank = lookup.find(n, m, p).map(a -> a.r).orElse(-1);
		try {
			Optional<RecursiveMaterialiser.Result> r = mat.materialise(n, m, p);
			if (r.isEmpty()) {
				results.incrementAndGet(3);
			} else {
				NonCubicBilinearAlgorithm alg = r.get().alg();
				if (r.get().fromDisk()) {
					results.incrementAndGet(2);
				} else if (catalogRank < 0 || alg.r < catalogRank) {
					results.incrementAndGet(0);
					System.out.printf("%-12s  %-8d  %-8s  WIN%n", shape, alg.r,
							catalogRank < 0 ? "—" : Integer.toString(catalogRank));
				} else {
					results.incrementAndGet(1);
				}
			}
		} catch (Throwable e) {
			// FAIL LOUD (was: `catch (Exception e)` with NO log, swallowing the error).
			// Widened to Throwable: an Error (OutOfMemoryError / StackOverflowError from
			// expanding a huge projection-parent stub to double[][]) previously escaped the
			// catch(Exception) here AND the future barrier, leaving every counter at 0 and
			// silently discarding an already-found recombination win. Surface the stack.
			results.incrementAndGet(6);
			log.error("materialise {} threw — failing loud (was silently swallowed)", shape, e);
			throw new IllegalStateException("materialise " + shape + " failed", e);
		}
		long processed = results.incrementAndGet(4);
		if (processed % 100 == 0) {
			long elapsedMs = (System.nanoTime() - results.get(5)) / 1_000_000L;
			System.out.printf("  [progress] %4d processed  (%d wins, %d ties, %d direct, %d no-better, %d err)  %dms%n",
					processed, results.get(0), results.get(1), results.get(2), results.get(3), results.get(6), elapsedMs);
		}
	}

	// ---------------------------------------------------------------------
	// CLOSURE mode (was MaterializeClosureLoop)
	// ---------------------------------------------------------------------

	private static void runClosure(RunSpec spec, NamedConfig config, List<Shape> shapes) throws Exception {
		long runStartNs = System.nanoTime();
		FieldAwareLookup lookup = new FieldAwareLookup(spec.field);
		List<BlockSplitSearch.NamedBase> pool = buildPoolFor(config.config, spec, lookup);

		// Dependency order: smaller max-axis first, so larger-shape splits
		// see discoveries at smaller shapes during the same round.
		List<Shape> ordered = new ArrayList<>(shapes);
		ordered.sort(BY_MAX_AXIS_ASC);

		// Per-shape phase logging threshold — single-shape probes need
		// visibility into search/materialise timing.
		boolean verbosePerShape = ordered.size() <= 4;

		// --projection-only: skip the upward search + materialise (Phases 1/2)
		// and run just the downward projection closure (Phase 3). Cheap way to
		// probe / apply the projection operator against an already-converged
		// catalog without paying for the full recombination search.
		if (spec.projectionOnly) {
			System.out.println("Projection-only run (--projection-only=true): "
					+ "skipping search + materialise phases.");
			runProjectionClosure(ordered, pool, lookup, spec);
			System.out.printf("%n=== runClosure complete: %s wall-clock ===%n",
					formatDuration(System.nanoTime() - runStartNs));
			return;
		}

		// The upward search phase is the recombination/Kron/concat operator
		// (findBestStrategy). Skip it entirely when 'recombination' is not selected
		// — e.g. a serendipity-only closure run. (Serendipity itself fires per-shape
		// in compose(); in closure mode it only sees overlay shapes, so a pure
		// serendipity sweep is better run via --mode=materialize.)
		Map<Shape, BlockSplitSearch.NonCubicStrategy> pendingOverlay =
				spec.strategies.contains(RecursiveMaterialiser.STRAT_RECOMBINATION)
						? runSearchRounds(ordered, pool, lookup, config, verbosePerShape)
						: new java.util.LinkedHashMap<>();
		if (pendingOverlay.isEmpty()
				&& !spec.strategies.contains(RecursiveMaterialiser.STRAT_RECOMBINATION)) {
			System.out.println();
			System.out.println("Skipping upward search (--strategies excludes 'recombination').");
		}

		// ── PHASE 2: materialise overlay wins (unless --dry-run). ─
		if (spec.skipMaterialise) {
			System.out.println();
			System.out.println("Skipping materialisation (--dry-run=true).");
			if (!pendingOverlay.isEmpty()) {
				System.out.println("Pending overlay (not written to disk):");
				List<Shape> orderedOverlay = new ArrayList<>(pendingOverlay.keySet());
				orderedOverlay.sort(BY_MAX_AXIS_ASC);
				for (Shape k : orderedOverlay) {
					BlockSplitSearch.NonCubicStrategy w = pendingOverlay.get(k);
					System.out.printf("  %s = %d  via %s%n", k, w.rank(), w.label());
				}
			}
		} else if (!pendingOverlay.isEmpty()) {
			materialiseOverlay(pendingOverlay, pool, lookup, spec, verbosePerShape);
		}

		// ── PHASE 3: downward projection closure (multipass). ──────────────
		// Unlike the upward operators (Kron / concat / recombination /
		// serendipitous), projection improves a SMALLER shape from a better
		// LARGER one — so it breaks the smallest-first acyclicity and must
		// iterate to a fixpoint. See paper §multipass.
		if (!spec.skipMaterialise
				&& spec.strategies.contains(RecursiveMaterialiser.STRAT_PROJECTION)) {
			runProjectionClosure(ordered, pool, lookup, spec);
		} else if (!spec.strategies.contains(RecursiveMaterialiser.STRAT_PROJECTION)) {
			System.out.println();
			System.out.println("Skipping projection closure (--strategies excludes 'projection').");
		}

		System.out.printf("%n=== runClosure complete: %s wall-clock ===%n",
				formatDuration(System.nanoTime() - runStartNs));
	}

	/**
	 * PHASE 1 of {@link #runClosure}: iterate search rounds against an in-memory
	 * overlay (per-round wins are visible to later rounds without disk I/O) until
	 * a round produces no new strict win. Returns the converged overlay.
	 */
	private static Map<Shape, BlockSplitSearch.NonCubicStrategy> runSearchRounds(
			List<Shape> ordered, List<BlockSplitSearch.NamedBase> pool,
			FieldAwareLookup lookup, NamedConfig config, boolean verbosePerShape) {
		// Per-round wins live in {@code pendingOverlay}; the resolver consults the
		// overlay first so subsequent rounds see prior-round discoveries WITHOUT
		// having to materialise + write + reload from disk.
		Map<Shape, BlockSplitSearch.NonCubicStrategy> pendingOverlay = new HashMap<>();

		Recombination.SotaResolver flat = (a, b, c) -> {
			BlockSplitSearch.NonCubicStrategy ov = pendingOverlay.get(new Shape(a, b, c));
			// Rank-only fast path: skips CompactScheme.expand() and the
			// double[][] U/V/W materialisation that find() would otherwise
			// do per call. See FieldAwareLookup.findRank.
			int fromCatalog = lookup.findRank(a, b, c);
			if (ov == null) return fromCatalog;
			return Math.min(fromCatalog, (int) ov.rank());
		};

		int round = 0;
		while (true) {
			round++;
			List<String> winsThisRound = new ArrayList<>();
			List<String> attemptDigest = new ArrayList<>();
			long bestDeltaThisRound = 0L;

			ProgressMonitor monitor = ProgressMonitor.builder()
					.total(ordered.size())
					.weight(i -> etaWeight(ordered.get((int) i)))
					.logEveryMillis(5_000)
					.label("SchemeSweep:search round " + round)
					.build();
			monitor.start();

			for (int idx = 0; idx < ordered.size(); idx++) {
				Shape s = ordered.get(idx);
				int n = s.n(), m = s.m(), p = s.p();
				try {
					Optional<eu.solven.matmul.catalog.FieldAwareLookup.WithSource> currentWithSource =
							lookup.findWithSource(n, m, p);
					Optional<Integer> catalogRank = currentWithSource.map(ws -> ws.alg().r);
					String currentFile = currentWithSource
							.map(ws -> ws.path() == null ? "?" : ws.path().getFileName().toString())
							.orElse("—");
					// Effective current = min(catalog, overlay)
					BlockSplitSearch.NonCubicStrategy overlayWin = pendingOverlay.get(s);
					Integer currentRank = catalogRank.orElse(null);
					String currentSource = currentFile;
					if (overlayWin != null && (currentRank == null || overlayWin.rank() < currentRank)) {
						currentRank = (int) overlayWin.rank();
						currentSource = "overlay@round" + round + " via " + overlayWin.label();
					}
					if (verbosePerShape) {
						System.out.printf("  [phase] ⟨%d,%d,%d⟩ search start: pool=%d balancedOnly=%s "
								+ "maxImbalance=%d maxCombos=%d maxPadding=%d current=%s%n",
								n, m, p, pool.size(), config.config.balancedOnly(),
								config.config.maxImbalance(), config.config.maxCombinations(),
								config.config.maxPadding(),
								currentRank == null ? "—" : currentRank.toString());
					}
					long searchStart = System.nanoTime();
					Optional<BlockSplitSearch.NonCubicStrategy> picked =
							BlockSplitSearch.findBestStrategy(n, m, p, pool, flat,
									config.config.balancedOnly(), config.config.maxImbalance(),
									config.config.maxCombinations(), config.config.maxPadding());
					if (verbosePerShape) {
						long searchMs = (System.nanoTime() - searchStart) / 1_000_000L;
						String pickedLabel = picked.map(q -> "rank=" + q.rank() + " via " + q.label())
								.orElse("no-candidate");
						System.out.printf("  [phase] ⟨%d,%d,%d⟩ search done in %d ms — %s%n",
								n, m, p, searchMs, pickedLabel);
					}
					if (picked.isEmpty()) {
						attemptDigest.add(String.format("⟨%d,%d,%d⟩  no-candidate  current=%s (%s)",
								n, m, p, currentRank == null ? "—" : currentRank.toString(), currentSource));
						continue;
					}
					long predicted = picked.get().rank();
					String curStr = currentRank == null ? "—" : currentRank.toString();
					if (currentRank != null) {
						long delta = predicted - currentRank;
						if (delta < bestDeltaThisRound) bestDeltaThisRound = delta;
						if (predicted >= currentRank) {
							String tag = (predicted == currentRank) ? "TIE " : "WORSE";
							attemptDigest.add(String.format(
									"⟨%d,%d,%d⟩  %s   predicted=%d via %s  |  current=%s (%s)  Δ=%+d",
									n, m, p, tag, predicted, picked.get().label(), curStr, currentSource, delta));
							continue;
						}
					}
					// New strict win — record in overlay only, do NOT materialise yet.
					pendingOverlay.put(s, picked.get());
					winsThisRound.add(String.format("⟨%d,%d,%d⟩=%d via %s  (was %s [%s])",
							n, m, p, predicted, picked.get().label(), curStr, currentSource));
				} finally {
					monitor.tick(idx);
				}
			}
			monitor.done();

			System.out.printf("%n=== Search round %d: %d new overlay wins (best Δ attempted: %+d) ===%n",
					round, winsThisRound.size(), bestDeltaThisRound);
			if (!winsThisRound.isEmpty()) {
				Collections.sort(winsThisRound);
				System.out.println("  -- new wins (overlay) --");
				for (String s : winsThisRound) System.out.println("    " + s);
			}
			if (!attemptDigest.isEmpty()) {
				Collections.sort(attemptDigest);
				int cap = Math.min(attemptDigest.size(), 30);
				System.out.println("  -- attempts (no win) --");
				for (int i = 0; i < cap; i++) System.out.println("    " + attemptDigest.get(i));
				if (attemptDigest.size() > cap) {
					System.out.printf("    ... and %d more (capped); rerun with --shape=NxMxP to drill in%n",
							attemptDigest.size() - cap);
				}
			}
			if (winsThisRound.isEmpty()) break;
		}
		System.out.printf("%nSearch converged in %d round(s). Total overlay wins: %d%n",
				round, pendingOverlay.size());
		return pendingOverlay;
	}

	/**
	 * PHASE 2 of {@link #runClosure}: materialise each overlay win in dependency
	 * order (smallest max-axis first), exposing just-materialised smaller winners
	 * as sub-products so larger winners can reuse them. Verifies each scheme and
	 * writes it under {@code spec.schemesRoot}.
	 */
	private static void materialiseOverlay(
			Map<Shape, BlockSplitSearch.NonCubicStrategy> pendingOverlay,
			List<BlockSplitSearch.NamedBase> pool, FieldAwareLookup lookup,
			RunSpec spec, boolean verbosePerShape) throws IOException {
		System.out.println();
		System.out.println("=== PHASE 2: materialise overlay wins (RecursiveMaterialiser) ===");
		// Drive the FULL RecursiveMaterialiser (same component MATERIALIZE mode uses):
		// it handles every composition op (Kronecker / concat / recombination),
		// records the construction lineage, writes lineage-only stubs above
		// MATERIALISE_MAX_DIM, and drops built matrices after persisting — so memory
		// stays bounded (the previous ad-hoc materialiser retained every built scheme
		// in a map → OOM, had no Kronecker support → NPE, and recorded lineage only
		// for recombination wins). It writes each improvement to disk itself.
		RecursiveClosureSota sota = new RecursiveClosureSota(lookup, pool, true, true);
		RecursiveMaterialiser mat = new RecursiveMaterialiser(
				lookup, pool, sota, spec.schemesRoot, true, true, /* improveExisting */ true);
		mat.setStrategies(spec.strategies);
		mat.setAllowCommutative(spec.allowCommutative());
		List<Shape> orderedOverlay = new ArrayList<>(pendingOverlay.keySet());
		orderedOverlay.sort(BY_MAX_AXIS_ASC);
		AtomicLongArray results = new AtomicLongArray(6);
		results.set(5, System.nanoTime());
		for (Shape s : orderedOverlay) {
			if (verbosePerShape) {
				BlockSplitSearch.NonCubicStrategy w = pendingOverlay.get(s);
				System.out.printf("  [phase] ⟨%d,%d,%d⟩ materialise (target rank=%d) via %s%n",
						s.n(), s.m(), s.p(), w.rank(), w.label());
			}
			materializeOne(mat, lookup, s.n(), s.m(), s.p(), results);
		}
		System.out.printf("%nMaterialised overlay: %d wins, %d ties, %d direct, %d failed (of %d).%n",
				results.get(0), results.get(1), results.get(2), results.get(3), pendingOverlay.size());
	}

	/** Safety cap on projection-closure rounds (each round is a full descending sweep). */
	private static final int PROJECTION_MAX_ROUNDS = 8;

	/**
	 * PHASE 3 of {@link #runClosure}: the downward projection closure. For every
	 * shape (largest max-axis first, so a freshly-improved parent is visible when
	 * its child projects in the same round), attempt to improve it by restricting
	 * a slightly-larger catalog parent down to it via
	 * {@link RecursiveMaterialiser#projectInto}. Because projection is the only
	 * <em>downward</em> operator, a gain at one shape can unlock a smaller one in
	 * the next round, so the pass iterates until a round writes no new win (or the
	 * round cap is hit).
	 */
	private static void runProjectionClosure(List<Shape> ordered,
			List<BlockSplitSearch.NamedBase> pool, FieldAwareLookup lookup, RunSpec spec) {
		System.out.println();
		System.out.println("=== PHASE 3: downward projection closure (multipass) ===");
		RecursiveClosureSota sota = new RecursiveClosureSota(lookup, pool, true, true);
		RecursiveMaterialiser mat = new RecursiveMaterialiser(
				lookup, pool, sota, spec.schemesRoot, true, true, /* improveExisting */ true);
		mat.setStrategies(spec.strategies);
		mat.setAllowCommutative(spec.allowCommutative());
		// Largest first: project a parent before the children that may project from it.
		List<Shape> desc = new ArrayList<>(ordered);
		desc.sort(BY_MAX_AXIS_ASC.reversed());

		int round = 0;
		while (round < PROJECTION_MAX_ROUNDS) {
			round++;
			int wins = 0;
			ProgressMonitor monitor = ProgressMonitor.builder()
					.total(desc.size())
					.weight(i -> etaWeight(desc.get((int) i)))
					.logEveryMillis(5_000)
					.label("SchemeSweep:projection round " + round)
					.build();
			monitor.start();
			for (int idx = 0; idx < desc.size(); idx++) {
				Shape s = desc.get(idx);
				try {
					int before = lookup.findRank(s.n(), s.m(), s.p());
					Optional<RecursiveMaterialiser.Result> r = mat.projectInto(s.n(), s.m(), s.p());
					if (r.isPresent()) {
						wins++;
						System.out.printf("  [projection win] ⟨%d,%d,%d⟩ %s → %d%n",
								s.n(), s.m(), s.p(), before <= 0 ? "—" : Integer.toString(before), r.get().alg().r);
					}
				} catch (RuntimeException e) {
					log.warn("projection failed at ⟨{},{},{}⟩: {}", s.n(), s.m(), s.p(), e.toString());
				} finally {
					monitor.tick(idx);
				}
			}
			monitor.done();
			System.out.printf("=== Projection round %d: %d new win(s) ===%n", round, wins);
			if (wins == 0) break;
		}
		System.out.printf("Projection closure finished after %d round(s).%n", round);
	}

	/** Format a nanosecond duration as {@code HhMmS.fs} (or smaller units). */
	private static String formatDuration(long nanos) {
		long ms = nanos / 1_000_000L;
		if (ms < 1000) return ms + "ms";
		long s = ms / 1000;
		ms = ms % 1000;
		if (s < 60) return String.format("%d.%03ds", s, ms);
		long m = s / 60;
		s = s % 60;
		if (m < 60) return String.format("%dm%02ds", m, s);
		long h = m / 60;
		m = m % 60;
		return String.format("%dh%02dm%02ds", h, m, s);
	}

	// ---------------------------------------------------------------------
	// CLI parsing
	// ---------------------------------------------------------------------

	private static RunSpec parseArgs(String[] args) {
		Mode mode = Mode.EVALUATE;
		// MANDATORY — no default. A sweep's results are field-specific (a scheme valid
		// over F₂ need not lift to Q/R; the rank landscape differs per field, e.g.
		// ⟨4,4,4⟩ is 47/F₂ vs 49/Z). "Best scheme whatever the field" is meaningless and
		// produces cross-field phantoms (an F₂ base leaking into a Q search). The caller
		// MUST declare the algebra (Q | Z | R | C | F2 | F3); we reject early otherwise.
		String field = null;
		List<String> configNames = new ArrayList<>();
		List<int[]> extraShapes = new ArrayList<>();
		int cubicMin = OFF, cubicMax = OFF;
		int noncubicMin = OFF, noncubicMax = OFF;
		boolean anyScopeFlag = false;
		Path out = null;
		Path schemesRoot = DEFAULT_SCHEMES_ROOT;
		int threads = 0;
		Boolean balancedOnlyOverride = null;
		Integer maxImbalanceOverride = null;
		Integer maxBaseDimOverride = null;
		eu.solven.matmul.SymmetryTransforms.InternalOrbitMode orbitModeOverride = null;
		Integer maxCombinationsOverride = null;
		Integer maxPaddingOverride = null;
		boolean skipMaterialise = false;
		boolean projectionOnly = false;
		List<int[]> baseShapes = new ArrayList<>();
		String baseFilter = null;
		String baseShape = null;
		int minMaxDim = 0;  // 0 = no band floor; else keep shapes with max(n,m,p) ≥ this
		int bandMin = OFF, bandMax = OFF;  // --band=lo-hi (or N): all shapes with maxDim in [lo,hi]
		// --strategies: restrict the compose phase to a subset of upward operators.
		// Default = all three; e.g. --strategies=serendipitous runs ONLY the
		// serendipitous-product operator (skips recombination/Kron/concat search
		// AND the downward projection closure).
		java.util.Set<String> strategies = java.util.Set.of(
				RecursiveMaterialiser.STRAT_RECOMBINATION,
				RecursiveMaterialiser.STRAT_SERENDIPITOUS,
				RecursiveMaterialiser.STRAT_PROJECTION);
		// Commutativity axis (orthogonal to the field): default NON-COMMUTATIVE — schemes
		// must lift to recursive matmul over an NC ring. --commutative=true allows
		// commutative-only schemes (Waksman/Rosowski/Makarov) as ingredients (scalar /
		// commutative-matmul sweeps).
		boolean allowCommutative = false;
		for (String raw : args) {
			if (raw.equals("--help") || raw.equals("-h")) {
				printUsage();
				System.exit(0);
			}
			if (!raw.startsWith("--")) {
				throw new IllegalArgumentException(
						"Unrecognised arg '" + raw + "' — use --help.");
			}
			int eq = raw.indexOf('=');
			if (eq < 0) {
				throw new IllegalArgumentException("Flag '" + raw + "' needs a value (--key=value).");
			}
			String key = raw.substring(2, eq);
			String value = raw.substring(eq + 1);
			switch (key) {
				case "mode" -> mode = Mode.valueOf(value.toUpperCase(java.util.Locale.ROOT));
				case "field" -> field = value;
				case "optimizer" -> {
					// flat (default) = AllocationOptimizer; assignment = partition+assignment B&B.
					switch (value.toLowerCase(java.util.Locale.ROOT)) {
						case "flat", "allocation" -> eu.solven.matmul.recombination.BlockSplitSearch.USE_ASSIGNMENT_OPTIMIZER = false;
						case "assignment", "partition" -> eu.solven.matmul.recombination.BlockSplitSearch.USE_ASSIGNMENT_OPTIMIZER = true;
						default -> throw new IllegalArgumentException(
								"--optimizer must be 'flat' or 'assignment', got '" + value + "'");
					}
				}
				case "assignmentMaxNodes" ->
						eu.solven.matmul.recombination.BlockSplitSearch.ASSIGNMENT_MAX_NODES = Long.parseLong(value);
				// Anytime budgets for the default flat AllocationOptimizer path. Make
				// large-base bands (≥11) terminate fast at a near-optimal rank instead
				// of paying the unbounded exact proof.
				case "maxNodes" ->
						eu.solven.matmul.recombination.BlockSplitSearch.ALLOC_MAX_NODES = Long.parseLong(value);
				case "stagnation" ->
						eu.solven.matmul.recombination.BlockSplitSearch.ALLOC_STAGNATION = Long.parseLong(value);
				// Write-gate for materialize (a TOP-LEVEL decision only — it does NOT
				// affect how sub-blocks are resolved inside a composition; see
				// RecursiveMaterialiser.resolveSubScheme). Two equivalent spellings:
				//   --only-if-missing[=true]  → fill: write a shape only when MISSING;
				//                               never touch an existing entry.
				//   --only-if-missing=false   → also rewrite an existing shape when a
				//                               Verifier-passing composition is strictly
				//                               better (closes re-derivable gaps, e.g.
				//                               ⟨3,3,18⟩=120, ⟨26,26,32⟩=11165).
				// `--improve` is the deprecated inverse alias (--improve=true == fill off).
				case "only-if-missing", "only_if_missing", "onlyifmissing" ->
						MATERIALIZE_IMPROVE = !(value.isBlank() || Boolean.parseBoolean(value));
				case "improve" -> MATERIALIZE_IMPROVE = value.isBlank() || Boolean.parseBoolean(value);
				case "buildable" -> {
					if (value.equalsIgnoreCase("required")) MATERIALIZE_BUILDABLE_REQUIRED = true;
					else if (value.equalsIgnoreCase("optional")) MATERIALIZE_BUILDABLE_REQUIRED = false;
					else throw new IllegalArgumentException("--buildable must be 'required' or 'optional', got '" + value + "'");
				}
				case "best-derived", "best_derived", "bestderived", "derive-best" ->
						MATERIALIZE_DERIVE_BEST = value.isBlank() || Boolean.parseBoolean(value);
				case "derive-all", "derive_all", "deriveall", "register-derived-anyway" -> {
					boolean on = value.isBlank() || Boolean.parseBoolean(value);
					MATERIALIZE_DERIVE_ALL = on;
					if (on) MATERIALIZE_DERIVE_BEST = true;  // derive-all implies best-derived
				}
				case "config" -> {
					for (String n : value.split(",")) {
						if (!n.isBlank()) configNames.add(n.trim());
					}
				}
				case "shape" -> {
					for (String s : value.split(",")) {
						if (s.isBlank()) continue;
						String[] parts = s.trim().split("x");
						if (parts.length != 3) {
							throw new IllegalArgumentException("Bad --shape '" + s + "' — expected NxMxP.");
						}
						extraShapes.add(new int[] {
								Integer.parseInt(parts[0]),
								Integer.parseInt(parts[1]),
								Integer.parseInt(parts[2]) });
					}
					anyScopeFlag = true;
				}
				case "shape-file", "shape_file", "shapefile" -> {
					try {
						for (String line : java.nio.file.Files.readAllLines(java.nio.file.Path.of(value))) {
							String s = line.trim();
							if (s.isBlank()) continue;
							// tolerate "NxMxP …extra…" rows — take the leading token
							String tok = s.split("\\s+")[0];
							String[] parts = tok.split("x");
							if (parts.length != 3) continue;
							extraShapes.add(new int[] {
									Integer.parseInt(parts[0]),
									Integer.parseInt(parts[1]),
									Integer.parseInt(parts[2]) });
						}
					} catch (java.io.IOException e) {
						throw new IllegalArgumentException("--shape-file unreadable: " + value, e);
					}
					anyScopeFlag = true;
				}
				case "cubic" -> {
					anyScopeFlag = true;
					if (value.equalsIgnoreCase("off")) {
						cubicMin = OFF; cubicMax = OFF;
					} else {
						int[] range = parseRange(value);
						cubicMin = range[0]; cubicMax = range[1];
					}
				}
				case "noncubic" -> {
					anyScopeFlag = true;
					if (value.equalsIgnoreCase("off")) {
						noncubicMin = OFF; noncubicMax = OFF;
					} else {
						int[] range = parseRange(value);
						noncubicMin = range[0]; noncubicMax = range[1];
					}
				}
				case "minMaxDim", "min-maxdim", "band-floor" -> minMaxDim = Integer.parseInt(value);
				case "band" -> {
					anyScopeFlag = true;
					int[] r = parseRange(value);  // "9" → 9-9, "2-9" → 2-9 (on maxDim)
					bandMin = r[0]; bandMax = r[1];
				}
				case "out" -> out = Path.of(value);
				case "schemes-root" -> schemesRoot = Path.of(value);
				case "threads" -> threads = Integer.parseInt(value);
				case "balancedOnly", "balanced-only" -> balancedOnlyOverride = Boolean.parseBoolean(value);
				case "maxImbalance", "max-imbalance" -> {
					maxImbalanceOverride = "off".equalsIgnoreCase(value)
							? PoolConfig.UNBOUNDED_IMBALANCE
							: Integer.parseInt(value);
				}
				case "maxBaseDim", "max-base-dim" -> maxBaseDimOverride = Integer.parseInt(value);
				case "orbitMode", "orbit-mode" -> orbitModeOverride =
						eu.solven.matmul.SymmetryTransforms.InternalOrbitMode.valueOf(
								value.toUpperCase(java.util.Locale.ROOT));
				case "maxCombinations", "max-combinations" -> maxCombinationsOverride =
						"off".equalsIgnoreCase(value) ? PoolConfig.UNBOUNDED_COMBINATIONS
								: Integer.parseInt(value);
				case "maxPadding", "max-padding" -> maxPaddingOverride = Integer.parseInt(value);
				case "dry-run", "dryRun", "skipMaterialise", "skip-materialise", "skip-materialize" ->
						skipMaterialise = Boolean.parseBoolean(value);
				case "projectionOnly", "projection-only" ->
						projectionOnly = Boolean.parseBoolean(value);
				case "strategies", "strategy" -> {
					java.util.Set<String> sel = new java.util.LinkedHashSet<>();
					for (String tok : value.split(",")) {
						String t = tok.trim().toLowerCase(java.util.Locale.ROOT);
						if (t.isBlank()) continue;
						// Accept a few friendly aliases.
						switch (t) {
							case "serendipitous", "serendipity", "serendip" ->
									sel.add(RecursiveMaterialiser.STRAT_SERENDIPITOUS);
							case "recombination", "recomb", "kron", "concat", "split" ->
									sel.add(RecursiveMaterialiser.STRAT_RECOMBINATION);
							case "projection", "project", "proj" ->
									sel.add(RecursiveMaterialiser.STRAT_PROJECTION);
							default -> throw new IllegalArgumentException(
									"--strategies: unknown token '" + t + "' (expected "
									+ "serendipitous|recombination|projection).");
						}
					}
					if (sel.isEmpty()) {
						throw new IllegalArgumentException("--strategies needs ≥1 token.");
					}
					strategies = java.util.Set.copyOf(sel);
				}
				case "commutative", "comm" -> allowCommutative = value.isBlank() || Boolean.parseBoolean(value);
				case "baseFilter", "base-filter" -> baseFilter = value;
				case "baseShape", "base-shape" -> baseShape = value;
				case "base", "bases" -> {
					// --base=NxMxP[,NxMxP…]: use the given catalog scheme(s) as the
					// outer recombination pool (each oriented over its axis-perm orbit)
					// instead of the default rootPool. Lets you sweep "with base 3x3x6".
					for (String s : value.split(",")) {
						if (s.isBlank()) continue;
						String[] parts = s.trim().split("x");
						if (parts.length != 3) {
							throw new IllegalArgumentException("Bad --base '" + s + "' — expected NxMxP.");
						}
						baseShapes.add(new int[] {
								Integer.parseInt(parts[0]),
								Integer.parseInt(parts[1]),
								Integer.parseInt(parts[2]) });
					}
				}
				default -> throw new IllegalArgumentException("Unknown flag '--" + key + "'. Use --help.");
			}
		}
		if (!anyScopeFlag && extraShapes.isEmpty()) {
			// Default scope depends on mode.
			switch (mode) {
				case EVALUATE -> { cubicMin = 2; cubicMax = MAX_CUBIC_DIM; noncubicMin = 2; noncubicMax = MAX_NONCUBIC_DIM; }
				case MATERIALIZE -> { cubicMin = 4; cubicMax = MAX_CUBIC_DIM; noncubicMin = 3; noncubicMax = MAX_CUBIC_DIM; }
				case CLOSURE -> { cubicMin = 4; cubicMax = MAX_CUBIC_DIM; /* non-cubic off for closure default */ }
			}
		}
		if (out == null) {
			out = (mode == Mode.EVALUATE) ? DEFAULT_REPORT_OUT
					: Path.of("target/scheme-sweep/" + mode.name().toLowerCase(java.util.Locale.ROOT) + ".log");
		}
		// Field is mandatory — reject the sweep BEFORE building any index.
		if (field == null || field.isBlank()) {
			throw new IllegalArgumentException("--field is MANDATORY (no default): a sweep is "
					+ "field-specific. Pass one of Q | Z | R | C | F2 | F3 (e.g. --field=Q). "
					+ "Looking for a 'best scheme whatever the field' is meaningless and leaks "
					+ "cross-field phantoms (an F₂ base into a Q search).");
		}
		// Validate it's a known algebra (throws with a clear message otherwise).
		eu.solven.matmul.algebra.Field.fromTag(field);
		return new RunSpec(mode, field, configNames, cubicMin, cubicMax,
				noncubicMin, noncubicMax, extraShapes, out, schemesRoot, threads,
				balancedOnlyOverride, maxImbalanceOverride,
				maxBaseDimOverride, orbitModeOverride,
				maxCombinationsOverride, maxPaddingOverride,
				skipMaterialise, projectionOnly, baseFilter, baseShape, baseShapes, minMaxDim, bandMin, bandMax,
				strategies, allowCommutative);
	}

	private static int[] parseRange(String s) {
		int dash = s.indexOf('-');
		if (dash < 0) {
			int n = Integer.parseInt(s);
			return new int[] { n, n };
		}
		return new int[] { Integer.parseInt(s.substring(0, dash)),
				Integer.parseInt(s.substring(dash + 1)) };
	}

	private static void printUsage() {
		System.out.println("SchemeSweep — unified shape-sweep driver");
		System.out.println();
		System.out.println("Flags (all optional):");
		System.out.println("  --mode=evaluate|materialize|closure   default: evaluate");
		System.out.println("    evaluate    : multi-config A/B; no writes; markdown report");
		System.out.println("    materialize : single config, parallel, RecursiveMaterialiser; writes catalog files");
		System.out.println("    closure     : single config, sequential, flat materialiser; iterates to fixed point");
		System.out.println();
		System.out.println("  --field=Q                       field for FieldAwareLookup (default Q)");
		System.out.println("  --config=simple,auditAxisFlip   comma-separated config names (or 'all')");
		System.out.println("                                  evaluate: many OK; materialize/closure: exactly one");
		System.out.println("  --shape=NxMxP[,NxMxP…]          explicit shape list");
		System.out.println("  --cubic=2-16  or  --cubic=8     cubic ⟨n,n,n⟩ range (or 'off')");
		System.out.println("  --noncubic=2-8                  non-cubic max-axis range (or 'off')");
		System.out.println("  --out=PATH                      report/log path (default under target/scheme-sweep/)");
		System.out.println("  --schemes-root=PATH             catalog write root (materialize/closure)");
		System.out.println("                                  default " + DEFAULT_SCHEMES_ROOT);
		System.out.println("  --threads=N                     materialize parallelism (default: cpus-1)");
		System.out.println("  --dry-run=true                  closure: run only Phase 1 search, no disk writes (alias: --skip-materialise)");
		System.out.println("  --base=NxMxP[,NxMxP…]           use these catalog scheme(s) as the outer");
		System.out.println("                                  recombination pool (each axis-oriented),");
		System.out.println("                                  instead of the default rootPool. e.g. 3x3x6");
		System.out.println("  --projection-only=true          closure: run ONLY the downward projection closure");
		System.out.println("                                  (Phase 3) — skips the upward search/materialise");
		System.out.println();
		System.out.println("Search-shape overrides (override the PoolConfig preset):");
		System.out.println("  --balancedOnly=true|false       force balanced allocations on each axis");
		System.out.println("                                  (true = old DIS09 heuristic; default false)");
		System.out.println("  --maxImbalance=N                cap on max(parts) − min(parts) per axis");
		System.out.println("                                  (default unbounded; e.g. =3 admits 9+8+0..3-imbalance)");
		System.out.println();
		System.out.println("Known configs: " + ALL_CONFIGS.stream().map(NamedConfig::name).toList());
	}

	private static List<NamedConfig> resolveConfigs(List<String> names) {
		if (names.isEmpty()) {
			return DEFAULT_CONFIGS;
		}
		if (names.size() == 1 && "all".equalsIgnoreCase(names.get(0))) {
			return ALL_CONFIGS;
		}
		Set<String> wanted = new LinkedHashSet<>(names);
		List<NamedConfig> out = new ArrayList<>();
		for (NamedConfig nc : ALL_CONFIGS) {
			if (wanted.remove(nc.name)) out.add(nc);
		}
		if (!wanted.isEmpty()) {
			throw new IllegalArgumentException("Unknown config name(s): " + wanted
					+ ". Known: " + ALL_CONFIGS.stream().map(NamedConfig::name).toList()
					+ " or 'all'.");
		}
		if (out.isEmpty()) {
			throw new IllegalArgumentException("No configs selected from --config=" + names);
		}
		return out;
	}

	private static List<Shape> buildShapeSet(RunSpec spec) {
		Set<Shape> seen = new LinkedHashSet<>();
		if (spec.cubicMin > 0) {
			for (int n = spec.cubicMin; n <= spec.cubicMax; n++) {
				seen.add(Shape.of(n, n, n));
			}
		}
		if (spec.noncubicMin > 0) {
			for (int n = spec.noncubicMin; n <= spec.noncubicMax; n++) {
				for (int m = n; m <= spec.noncubicMax; m++) {
					for (int p = m; p <= spec.noncubicMax; p++) {
						if (n == m && m == p) continue;
						seen.add(Shape.of(n, m, p));
					}
				}
			}
		}
		// --band=lo-hi : every shape (cubic + noncubic) whose maxDim is in [lo,hi],
		// minDim ≥ 2. The expressive way to sweep a max-dim band in one flag
		// (replaces the noncubic-range + minMaxDim-floor combo).
		if (spec.bandMin > 0) {
			for (int n = 2; n <= spec.bandMax; n++) {
				for (int m = n; m <= spec.bandMax; m++) {
					for (int p = m; p <= spec.bandMax; p++) {
						if (p >= spec.bandMin) seen.add(Shape.of(n, m, p));
					}
				}
			}
		}
		for (int[] xs : spec.extraShapes) {
			seen.add(Shape.of(xs[0], xs[1], xs[2]));
		}
		// Band floor (--minMaxDim): keep only shapes whose largest axis is ≥ the
		// floor. Combined with the cubic/noncubic upper bound this carves out a
		// max-dim band, e.g. --noncubic=2-9 --cubic=9-9 --minMaxDim=9 → exactly
		// the maxDim-9 shapes, so successive bands can be swept independently.
		if (spec.minMaxDim > 0) {
			seen.removeIf(s -> Math.max(s.n(), Math.max(s.m(), s.p())) < spec.minMaxDim);
		}
		return new ArrayList<>(seen);
	}

	// ---------------------------------------------------------------------
	// EVALUATE-mode report rendering (markdown)
	// ---------------------------------------------------------------------

	private static String renderReport(List<Shape> shapes,
			Map<Shape, Map<String, Long>> results,
			List<NamedConfig> configs) {
		boolean hasSimple = configs.stream().anyMatch(nc -> "simple".equals(nc.name));
		StringWriter buf = new StringWriter();
		try (PrintWriter pw = new PrintWriter(buf)) {
			pw.println("# SchemeSweep evaluate coverage summary");
			pw.println();
			pw.println("Generated by `eu.solven.matmul.docs.SchemeSweep --mode=evaluate`.");
			pw.println();
			pw.println("Compares the rank that `BlockSplitSearch.findBestStrategy`");
			pw.println("returns under each named PoolConfig across the shape sweep.");
			pw.println("Rows where a non-`simple` config beats `simple` are evidence");
			pw.println("that the richer config genuinely closes catalog gaps.");
			pw.println();
			pw.println("**Columns**: `simple` is the baseline; subsequent");
			pw.println("columns report (rank, Δ from simple). A negative Δ means");
			pw.println("the config found a strictly better recipe than `simple`.");
			pw.println();

			pw.println("## Headline");
			pw.println();
			if (!hasSimple) {
				pw.println("_(headline skipped — `simple` baseline not in this run)_");
				pw.println();
			} else {
				pw.println("| config | shapes-improved | average Δ on improved | best Δ |");
				pw.println("|---|---:|---:|---:|");
				for (NamedConfig nc : configs) {
					if (nc.name.equals("simple")) continue;
					int improved = 0;
					long sumDelta = 0;
					long bestDelta = 0;
					for (Shape s : shapes) {
						Long simple = results.get(s).get("simple");
						Long here = results.get(s).get(nc.name);
						if (simple == null || here == null || simple < 0 || here < 0) continue;
						long delta = here - simple;
						if (delta < 0) {
							improved++;
							sumDelta += delta;
							if (delta < bestDelta) bestDelta = delta;
						}
					}
					String avg = improved == 0 ? "—" : String.format("%.1f", sumDelta / (double) improved);
					pw.printf("| `%s` | %d | %s | %d |%n",
							nc.name, improved, avg, bestDelta);
				}
				pw.println();
			}

			if (hasSimple) {
				pw.println("## Per-shape detail (only rows where some config beats `simple`)");
			} else {
				pw.println("## Per-shape detail (all rows; no `simple` baseline)");
			}
			pw.println();
			pw.print("| shape |");
			for (NamedConfig nc : configs) pw.print(" " + nc.name + " |");
			pw.println();
			pw.print("|---|");
			for (int i = 0; i < configs.size(); i++) pw.print("---:|");
			pw.println();
			int rowsShown = 0;
			for (Shape s : shapes) {
				Long simple = results.get(s).get("simple");
				if (hasSimple) {
					boolean anyBetter = false;
					for (NamedConfig nc : configs) {
						if (nc.name.equals("simple")) continue;
						Long here = results.get(s).get(nc.name);
						if (here != null && simple != null && here >= 0 && simple >= 0 && here < simple) {
							anyBetter = true;
							break;
						}
					}
					if (!anyBetter) continue;
				}
				rowsShown++;
				pw.printf("| %s |", s);
				for (NamedConfig nc : configs) {
					Long v = results.get(s).get(nc.name);
					if (v == null || v < 0) {
						pw.print(" — |");
					} else if (!hasSimple || nc.name.equals("simple")) {
						pw.print(" " + v + " |");
					} else {
						long delta = v - simple;
						String sgn = delta < 0 ? "" : (delta > 0 ? "+" : "±");
						pw.printf(" %d (%s%d) |", v, sgn, delta);
					}
				}
				pw.println();
			}
			pw.println();
			if (hasSimple) {
				pw.printf("_%d shapes with at least one improvement._%n", rowsShown);
			} else {
				pw.printf("_%d shapes._%n", rowsShown);
			}
		}
		return buf.toString();
	}
}
