package eu.solven.matmul.docs;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.search.flip.FlipGraphWalk;
import eu.solven.matmul.search.flip.FlipObjective;
import eu.solven.matmul.search.flip.FlipObjectives;
import eu.solven.matmul.search.flip.FlipScheme;
import eu.solven.matmul.search.flip.MetaFlipWalk;
import lombok.extern.slf4j.Slf4j;

/**
 * Flip-graph walk driver (in-Java port of the Kauers–Moosbauer / Perminov
 * flip-graph engines, {@code search.flip}). Unlike the published engines, the
 * default objectives target the catalog's STRUCTURE metrics — bud-richness
 * (serendipitous-product fuel) and projection margin (downward-parent
 * strength) — under a rank-first lexicographic cost; pure rank descent is also
 * available.
 *
 * <pre>
 *   mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.FlipGraphSweep \
 *       -Dexec.args="--shape=3x3x3 --objective=bud --steps=500000 --walks=4 --write"
 * </pre>
 *
 * <ul>
 *   <li>{@code --shape=NxMxP} (required) — target format.</li>
 *   <li>{@code --objective=bud|margin|rank} (default {@code bud}).</li>
 *   <li>{@code --seed-file=path} — explicit seed scheme JSON; default is the
 *       catalog-best integer ({@code Z}) scheme at the shape. Integer seeds
 *       only (the walk is exact-integer).</li>
 *   <li>{@code --steps} (default 200000), {@code --walks} (default 4,
 *       rng seeds {@code rng..rng+walks−1}), {@code --rng} (default 42).</li>
 *   <li>{@code --cap} (default 1 = ternary/ZT alphabet; 0 = unbounded Z),
 *       {@code --worse} (worse-accept prob, default 0.01),
 *       {@code --plateau} (default 5000), {@code --rank-above} (split rank
 *       headroom, default 2).</li>
 *   <li>{@code --objective=weighted --w-rank=10 --w-bud=1 --w-margin=0} —
 *       tradeoff objective where rank is NOT lexicographic (rank+1 is worth
 *       paying when it buys ≥ w-rank/w-bud structure points). Pair with
 *       {@code --split-prob=0.01} (lets the walk buy rank) and optionally
 *       {@code --no-full-reduce} (merges stay cost-gated instead of
 *       automatic, so they don't eat freshly-built bud pairs).</li>
 *   <li>{@code --meta} — cross-format walk ({@link MetaFlipWalk}): extend /
 *       project meta-moves with {@code --extend-prob}, {@code --project-prob},
 *       {@code --min-dim}, {@code --max-dim}, {@code --max-rank}; reports the
 *       best-found per visited shape, write-gated against the catalog-best
 *       integer scheme under the same objective.</li>
 *   <li>{@code --write} — persist the best walk result under
 *       {@code schemes/derived/sectionN/} when it strictly improves the
 *       objective vs the seed (fixed-format) / vs the catalog (meta), verified
 *       at the write boundary.</li>
 * </ul>
 *
 * <p>Honesty: every reported value is a <strong>bound</strong> (best found by
 * a heuristic anytime walk), never an optimum. Phase 2 (meta-moves across
 * formats: extend / project / merge / product, Kauers–Wood arXiv:2510.19787)
 * plugs into the same walk as kick-time moves.</p>
 */
@Slf4j
public class FlipGraphSweep {

	public static void main(String[] args) throws IOException {
		Map<String, String> opts = parse(args);
		String shape = require(opts, "shape");
		String[] dims = shape.split("x");
		int n = Integer.parseInt(dims[0]);
		int m = Integer.parseInt(dims[1]);
		int p = Integer.parseInt(dims[2]);

		FlipObjective objective = switch (opts.getOrDefault("objective", "bud")) {
			case "bud" -> FlipObjectives.maxBudScore();
			case "margin" -> FlipObjectives.maxProjectionMargin();
			case "rank" -> FlipObjectives.minRank();
			// Tradeoff regime: rank is NOT lexicographic. weighted(10,1,0) means
			// "rank+1 is worth paying when it buys ≥10 bud points". Pair with
			// --split-prob so the walk can actually buy rank, and consider
			// --no-full-reduce so merges don't eat the bud pairs it builds.
			case "weighted" -> FlipObjectives.weighted(
					Long.parseLong(opts.getOrDefault("w-rank", "10")),
					Long.parseLong(opts.getOrDefault("w-bud", "1")),
					Long.parseLong(opts.getOrDefault("w-margin", "0")));
			// Direct predicted-product cost against the catalog rank oracle —
			// use this (not bud) for harvest walks; see ProbeFlipBudHarvest.
			case "seren" -> {
				String[] inner = opts.getOrDefault("inner", "2x2x2").split("x");
				yield FlipObjectives.serendipitous(new FieldAwareLookup(Field.Q),
						Integer.parseInt(inner[0]), Integer.parseInt(inner[1]),
						Integer.parseInt(inner[2]));
			}
			// Parameter-free serendipity potential: inner = the base's own shape
			// (predicts base⊗base — the recursive-squaring / ω quantity).
			case "seren-self" -> FlipObjectives.selfSerendipitous(new FieldAwareLookup(Field.Q));
			// Direct projected cost against a concrete target — exchange rate of
			// rank vs margin is exactly 1, so pair with --max-rank-above > 0 and
			// --split-prob: rank+δ pays whenever it buys margin > δ.
			case "project" -> {
				String[] t = opts.getOrDefault("target", "").split("x");
				if (t.length != 3) {
					throw new IllegalArgumentException("--objective=project needs --target=NxMxP");
				}
				yield FlipObjectives.projectedTo(Integer.parseInt(t[0]),
						Integer.parseInt(t[1]), Integer.parseInt(t[2]),
						Integer.parseInt(opts.getOrDefault("project-delta", "2")));
			}
			default -> throw new IllegalArgumentException(
					"--objective must be bud|margin|rank|weighted|seren|seren-self|project, got "
							+ opts.get("objective"));
		};
		long steps = Long.parseLong(opts.getOrDefault("steps", "200000"));
		int walks = Integer.parseInt(opts.getOrDefault("walks", "4"));
		long rng = Long.parseLong(opts.getOrDefault("rng", "42"));
		int cap = Integer.parseInt(opts.getOrDefault("cap", "1"));
		double worse = Double.parseDouble(opts.getOrDefault("worse", "0.01"));
		long plateau = Long.parseLong(opts.getOrDefault("plateau", "5000"));
		int rankAbove = Integer.parseInt(opts.getOrDefault("rank-above", "2"));
		double splitProb = Double.parseDouble(opts.getOrDefault("split-prob", "0"));
		boolean fullReduce = !opts.containsKey("no-full-reduce");

		NonCubicBilinearAlgorithm seedAlg = loadSeed(opts, n, m, p);
		FlipScheme seed = FlipScheme.of(seedAlg);
		log.info("seed ⟨{},{},{}⟩: rank={} budScore={} margin={} maxAbsCoef={} ({})",
				n, m, p, seed.rank(), FlipObjectives.budScore(seed),
				FlipObjectives.projectionMargin(seed), seed.maxAbsCoefficient(),
				objective.describe());

		if (opts.containsKey("meta")) {
			runMeta(opts, seed, seedAlg, objective, steps, walks, rng, cap, worse, splitProb,
					fullReduce);
			return;
		}

		FlipGraphWalk.Result best = null;
		for (int wlk = 0; wlk < walks; wlk++) {
			FlipGraphWalk.Config cfg = new FlipGraphWalk.Config(
					steps, rng + wlk, cap, worse, plateau, rankAbove, Math.max(steps / 10, 1),
					splitProb, fullReduce);
			log.info("walk {}/{} (rngSeed={}) starting", wlk + 1, walks, rng + wlk);
			FlipGraphWalk.Result r = FlipGraphWalk.walk(seed, objective, cfg);
			log.info("walk {}/{}: bestCost={} rank={} budScore={} margin={} (steps={} accepted={}"
					+ " kicks={} restarts={})",
					wlk + 1, walks, r.bestCost(), r.best().rank(),
					FlipObjectives.budScore(r.best()), FlipObjectives.projectionMargin(r.best()),
					r.steps(), r.accepted(), r.kicks(), r.restarts());
			if (best == null || r.bestCost() < best.bestCost()) {
				best = r;
			}
		}

		FlipScheme winner = best.best();
		FieldAwareLookup qLookup = new FieldAwareLookup(Field.Q);
		long[] savings = FlipObjectives.serendipitySavingByAxis(
				winner, qLookup, winner.n, winner.m, winner.p);
		log.info("BEST-FOUND (bound, not an optimum): rank {}→{} budScore {}→{} margin {}→{}"
				+ " self-savings(σU,σV,σW)=({},{},{})",
				seed.rank(), winner.rank(),
				FlipObjectives.budScore(seed), FlipObjectives.budScore(winner),
				FlipObjectives.projectionMargin(seed), FlipObjectives.projectionMargin(winner),
				savings[0], savings[1], savings[2]);

		if (opts.containsKey("write")) {
			if (!best.improvedOverSeed()) {
				log.info("no strict objective improvement over the seed — nothing written");
				return;
			}
			// Routing by objective family: pure rank winners go to the standard
			// derived tree; serendipity-structure winners are BUD BASES; projection
			// winners are MARGIN BASES (possibly worse-rank but projection-rich —
			// margin is INTRINSIC to the scheme, so the stamp is the target-free
			// per-axis triple; only the win question is target/catalog-relative).
			String obj = opts.getOrDefault("objective", "bud");
			String category = switch (obj) {
				case "rank" -> "derived";
				case "project", "margin" -> "margin-bases";
				default -> "bud-bases";
			};
			write(winner, seedAlg, category, savings);
		}
	}

	/**
	 * Meta mode: walks wander across formats (extend/project meta-moves), with
	 * per-shape bests merged over all walks. The write gate is honest by
	 * construction: a shape's result is persisted only when it strictly beats
	 * the catalog-best INTEGER scheme at that shape under the SAME objective.
	 */
	private static void runMeta(Map<String, String> opts, FlipScheme seed,
			NonCubicBilinearAlgorithm seedAlg, FlipObjective objective, long steps, int walks,
			long rng, int cap, double worse, double splitProb, boolean fullReduce)
			throws IOException {
		double extendProb = Double.parseDouble(opts.getOrDefault("extend-prob", "0.002"));
		double projectProb = Double.parseDouble(opts.getOrDefault("project-prob", "0.004"));
		int minDim = Integer.parseInt(opts.getOrDefault("min-dim", "2"));
		int maxDim = Integer.parseInt(opts.getOrDefault("max-dim", "6"));
		int maxRank = Integer.parseInt(opts.getOrDefault("max-rank", "400"));

		Map<String, MetaFlipWalk.ShapeBest> merged = new java.util.TreeMap<>();
		for (int wlk = 0; wlk < walks; wlk++) {
			MetaFlipWalk.Config cfg = new MetaFlipWalk.Config(steps, rng + wlk, cap, worse,
					splitProb, extendProb, projectProb, minDim, maxDim, maxRank, fullReduce,
					Math.max(steps / 10, 1));
			log.info("meta walk {}/{} (rngSeed={}) starting", wlk + 1, walks, rng + wlk);
			MetaFlipWalk.Result r = MetaFlipWalk.walk(seed, objective, cfg);
			log.info("meta walk {}/{}: {} shapes (accepted={} extends={} projects={} restarts={})",
					wlk + 1, walks, r.bestByShape().size(), r.accepted(), r.extendsTaken(),
					r.projectsTaken(), r.restarts());
			r.bestByShape().forEach((shape, sb) -> merged.merge(shape, sb,
					(a, b) -> a.cost() <= b.cost() ? a : b));
		}

		FieldAwareLookup zLookup = new FieldAwareLookup(Field.Z);
		FieldAwareLookup qLookup = new FieldAwareLookup(Field.Q);
		boolean write = opts.containsKey("write");
		String structureObjective = switch (opts.getOrDefault("objective", "bud")) {
			case "rank" -> "derived";
			case "project", "margin" -> "margin-bases";
			default -> "bud-bases";
		};
		for (Map.Entry<String, MetaFlipWalk.ShapeBest> e : merged.entrySet()) {
			FlipScheme s = e.getValue().scheme();
			Long catalogCost = zLookup.find(s.n, s.m, s.p)
					.map(alg -> objective.cost(FlipScheme.of(alg)))
					.orElse(null);
			boolean beatsCatalog = catalogCost == null || e.getValue().cost() < catalogCost;
			long[] savings = FlipObjectives.serendipitySavingByAxis(s, qLookup, s.n, s.m, s.p);
			log.info("⟨{}⟩ best-found (bound): rank={} budScore={} margin={} σ=({},{},{})"
					+ " cost={} catalogCost={}{}",
					e.getKey(), s.rank(), FlipObjectives.budScore(s),
					FlipObjectives.projectionMargin(s), savings[0], savings[1], savings[2],
					e.getValue().cost(), catalogCost == null ? "∅" : catalogCost,
					beatsCatalog ? "  ← BEATS catalog under this objective" : "");
			if (write && beatsCatalog) {
				write(s, seedAlg, structureObjective, savings);
			}
		}
	}

	private static NonCubicBilinearAlgorithm loadSeed(Map<String, String> opts,
			int n, int m, int p) throws IOException {
		if (opts.containsKey("seed-file")) {
			return SchemeIO.readBilinear(new File(opts.get("seed-file")));
		}
		// Catalog-best INTEGER scheme: the walk is exact-integer, so the Z chain
		// (not Q) is the right default lookup.
		return new FieldAwareLookup(Field.Z).find(n, m, p)
				.orElseThrow(() -> new IllegalStateException(
						"no integer catalog scheme at ⟨" + n + "," + m + "," + p
								+ "⟩ — pass --seed-file"));
	}

	private static void write(FlipScheme winner, NonCubicBilinearAlgorithm seedAlg,
			String category, long[] selfSavings) throws IOException {
		NonCubicBilinearAlgorithm alg = winner.toAlgorithm();
		// Write-boundary verification per catalog discipline: exact integer Brent
		// check (immune to coefficient growth) + the standard double-based oracle.
		if (!winner.isExactIntTensor() || !Verifier.isExactNonCubic(alg)) {
			throw new IllegalStateException(
					"walk result failed write-boundary verification — refusing to persist");
		}
		int maxDim = Math.max(alg.n, Math.max(alg.m, alg.p));
		// Bud bases (serendipity-structure winners) and margin bases (projection
		// winners) — both possibly worse-rank but structure-rich — get their own
		// TOP-LEVEL categories, siblings of known/derived/curated/constructed:
		// derived/ connotes exhaustive or deterministic derivation, while metaflip
		// output is heuristic (user 2026-06-11/12). Readers recurse the whole
		// schemes root, so both subtrees stay in every pool; filenames keep the
		// metaflip provenance note.
		boolean heuristicBase = !"derived".equals(category);
		Path dir = Path.of("src/main/resources/schemes/" + category + "/section" + maxDim);
		dir.toFile().mkdirs();
		int adds = Verifier.additionCount(alg);
		// Filenames are pure labels (content-driven catalog), but the note token
		// carries provenance: plain "metaflip" for the heuristic-base subtrees,
		// "derived_metaflip" only under derived/.
		String fn = SchemeIO.canonicalName(alg, heuristicBase ? "metaflip" : "derived_metaflip");
		File out = dir.resolve(fn).toFile();
		if (out.exists()) {
			log.info("already on disk: {}", out.getName());
			return;
		}
		SchemeIO.write(alg, out);
		// Integer walk output: Z-exact, hence Q/R/C by inclusion and F2/F3 by
		// reduction — same field set every integer derived scheme carries. The
		// seed's fields are NOT copied: the walk may only widen within integers.
		Map<String, Object> meta = new LinkedHashMap<>();
		meta.put("fields", List.of("F2", "F3", "Z", "Q", "R", "C"));
		meta.put("source", "metaflip");
		meta.put("seed_hash", SchemeIO.contentHash(seedAlg));
		if ("bud-bases".equals(category)) {
			// Per-axis priced savings (σU,σV,σW) at inner = own shape — the
			// multi-dimensional serendipity profile (a greedy bound, like every
			// bud figure; see FlipObjectives.serendipitySavingByAxis).
			meta.put("self_serendipity_savings", List.of(selfSavings[0], selfSavings[1],
					selfSavings[2]));
		}
		if ("margin-bases".equals(category)) {
			// Per-axis single-drop margins [μn, μm, μp] — INTRINSIC to the scheme
			// (supports only, no partner pricing), valid towards every target
			// below it; multi-drop margins are recomputable from the factors
			// (ProjectionSearch.projectedRank). Win-or-not is catalog-relative
			// and evaluated at query time, never stamped.
			int[] mu = eu.solven.matmul.catalog.ProjectionSearch.axisMargins(alg);
			meta.put("projection_margins", List.of(mu[0], mu[1], mu[2]));
		}
		SchemeIO.addFields(out, meta, true);
		log.info("wrote {} (rank={}, budScore={}, margin={}, σ=({},{},{}))", out.getName(),
				alg.r, FlipObjectives.budScore(winner), FlipObjectives.projectionMargin(winner),
				selfSavings[0], selfSavings[1], selfSavings[2]);
	}

	private static Map<String, String> parse(String[] args) {
		Map<String, String> opts = new LinkedHashMap<>();
		for (String a : args) {
			if (!a.startsWith("--")) {
				throw new IllegalArgumentException("unexpected arg: " + a);
			}
			int eq = a.indexOf('=');
			if (eq < 0) {
				opts.put(a.substring(2), "");
			} else {
				opts.put(a.substring(2, eq), a.substring(eq + 1));
			}
		}
		return opts;
	}

	private static String require(Map<String, String> opts, String key) {
		String v = opts.get(key);
		if (v == null || v.isEmpty()) {
			throw new IllegalArgumentException("--" + key + "=... is required");
		}
		return v;
	}
}
