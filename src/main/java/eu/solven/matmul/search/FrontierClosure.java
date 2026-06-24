package eu.solven.matmul.search;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Set;

import tools.jackson.databind.JsonNode;
import tools.jackson.databind.json.JsonMapper;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.recombination.Recombination;
import eu.solven.matmul.catalog.SchemeIO;

/**
 * Frontier closure search: propagate from low-ω catalog seeds via Kronecker
 * products to discover catalog improvements without brute-forcing every
 * (shape, allocation, base) triple.
 *
 * <p><b>Algorithm</b>:
 * <ol>
 *   <li>Read {@code docs/catalog.json} to get the catalog's best-known
 *       rank per shape (sorted canonical form ⟨a,b,c⟩ with a ≤ b ≤ c).</li>
 *   <li>Compute ω_eff = 3·ln(R) / ln(a·b·c) per shape.</li>
 *   <li>Push every catalog leaf into a {@link PriorityQueue} keyed by
 *       ω ascending (low ω = "high quality seed" = explored first).</li>
 *   <li>Pop the lowest-ω seed; try Kronecker with each base in a fixed
 *       expansion pool (Strassen ⟨2,2,2⟩=7, axis-split ⟨1,2,1⟩=2, Laderman
 *       ⟨3,3,3⟩=23, AT-Z ⟨2,2,3⟩=11, Smirnov ⟨3,3,6⟩=40). Compute the
 *       Kronecker target shape (a·n, b·m, c·p) sorted to canonical form
 *       and the predicted rank R_seed · R_base.</li>
 *   <li>If the target shape has max-dim ≤ {@link #maxDim} AND predicted
 *       rank is strictly less than the catalog SOTA there, record the
 *       improvement and push the new (target, predicted) onto the queue
 *       so it can seed further propagation.</li>
 *   <li>Stop when the queue is empty.</li>
 * </ol>
 *
 * <p><b>Use case</b>: re-discover known FMM-Lille / Smirnov / Pan-TA
 * catalog entries automatically by propagating from low-ω leaves, building
 * a derivation chain in the process. If the closure produces many of the
 * existing FMM entries from a small seed set, that validates the approach
 * for closing further gaps.
 *
 * <p><b>Limitations of this first cut</b>:
 * <ul>
 *   <li>Only Kronecker products — no Recombination (BlockSplit) or
 *       Pan-TA fusion yet.</li>
 *   <li>No materialisation — only rank comparison. Successful predictions
 *       need a separate materialiser pass to write factor matrices.</li>
 *   <li>Single field (Q) — F2 / commutative variants would need
 *       per-field SOTA tables.</li>
 * </ul>
 *
 * <p>Run via:
 * <pre>{@code
 *   mvn -q exec:java -Dexec.mainClass=eu.solven.matmul.search.FrontierClosure
 * }</pre>
 */
public final class FrontierClosure {

	private static final File CATALOG_JSON = new File("docs/catalog.json");

	/** Cap on max(n,m,p) for target shapes. Avoids exploding into very large shapes. */
	public final int maxDim;
	/** Cap on the seed's ω — seeds above this are not popped. */
	public final double maxOmega;
	/** When > 0: use dynamic base pool from catalog up to this max base-dim. */
	public final int maxBaseDim;
	/** Enable recombination operator (Strassen-style allocation + peel). */
	public boolean enableRecombination = false;
	/** Max imbalance for recombination allocations. */
	public int recombMaxImbalance = 3;
	/** Max output peel (padding before peel) per axis. */
	public int recombMaxPadding = 1;
	/** Cap recombination rounds. */
	public int maxRecombRounds = 5;

	public FrontierClosure(int maxDim, double maxOmega) {
		this(maxDim, maxOmega, 0);
	}

	public FrontierClosure(int maxDim, double maxOmega, int maxBaseDim) {
		this.maxDim = maxDim;
		this.maxOmega = maxOmega;
		this.maxBaseDim = maxBaseDim;
	}

	// ── Expansion pool: small bases used as Kronecker factors ──
	//
	// Each is a (shape, rank, label) triple. We don't load JSON files for
	// these — they're hard-coded constants that every catalog has.
	public record Base(int n, int m, int p, long rank, String label) {}

	public static List<Base> defaultBasePool() {
		List<Base> pool = new ArrayList<>();
		pool.add(new Base(1, 2, 1, 2, "axis-split⟨1,2,1⟩"));
		pool.add(new Base(2, 1, 1, 2, "axis-split⟨2,1,1⟩"));
		pool.add(new Base(1, 1, 2, 2, "axis-split⟨1,1,2⟩"));
		pool.add(new Base(2, 2, 2, 7, "Strassen⟨2,2,2⟩=7"));
		pool.add(new Base(2, 2, 3, 11, "AT-Z⟨2,2,3⟩=11"));
		pool.add(new Base(2, 3, 3, 15, "HK⟨2,3,3⟩=15"));
		pool.add(new Base(3, 3, 3, 23, "Laderman⟨3,3,3⟩=23"));
		pool.add(new Base(3, 3, 6, 40, "Smirnov⟨3,3,6⟩=40"));
		return pool;
	}

	/**
	 * Dynamic base pool: every leaf-quality (non-Composed) catalog scheme up
	 * to {@code maxBaseDim}, including rectangular shapes. The 3 axis-split
	 * bases are added regardless (they're algebraically trivial but
	 * structurally essential — many derivations route through them).
	 *
	 * <p>A "leaf" is a scheme whose source string does NOT start with
	 * "Composed-" — these are hand-imported / hand-constructed primaries,
	 * not the output of our materialiser. Includes Strassen, Winograd,
	 * Laderman, AT-Z, AlphaEvolve, Smirnov ⟨3,P,Q⟩, HK ⟨2,b,c⟩, Pan TA,
	 * Schwartz-Zwecher, etc.
	 *
	 * <p>S₃-orbit: each leaf shape contributes all permutations of (n,m,p)
	 * at the same rank, so we cover e.g. ⟨3,3,6⟩, ⟨3,6,3⟩, ⟨6,3,3⟩.
	 */
	public List<Base> dynamicBasePool(int maxBaseDim) throws IOException {
		List<Base> pool = new ArrayList<>();
		// Axis-splits always included
		pool.add(new Base(1, 2, 1, 2, "axis-split⟨1,2,1⟩"));
		pool.add(new Base(2, 1, 1, 2, "axis-split⟨2,1,1⟩"));
		pool.add(new Base(1, 1, 2, 2, "axis-split⟨1,1,2⟩"));

		JsonMapper M = JsonMapper.builder().build();
		JsonNode root = M.readTree(CATALOG_JSON);
		JsonNode schemes = root.get("schemes");
		// Best rank per canonical (sorted) shape, restricted to leaf NC R/Q/Z.
		Map<String, long[]> best = new HashMap<>();
		Map<String, String> bestLabel = new HashMap<>();
		for (JsonNode s : schemes) {
			JsonNode fmt = s.get("format");
			if (fmt == null || fmt.size() != 3) continue;
			int n = fmt.get(0).asInt(), m = fmt.get(1).asInt(), p = fmt.get(2).asInt();
			if (n < 2 || m < 2 || p < 2) continue;
			if (Math.max(n, Math.max(m, p)) > maxBaseDim) continue;
			String field = s.has("field") ? s.get("field").asText() : "";
			if (!field.equals("R/Q/Z") && !field.equals("Q") && !field.equals("Z")
					&& !field.equals("R") && !field.equals("ZT")) continue;
			if (s.has("commutative") && s.get("commutative").asBoolean(false)) continue;
			String file = s.has("file") ? s.get("file").asText() : "";
			if (file.contains("_commutative")) continue;
			String src = s.has("source") ? s.get("source").asText() : "";
			String srcLow = src.toLowerCase();
			if (srcLow.startsWith("derived-")) continue; // derived, not leaf
			if (srcLow.startsWith("waksman") || srcLow.startsWith("rosowski")
					|| srcLow.startsWith("makarov") || srcLow.startsWith("islam")
					|| srcLow.startsWith("smith ") || srcLow.startsWith("probert")) continue;
			long rank = s.get("rank").asLong();
			int[] sorted = new int[] { n, m, p };
			Arrays.sort(sorted);
			String c = canon(sorted);
			long[] prev = best.get(c);
			if (prev == null || rank < prev[3]) {
				best.put(c, new long[] { sorted[0], sorted[1], sorted[2], rank });
				bestLabel.put(c, String.format("%s⟨%d,%d,%d⟩=%d",
						src.split(" ")[0], sorted[0], sorted[1], sorted[2], rank));
			}
		}
		// Emit each leaf in canonical (sorted) form. Permutations are not
		// added explicitly — the propagation step itself canonicalises target
		// shapes via canon(), so ⟨3,6,12⟩ via (⟨3,3,6⟩, axis-split) and via
		// (⟨3,6,3⟩, …) etc. all hit the same map key.
		for (var e : best.entrySet()) {
			long[] v = e.getValue();
			pool.add(new Base((int) v[0], (int) v[1], (int) v[2], v[3], bestLabel.get(e.getKey())));
		}
		return pool;
	}

	// ── Internal state ──

	/** Best known rank per canonical (sorted) shape (LIVE; updated by improvements). */
	private final Map<int[], Long> sotaByShape = new HashMap<>();
	/** Original catalog SOTA snapshot — captured at loadCatalog() time, never updated. */
	private final Map<String, Long> originalCatalogSota = new HashMap<>();
	/** Best frontier-predicted rank per canonical shape, regardless of whether it beat catalog. */
	private final Map<String, Long> bestFrontierByCanon = new HashMap<>();
	/** Set of shapes we've already pushed into the queue (to avoid infinite revisit). */
	private final Set<int[]> seen = new HashSet<>();

	private static class IntArrayKey {
		// We use int[] but need value-based equals/hashCode for the maps above.
		// Wrap to make this explicit if needed; for now we accept that the
		// raw int[] keys WILL NOT dedupe via equals — we use a string canonical
		// form instead.
	}

	private static String canon(int[] s) {
		int[] c = s.clone();
		Arrays.sort(c);
		return c[0] + "x" + c[1] + "x" + c[2];
	}

	/** Compute ω_eff = 3·ln(R) / ln(n·m·p). */
	public static double omegaOf(int n, int m, int p, long rank) {
		double denom = Math.log((double) n * m * p);
		if (denom <= 0) return 3.0;
		return 3.0 * Math.log(rank) / denom;
	}

	public record Seed(int n, int m, int p, long rank, double omega, String lineage)
			implements Comparable<Seed> {
		@Override
		public int compareTo(Seed o) {
			int c = Double.compare(omega, o.omega);
			if (c != 0) return c;
			return Long.compare(rank, o.rank);
		}
	}

	public record Improvement(int n, int m, int p, long predicted, long previousSota,
			String viaSeed, String viaBase) {
		@Override
		public String toString() {
			return String.format("⟨%d,%d,%d⟩=%d (was %d, Δ=%+d) via %s ⊗ %s",
					n, m, p, predicted, previousSota, predicted - previousSota,
					viaSeed, viaBase);
		}
	}

	/** Load the catalog SOTA from {@code docs/catalog.json}. */
	private void loadCatalog() throws IOException {
		JsonMapper M = JsonMapper.builder().build();
		JsonNode root = M.readTree(CATALOG_JSON);
		JsonNode schemes = root.get("schemes");
		Map<String, Long> bestByCanon = new HashMap<>();
		Map<String, int[]> shapeByCanon = new HashMap<>();
		for (JsonNode s : schemes) {
			JsonNode fmt = s.get("format");
			if (fmt == null || fmt.size() != 3) continue;
			int n = fmt.get(0).asInt(), m = fmt.get(1).asInt(), p = fmt.get(2).asInt();
			if (n < 2 || m < 2 || p < 2) continue;
			if (Math.max(n, Math.max(m, p)) > maxDim) continue;
			// Field filter: only R/Q/Z (broad characteristic-0 bucket).
			// F2-only and C-only ranks DO NOT lift to recursive matmul over
			// rationals, so they can't be used as Kronecker seeds.
			String field = s.has("field") ? s.get("field").asText() : "";
			if (!field.equals("R/Q/Z") && !field.equals("Q") && !field.equals("Z")
					&& !field.equals("R") && !field.equals("ZT")) {
				continue;
			}
			// Skip commutative-only schemes — they do NOT lift to NC matmul
			// under Kronecker.
			if (s.has("commutative") && s.get("commutative").asBoolean(false)) continue;
			String file = s.has("file") ? s.get("file").asText() : "";
			if (file.contains("_commutative")) continue;
			String src = s.has("source") ? s.get("source").asText().toLowerCase() : "";
			if (src.startsWith("waksman") || src.startsWith("rosowski")
					|| src.startsWith("makarov") || src.startsWith("islam")
					|| src.startsWith("smith ") || src.startsWith("probert")) {
				continue; // commutative-only families
			}
			long rank = s.get("rank").asLong();
			String c = canon(new int[] { n, m, p });
			Long prev = bestByCanon.get(c);
			if (prev == null || rank < prev) {
				bestByCanon.put(c, rank);
				int[] sorted = new int[] { n, m, p };
				Arrays.sort(sorted);
				shapeByCanon.put(c, sorted);
			}
		}
		for (var e : bestByCanon.entrySet()) {
			sotaByShape.put(shapeByCanon.get(e.getKey()), e.getValue());
			originalCatalogSota.put(e.getKey(), e.getValue());
		}
	}

	private Long getSota(int[] shape) {
		int[] sorted = shape.clone();
		Arrays.sort(sorted);
		String c = canon(sorted);
		for (var entry : sotaByShape.entrySet()) {
			if (canon(entry.getKey()).equals(c)) return entry.getValue();
		}
		return null;
	}

	private void setSota(int[] shape, long rank) {
		int[] sorted = shape.clone();
		Arrays.sort(sorted);
		// Remove any existing key for this canon
		String c = canon(sorted);
		sotaByShape.entrySet().removeIf(e -> canon(e.getKey()).equals(c));
		sotaByShape.put(sorted, rank);
	}

	// ── Recombination (outer-base + allocation + peel) operator ──
	//
	// For each candidate target shape ⟨n,m,p⟩ within maxDim, try each small
	// outer base in {Strassen, Winograd, axflip-cousins, Laderman, AT-Z}
	// with allocations splitting n,m,p across the base's grid. Each outer
	// mult M_k becomes a sub-matmul whose rank is looked up in the FRONTIER's
	// current best-known map. Sum the sub-ranks. Optionally over-allocate
	// (padding) and peel zero-outputs (DIS09 γ5 reduction).
	//
	// This subsumes pure Kronecker (set all alloc slots equal = standard
	// Kronecker) AND adds asymmetric allocations + peel, which is where
	// most Strassen-style catalog wins come from.

	/** Small outer bases for recombination. Each loaded from disk via SchemeIO. */
	private List<NonCubicBilinearAlgorithm> outerBases = null;
	/** Cache of per-(scheme, allocation, peel) shape multisets. */
	private final CachedMultisets multisetCache = new CachedMultisets();

	private void loadOuterBases() {
		if (outerBases != null) return;
		outerBases = new ArrayList<>();
		String[] paths = {
				"src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json",
				"src/main/resources/schemes/known/section2/winograd_1971-2x2x2_m7_a24.json",
				"src/main/resources/schemes/curated/section2/solven_winograd_cousin_axflip1-2x2x2_m7_a24.json",
				"src/main/resources/schemes/known/section3/alphatensor_Z-2x2x3_m11_a25.json",
				"src/main/resources/schemes/known/section3/laderman_1976-3x3x3_m23_a98.json",
		};
		// Resolve by content (byHint parses shape+source token) — the 2026-06 rename
		// broke these literal paths; without this every outer base was silently dropped.
		for (String p : paths) {
			try {
				outerBases.add(SchemeIO.read(eu.solven.matmul.catalog.SchemeResolver.byHint(p)));
			} catch (Exception e) {
				System.err.println("[recombination] skip missing base " + p + ": " + e.getMessage());
			}
		}
	}

	/** SotaResolver backed by the frontier's current bestFrontierByCanon ∪ catalog SOTA. */
	private Recombination.SotaResolver frontierResolver() {
		return (a, b, c) -> {
			if (a == 0 || b == 0 || c == 0) return 0;
			int[] sorted = new int[] { a, b, c };
			Arrays.sort(sorted);
			String key = canon(sorted);
			Long frontier = bestFrontierByCanon.get(key);
			Long cat = originalCatalogSota.get(key);
			long best = Long.MAX_VALUE / 4;
			if (frontier != null) best = Math.min(best, frontier);
			if (cat != null) best = Math.min(best, cat);
			if (best >= Long.MAX_VALUE / 4) return a * b * c; // naive fallback
			return (int) best;
		};
	}

	/** Enumerate compositions of {@code budget} into {@code blocks} non-negative parts. */
	private static List<int[]> compositions(int blocks, int budget) {
		List<int[]> out = new ArrayList<>();
		int[] cur = new int[blocks];
		compRec(blocks, budget, cur, 0, out);
		return out;
	}

	private static void compRec(int blocks, int budget, int[] cur, int idx, List<int[]> out) {
		if (idx == blocks - 1) {
			cur[idx] = budget;
			out.add(cur.clone());
			return;
		}
		for (int i = 0; i <= budget; i++) {
			cur[idx] = i;
			compRec(blocks, budget - i, cur, idx + 1, out);
		}
	}

	/** Try recombining {@code outer} onto target ⟨n,m,p⟩. Returns the best predicted rank found. */
	private long recombineBest(NonCubicBilinearAlgorithm outer, int n, int m, int p,
			int maxImbalance, int maxPadding, Recombination.SotaResolver sota) {
		long best = Long.MAX_VALUE / 4;
		int baseA = outer.n, baseB = outer.m, baseC = outer.p;
		// Skip if target dim is less than base dim (no valid alloc).
		if (n < baseA || m < baseB || p < baseC) return best;
		for (int dN = 0; dN <= maxPadding; dN++) {
			for (int dM = 0; dM <= maxPadding; dM++) {
				for (int dP = 0; dP <= maxPadding; dP++) {
					int padN = n + dN, padM = m + dM, padP = p + dP;
					if (padN < baseA || padM < baseB || padP < baseC) continue;
					List<int[]> allocsA = compositions(baseA, padN);
					List<int[]> allocsB = compositions(baseB, padM);
					List<int[]> allocsC = compositions(baseC, padP);
					for (int[] aA : allocsA) {
						if (imbalance(aA) > maxImbalance) continue;
						if (aA[aA.length - 1] < dN) continue;
						int[] peelA = (dN == 0) ? null : tailPeel(aA.length, dN);
						for (int[] aB : allocsB) {
							if (imbalance(aB) > maxImbalance) continue;
							if (aB[aB.length - 1] < dM) continue;
							int[] peelB = (dM == 0) ? null : tailPeel(aB.length, dM);
							for (int[] aC : allocsC) {
								if (imbalance(aC) > maxImbalance) continue;
								if (aC[aC.length - 1] < dP) continue;
								int[] peelC = (dP == 0) ? null : tailPeel(aC.length, dP);
								// Use cached per-scheme symbolic templates instead of
								// re-scanning factor matrices on every call (task #125).
								SchemeTemplates templates = SchemeTemplates.forScheme(outer);
								long total = templates.totalRank(aA, aB, aC, peelA, peelB, peelC, sota);
								if (total < best) best = total;
							}
						}
					}
				}
			}
		}
		return best;
	}

	private static int imbalance(int[] a) {
		int min = a[0], max = a[0];
		for (int v : a) { if (v < min) min = v; if (v > max) max = v; }
		return max - min;
	}

	private static int[] tailPeel(int blocks, int d) {
		int[] p = new int[blocks];
		p[blocks - 1] = d;
		return p;
	}

	/**
	 * One pass over all target shapes ⟨n,m,p⟩ ≤ maxDim, trying recombination
	 * with each outer base + bounded allocations + peel. Returns the list of
	 * improvements found in this pass.
	 */
	private List<Improvement> recombinationPass(int maxImbalance, int maxPadding) {
		loadOuterBases();
		Recombination.SotaResolver sota = frontierResolver();
		List<Improvement> found = new ArrayList<>();
		for (int n = 2; n <= maxDim; n++) {
			for (int m = n; m <= maxDim; m++) {
				for (int p = m; p <= maxDim; p++) {
					int[] sorted = new int[] { n, m, p };
					String c = canon(sorted);
					long currentBest = Long.MAX_VALUE / 4;
					Long fr = bestFrontierByCanon.get(c);
					if (fr != null) currentBest = fr;
					Long cat = originalCatalogSota.get(c);
					if (cat != null) currentBest = Math.min(currentBest, cat);
					long bestRecomb = Long.MAX_VALUE / 4;
					NonCubicBilinearAlgorithm bestBase = null;
					for (NonCubicBilinearAlgorithm outer : outerBases) {
						long r = recombineBest(outer, n, m, p, maxImbalance, maxPadding, sota);
						if (r < bestRecomb) { bestRecomb = r; bestBase = outer; }
					}
					if (bestRecomb < currentBest && bestBase != null) {
						long was = currentBest >= Long.MAX_VALUE / 4 ? Long.MAX_VALUE : currentBest;
						Improvement i = new Improvement(n, m, p, bestRecomb, was,
								String.format("frontier(⟨%d,%d,%d⟩)", n, m, p),
								String.format("recombine[%dx%dx%d]", bestBase.n, bestBase.m, bestBase.p));
						found.add(i);
						bestFrontierByCanon.put(c, bestRecomb);
						setSota(sorted, bestRecomb);
					}
				}
			}
		}
		return found;
	}

	public List<Improvement> run() throws IOException {
		loadCatalog();
		System.out.printf("Catalog loaded: %d distinct shapes (maxDim ≤ %d)%n",
				sotaByShape.size(), maxDim);

		PriorityQueue<Seed> queue = new PriorityQueue<>();
		Set<String> queuedCanon = new HashSet<>();
		for (var e : sotaByShape.entrySet()) {
			int[] s = e.getKey();
			long r = e.getValue();
			double w = omegaOf(s[0], s[1], s[2], r);
			if (w > maxOmega) continue;
			queue.add(new Seed(s[0], s[1], s[2], r, w, "catalog"));
			queuedCanon.add(canon(s));
		}
		System.out.printf("Initial seeds: %d (ω ≤ %.3f)%n", queue.size(), maxOmega);

		List<Base> bases = (maxBaseDim > 0) ? dynamicBasePool(maxBaseDim) : defaultBasePool();
		System.out.printf("Base pool: %d entries (maxBaseDim=%d)%n", bases.size(), maxBaseDim);
		// Build a flat snapshot of all known shape→rank pairs for the
		// axis-concat operator. We need to iterate this in the propagation
		// step alongside the multiplicative `bases`.
		List<int[]> snapshot = new ArrayList<>();
		for (var e : sotaByShape.entrySet()) {
			int[] s = e.getKey();
			long r = e.getValue();
			snapshot.add(new int[] { s[0], s[1], s[2], (int) r });
		}
		List<Improvement> improvements = new ArrayList<>();
		int popped = 0;
		while (!queue.isEmpty()) {
			Seed seed = queue.poll();
			popped++;
			// Re-check: the seed's rank may have been beaten by a later
			// improvement while it was in the queue.
			Long currentBest = getSota(new int[] { seed.n, seed.m, seed.p });
			if (currentBest == null || currentBest < seed.rank) continue;

			// --- Multiplicative (Kronecker) propagation ---
			for (Base b : bases) {
				int tn = seed.n * b.n();
				int tm = seed.m * b.m();
				int tp = seed.p * b.p();
				if (Math.max(tn, Math.max(tm, tp)) > maxDim) continue;
				long predicted = seed.rank * b.rank();
				maybeImprove(improvements, queue, queuedCanon, tn, tm, tp, predicted,
						"⟨" + seed.n + "," + seed.m + "," + seed.p + "⟩=" + seed.rank,
						"⊗ " + b.label());
			}
			// --- Axis-concat propagation (plain axis-split) ---
			// For each other catalog entry sharing 2 of the 3 dimensions with
			// the seed, the third axis can be concatenated:
			//   R(⟨n,m,p1+p2⟩) ≤ R(⟨n,m,p1⟩) + R(⟨n,m,p2⟩)
			// This is the trivial upper bound from running two independent
			// matmuls that share the unsplit axes — not τ-theorem (which
			// requires independent tensor sums with shared/cancelling
			// structure).
			//
			// Example: ⟨3,3,16⟩ ≤ ⟨3,3,4⟩=29 + ⟨3,3,12⟩=80 = 109. The matrices
			// share A (3×3) and split B (3×16) into B1 (3×4) + B2 (3×12);
			// output C = [A·B1 | A·B2] (column-concatenated). Cost is just
			// the two sub-matmuls.
			for (int[] other : snapshot) {
				int on = other[0], om = other[1], op = other[2];
				long or_ = other[3];
				// Concat along p (n,m must match)
				if (seed.n == on && seed.m == om) {
					int tn = seed.n, tm = seed.m, tp = seed.p + op;
					if (tp >= 2 && Math.max(tn, Math.max(tm, tp)) <= maxDim) {
						maybeImprove(improvements, queue, queuedCanon, tn, tm, tp,
								seed.rank + or_,
								"⟨" + seed.n + "," + seed.m + "," + seed.p + "⟩=" + seed.rank,
								"⊕ ⟨" + on + "," + om + "," + op + "⟩=" + or_ + " (concat-p)");
					}
				}
				// Concat along m
				if (seed.n == on && seed.p == op) {
					int tn = seed.n, tm = seed.m + om, tp = seed.p;
					if (tm >= 2 && Math.max(tn, Math.max(tm, tp)) <= maxDim) {
						maybeImprove(improvements, queue, queuedCanon, tn, tm, tp,
								seed.rank + or_,
								"⟨" + seed.n + "," + seed.m + "," + seed.p + "⟩=" + seed.rank,
								"⊕ ⟨" + on + "," + om + "," + op + "⟩=" + or_ + " (concat-m)");
					}
				}
				// Concat along n
				if (seed.m == om && seed.p == op) {
					int tn = seed.n + on, tm = seed.m, tp = seed.p;
					if (tn >= 2 && Math.max(tn, Math.max(tm, tp)) <= maxDim) {
						maybeImprove(improvements, queue, queuedCanon, tn, tm, tp,
								seed.rank + or_,
								"⟨" + seed.n + "," + seed.m + "," + seed.p + "⟩=" + seed.rank,
								"⊕ ⟨" + on + "," + om + "," + op + "⟩=" + or_ + " (concat-n)");
					}
				}
			}
		}
		System.out.printf("Kronecker+concat pass: popped %d seeds, found %d improvements%n",
				popped, improvements.size());

		// ── Recombination passes ──
		// After Kronecker+concat converges, sweep targets with outer-base
		// recombination + allocations + peel. Each improvement is fed back
		// as a new seed; we iterate until no recombination improvement.
		if (enableRecombination) {
			int round = 0;
			while (round < maxRecombRounds) {
				round++;
				long t0 = System.nanoTime();
				List<Improvement> recombFound = recombinationPass(recombMaxImbalance, recombMaxPadding);
				long ms = (System.nanoTime() - t0) / 1_000_000;
				System.out.printf("Recombination pass %d: %d improvements in %d ms%n",
						round, recombFound.size(), ms);
				if (recombFound.isEmpty()) break;
				improvements.addAll(recombFound);
				// Push each as a new seed so Kronecker+concat can propagate further
				for (Improvement i : recombFound) {
					int[] sortedT = new int[] { i.n, i.m, i.p };
					Arrays.sort(sortedT);
					String tc = canon(sortedT);
					if (queuedCanon.add(tc)) {
						double w = omegaOf(sortedT[0], sortedT[1], sortedT[2], i.predicted);
						if (w <= maxOmega) {
							queue.add(new Seed(sortedT[0], sortedT[1], sortedT[2],
									i.predicted, w, "recomb"));
						}
					}
				}
				// Re-drain the queue for any cascading Kronecker+concat wins
				int reKronecker = 0;
				while (!queue.isEmpty()) {
					Seed seed = queue.poll();
					reKronecker++;
					Long currentBest = getSota(new int[] { seed.n, seed.m, seed.p });
					if (currentBest == null || currentBest < seed.rank) continue;
					for (Base b : bases) {
						int tn = seed.n * b.n();
						int tm = seed.m * b.m();
						int tp = seed.p * b.p();
						if (Math.max(tn, Math.max(tm, tp)) > maxDim) continue;
						long predicted = seed.rank * b.rank();
						maybeImprove(improvements, queue, queuedCanon, tn, tm, tp, predicted,
								"⟨" + seed.n + "," + seed.m + "," + seed.p + "⟩=" + seed.rank,
								"⊗ " + b.label());
					}
					for (int[] other : snapshot) {
						int on = other[0], om = other[1], op = other[2];
						long or_ = other[3];
						if (seed.n == on && seed.m == om) {
							int tn = seed.n, tm = seed.m, tp = seed.p + op;
							if (tp >= 2 && Math.max(tn, Math.max(tm, tp)) <= maxDim) {
								maybeImprove(improvements, queue, queuedCanon, tn, tm, tp,
										seed.rank + or_,
										"⟨" + seed.n + "," + seed.m + "," + seed.p + "⟩=" + seed.rank,
										"⊕ ⟨" + on + "," + om + "," + op + "⟩=" + or_ + " (concat-p)");
							}
						}
						if (seed.n == on && seed.p == op) {
							int tn = seed.n, tm = seed.m + om, tp = seed.p;
							if (tm >= 2 && Math.max(tn, Math.max(tm, tp)) <= maxDim) {
								maybeImprove(improvements, queue, queuedCanon, tn, tm, tp,
										seed.rank + or_,
										"⟨" + seed.n + "," + seed.m + "," + seed.p + "⟩=" + seed.rank,
										"⊕ ⟨" + on + "," + om + "," + op + "⟩=" + or_ + " (concat-m)");
							}
						}
						if (seed.m == om && seed.p == op) {
							int tn = seed.n + on, tm = seed.m, tp = seed.p;
							if (tn >= 2 && Math.max(tn, Math.max(tm, tp)) <= maxDim) {
								maybeImprove(improvements, queue, queuedCanon, tn, tm, tp,
										seed.rank + or_,
										"⟨" + seed.n + "," + seed.m + "," + seed.p + "⟩=" + seed.rank,
										"⊕ ⟨" + on + "," + om + "," + op + "⟩=" + or_ + " (concat-n)");
							}
						}
					}
				}
				System.out.printf("  cascaded Kronecker+concat: %d additional seed pops%n", reKronecker);
			}
		}
		System.out.printf("TOTAL: %d improvements%n", improvements.size());
		return improvements;
	}

	/** Record an improvement at target (tn,tm,tp) if {@code predicted} beats current SOTA. */
	private void maybeImprove(List<Improvement> improvements, PriorityQueue<Seed> queue,
			Set<String> queuedCanon, int tn, int tm, int tp, long predicted,
			String viaSeed, String viaBaseOrConcat) {
		// Track best frontier prediction regardless of whether it improves catalog.
		int[] sortedTargetA = new int[] { tn, tm, tp };
		Arrays.sort(sortedTargetA);
		String tcA = canon(sortedTargetA);
		Long prevFrontier = bestFrontierByCanon.get(tcA);
		if (prevFrontier == null || predicted < prevFrontier) {
			bestFrontierByCanon.put(tcA, predicted);
		}
		Long prevSota = getSota(new int[] { tn, tm, tp });
		if (prevSota != null && predicted >= prevSota) return;
		long was = prevSota == null ? Long.MAX_VALUE : prevSota;
		improvements.add(new Improvement(tn, tm, tp, predicted, was, viaSeed, viaBaseOrConcat));
		setSota(new int[] { tn, tm, tp }, predicted);
		double w = omegaOf(tn, tm, tp, predicted);
		if (w <= maxOmega) {
			int[] sortedTarget = new int[] { tn, tm, tp };
			Arrays.sort(sortedTarget);
			String tc = canon(sortedTarget);
			if (queuedCanon.add(tc)) {
				queue.add(new Seed(sortedTarget[0], sortedTarget[1], sortedTarget[2],
						predicted, w, viaSeed + " " + viaBaseOrConcat));
			}
		}
	}

	public static void main(String[] args) throws Exception {
		int maxDim = eu.solven.matmul.catalog.CatalogPolicy.MATERIALISE_MAX_DIM;
		double maxOmega = 3.0; // no ω cap by default
		int maxBaseDim = 0; // 0 = use hardcoded defaultBasePool
		boolean enableRecomb = false;
		int recombMaxImbalance = 3;
		int recombMaxPadding = 1;
		for (int i = 0; i < args.length; i++) {
			if (i + 1 < args.length) {
				if (args[i].equals("--max-dim")) { maxDim = Integer.parseInt(args[i + 1]); continue; }
				if (args[i].equals("--max-omega")) { maxOmega = Double.parseDouble(args[i + 1]); continue; }
				if (args[i].equals("--max-base-dim")) { maxBaseDim = Integer.parseInt(args[i + 1]); continue; }
				if (args[i].equals("--recomb-imbalance")) { recombMaxImbalance = Integer.parseInt(args[i + 1]); continue; }
				if (args[i].equals("--recomb-padding")) { recombMaxPadding = Integer.parseInt(args[i + 1]); continue; }
			}
			if (args[i].equals("--recombination")) enableRecomb = true;
		}
		FrontierClosure fc = new FrontierClosure(maxDim, maxOmega, maxBaseDim);
		fc.enableRecombination = enableRecomb;
		fc.recombMaxImbalance = recombMaxImbalance;
		fc.recombMaxPadding = recombMaxPadding;
		long t0 = System.nanoTime();
		List<Improvement> impr = fc.run();
		long ms = (System.nanoTime() - t0) / 1_000_000;
		System.out.printf("Done in %d ms.%n%n", ms);

		impr.sort(Comparator.<Improvement, Integer>comparing(
						i -> Math.max(i.n, Math.max(i.m, i.p)))
				.thenComparing(i -> i.n)
				.thenComparing(i -> i.m)
				.thenComparing(i -> i.p));

		System.out.println("Top improvements (sorted by max-dim, then shape):");
		for (Improvement i : impr) {
			if (i.previousSota == Long.MAX_VALUE) {
				System.out.printf("  ⟨%d,%d,%d⟩=%d  NEW  via %s ⊗ %s%n",
						i.n, i.m, i.p, i.predicted, i.viaSeed, i.viaBase);
			} else {
				System.out.printf("  ⟨%d,%d,%d⟩=%d  (was %d, Δ=%+d)  via %s ⊗ %s%n",
						i.n, i.m, i.p, i.predicted, i.previousSota,
						i.predicted - i.previousSota, i.viaSeed, i.viaBase);
			}
		}

		// ── Coverage summary ──
		// For every shape in the original catalog, classify what frontier
		// produced at that shape.
		int better = 0, tie = 0, worse = 0, missed = 0;
		long sumWorseGap = 0;
		long maxWorseGap = 0;
		List<String> worstMissedShapes = new ArrayList<>();
		for (var e : fc.originalCatalogSota.entrySet()) {
			String canonKey = e.getKey();
			long sota = e.getValue();
			Long frontierBest = fc.bestFrontierByCanon.get(canonKey);
			if (frontierBest == null) {
				missed++;
				continue;
			}
			if (frontierBest < sota) better++;
			else if (frontierBest.longValue() == sota) tie++;
			else {
				worse++;
				long gap = frontierBest - sota;
				sumWorseGap += gap;
				maxWorseGap = Math.max(maxWorseGap, gap);
				worstMissedShapes.add(String.format("⟨%s⟩ sota=%d frontier=%d (gap=%+d)",
						canonKey, sota, frontierBest, gap));
			}
		}
		int total = fc.originalCatalogSota.size();
		System.out.println();
		System.out.println("─── Coverage vs catalog SOTA ───");
		System.out.printf("  total catalog shapes: %d%n", total);
		System.out.printf("  BETTER (frontier beat catalog):  %d (%.1f%%)%n",
				better, 100.0 * better / total);
		System.out.printf("  TIE    (frontier matched):       %d (%.1f%%)%n",
				tie, 100.0 * tie / total);
		System.out.printf("  WORSE  (frontier above catalog): %d (%.1f%%)  avg gap=%+.1f  max gap=%+d%n",
				worse, 100.0 * worse / total,
				worse > 0 ? (double) sumWorseGap / worse : 0.0,
				maxWorseGap);
		System.out.printf("  MISSED (frontier never visited): %d (%.1f%%)%n",
				missed, 100.0 * missed / total);

		System.out.println();
		System.out.println("Top 20 worst frontier gaps (where frontier most degraded vs catalog):");
		worstMissedShapes.sort(Comparator.comparing(s ->
				-Long.parseLong(s.replaceAll(".*gap=\\+", "").replaceAll("\\).*", ""))));
		for (int i = 0; i < Math.min(20, worstMissedShapes.size()); i++) {
			System.out.println("  " + worstMissedShapes.get(i));
		}
	}
}
