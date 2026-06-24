package eu.solven.matmul.search.als;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.junit.jupiter.api.Tag;
import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.catalog.SchemeIO;

/**
 * Drives {@link Ternary2x2x2MultisetEnumerator}: brute-force the ternary
 * {@code ⟨2,2,2⟩} rank-7 space and bucket schemes by their recombination
 * multiset, cross-referenced against the catalog's named schemes (Strassen,
 * Winograd, …).
 */
class TestTernary2x2x2Multisets {

	private static final String SECTION2 = "src/main/resources/schemes/known/section2/";

	/** Compute the (9,8)³ multiset of every catalog ⟨2,2,2⟩ rank-7 scheme. */
	private static Map<String, String> catalogMultisets() throws Exception {
		Map<String, String> nameToMultiset = new LinkedHashMap<>();
		File[] files = new File(SECTION2).listFiles((d, n) -> n.endsWith(".json"));
		if (files == null) return nameToMultiset;
		java.util.Arrays.sort(files, java.util.Comparator.comparing(File::getName));
		for (File f : files) {
			NonCubicBilinearAlgorithm alg;
			try {
				alg = SchemeIO.readBilinear(f);
			} catch (Exception e) {
				continue;
			}
			if (alg.n != 2 || alg.m != 2 || alg.p != 2 || alg.r != 7) continue;
			int[][] terms = termsFrom(alg);
			nameToMultiset.put(f.getName(), Ternary2x2x2MultisetEnumerator.multisetKey(terms));
		}
		return nameToMultiset;
	}

	/** Pack a ⟨2,2,2⟩ scheme's U/V/W columns into the enumerator's term layout. */
	private static int[][] termsFrom(NonCubicBilinearAlgorithm alg) {
		double[][] U = alg.denseU(); // [4][7]
		double[][] V = alg.denseV();
		double[][] W = alg.denseW();
		int[][] terms = new int[7][12];
		for (int k = 0; k < 7; k++) {
			for (int a = 0; a < 4; a++) terms[k][a] = U[a][k] != 0 ? 1 : 0;
			for (int b = 0; b < 4; b++) terms[k][4 + b] = V[b][k] != 0 ? 1 : 0;
			for (int c = 0; c < 4; c++) terms[k][8 + c] = W[c][k] != 0 ? 1 : 0;
		}
		return terms;
	}

	/** Load a catalog ⟨2,2,2⟩ scheme as the orbit-enumerator's {@code [7][3][4]} layout (exact int coeffs). */
	private static int[][][] seedFrom(NonCubicBilinearAlgorithm alg) {
		double[][] U = alg.denseU(), V = alg.denseV(), W = alg.denseW();
		int[][][] seed = new int[7][3][4];
		for (int k = 0; k < 7; k++) {
			for (int a = 0; a < 4; a++) seed[k][0][a] = (int) Math.round(U[a][k]);
			for (int b = 0; b < 4; b++) seed[k][1][b] = (int) Math.round(V[b][k]);
			for (int c = 0; c < 4; c++) seed[k][2][c] = (int) Math.round(W[c][k]);
		}
		return seed;
	}

	/** Build a NonCubicBilinearAlgorithm from a {@code [7][3][4]} scheme (for verification). */
	private static NonCubicBilinearAlgorithm algFrom(int[][][] s) {
		double[][] U = new double[4][7], V = new double[4][7], W = new double[4][7];
		for (int k = 0; k < 7; k++) {
			for (int a = 0; a < 4; a++) U[a][k] = s[k][0][a];
			for (int b = 0; b < 4; b++) V[b][k] = s[k][1][b];
			for (int c = 0; c < 4; c++) W[c][k] = s[k][2][c];
		}
		return new NonCubicBilinearAlgorithm(2, 2, 2, U, V, W);
	}

	private static List<int[][][]> catalogSeeds() throws Exception {
		List<int[][][]> seeds = new ArrayList<>();
		File[] files = new File(SECTION2).listFiles((d, n) -> n.endsWith(".json"));
		if (files == null) return seeds;
		java.util.Arrays.sort(files, Comparator.comparing(File::getName));
		for (File f : files) {
			NonCubicBilinearAlgorithm alg;
			try {
				alg = SchemeIO.readBilinear(f);
			} catch (Exception e) {
				continue;
			}
			if (alg.n != 2 || alg.m != 2 || alg.p != 2 || alg.r != 7) continue;
			if (!Verifier.isExactNonCubic(alg)) continue; // skip commutative / non-NC (e.g. Waksman)
			seeds.add(seedFrom(alg));
		}
		return seeds;
	}

	@Test
	@Tag("slow")
	void orbit_sweep() throws Exception {
		// Map canonical multiset -> a catalog name realising it (for labelling).
		Map<String, String> canonToCatalog = new LinkedHashMap<>();
		Map<String, String> nameToCanon = new LinkedHashMap<>();
		File[] files = new File(SECTION2).listFiles((d, n) -> n.endsWith(".json"));
		java.util.Arrays.sort(files, Comparator.comparing(File::getName));
		for (File f : files) {
			NonCubicBilinearAlgorithm alg;
			try { alg = SchemeIO.readBilinear(f); } catch (Exception e) { continue; }
			if (alg.n != 2 || alg.m != 2 || alg.p != 2 || alg.r != 7) continue;
			if (!Verifier.isExactNonCubic(alg)) continue;
			String canon = Ternary2x2x2Orbit.canonicalMultisetKey(seedFrom(alg));
			canonToCatalog.putIfAbsent(canon, f.getName());
			nameToCanon.put(f.getName(), canon);
		}

		List<int[][][]> seeds = catalogSeeds();
		System.out.println();
		System.out.printf("=== Ternary-isotropy orbit sweep from %d catalog seeds (X,Y,Z ∈ {-1,0,1}) ===%n",
				seeds.size());

		long t0 = System.nanoTime();
		Ternary2x2x2Orbit.Result res = Ternary2x2x2Orbit.sweep(seeds, -1, 1);
		long ms = (System.nanoTime() - t0) / 1_000_000;

		// Correctness: every representative must be an exact ⟨2,2,2⟩ matmul.
		for (int[][][] rep : res.representatives.values()) {
			assertThat(Verifier.isExactNonCubic(algFrom(rep)))
					.as("orbit image must be exact matmul").isTrue();
		}

		System.out.printf("transforms=%,d  ternary-hits=%,d  distinct-schemes=%,d  distinct-canonical-multisets=%d  %,d ms%n",
				res.transformsTried, res.ternaryHits,
				res.schemeCounts.values().stream().mapToLong(Long::longValue).sum(),
				res.representatives.size(), ms);
		System.out.println();

		List<Map.Entry<String, Long>> ordered = new ArrayList<>(res.schemeCounts.entrySet());
		ordered.sort(Comparator.comparing(Map.Entry::getKey));
		System.out.println("Distinct CANONICAL recombination multisets (at (9,8)³ reference split):");
		int i = 0;
		for (var e : ordered) {
			i++;
			String label = canonToCatalog.getOrDefault(e.getKey(), "(not in catalog)");
			System.out.printf("  #%-2d  ×%-9d  %-58s  [%s]%n",
					i, e.getValue(), Ternary2x2x2Orbit.pretty(e.getKey()), label);
		}
		System.out.println();
		System.out.println("Catalog schemes → canonical multiset:");
		for (var e : nameToCanon.entrySet()) {
			System.out.printf("  %-52s  %s%n", e.getKey(), Ternary2x2x2Orbit.pretty(e.getValue()));
		}

		// Note: per-axis block-flip (reversing a 2-part split's order) is a COST
		// optimisation the recombination search performs — NOT a symmetry that
		// identifies distinct multisets of a *given* decomposition. We deliberately
		// do NOT mod it out: doing so is ill-posed here because the (9,8) reading
		// collapses "touches-big-block" and "touches-both-blocks" to size 9, and it
		// would wrongly merge Strassen with Winograd (which reach distinct bags —
		// e.g. Winograd reaches a 2·⟨8,8,8⟩ bag, Strassen never does).

		// Seed-independence: the orbit from Strassen alone should reach the same set.
		List<int[][][]> strassenOnly = new ArrayList<>();
		for (int[][][] s : seeds) {
			if (Ternary2x2x2Orbit.canonicalMultisetKey(s).equals(
					Ternary2x2x2Orbit.canonicalMultisetKey(seedFrom(SchemeIO.readBilinear(
							new File(SECTION2, "2x2x2-r7-strassen-db11bcc.json")))))) {
				strassenOnly.add(s);
				break;
			}
		}
		Ternary2x2x2Orbit.Result strassenRes = Ternary2x2x2Orbit.sweep(strassenOnly, -1, 1);
		System.out.printf("%nSeed-independence: Strassen-only orbit → %d S₃-canonical multisets "
				+ "(vs %d from all %d seeds)%n",
				strassenRes.representatives.size(), res.representatives.size(), seeds.size());

		// Pin the headline count cited by the paper (multisets.tex): 6 ternary
		// multisets within the ternary change-of-basis scope. A failure here means
		// either a regression or a genuine discovery — investigate, don't bump.
		assertThat(res.representatives).hasSize(6);
		assertThat(strassenRes.representatives.keySet())
				.as("seed-independence: Strassen-only orbit reaches the same multisets")
				.isEqualTo(res.representatives.keySet());
	}

	/**
	 * Robustness to the change-of-basis alphabet: widen {@code X,Y,Z} entries to
	 * {@code {-2..2}}. If the canonical-multiset set is unchanged, that is strong
	 * evidence the ternary {@code ⟨2,2,2⟩} space realises exactly these multisets.
	 * Heavier (~10⁸ transforms) — tagged slow, run on demand.
	 */
	@Test
	@Tag("slow")
	void orbit_robustness_wide() throws Exception {
		List<int[][][]> seeds = catalogSeeds();
		System.out.printf("%n=== Wide orbit sweep (X,Y,Z ∈ {-2..2}) from %d seeds ===%n", seeds.size());
		long t0 = System.nanoTime();
		Ternary2x2x2Orbit.Result res = Ternary2x2x2Orbit.sweep(seeds, -2, 2, false);
		long ms = (System.nanoTime() - t0) / 1_000_000;
		System.out.printf("transforms=%,d  ternary-hits=%,d  distinct-canonical-multisets=%d  %,d ms%n",
				res.transformsTried, res.ternaryHits, res.representatives.size(), ms);
		for (var e : res.representatives.entrySet()) {
			System.out.printf("  %s%n", Ternary2x2x2Orbit.pretty(e.getKey()));
		}
		assertThat(res.representatives).isNotEmpty();
	}

	@Test
	@Tag("slow")
	void symbolic_certified_Q_orbit() throws Exception {
		NonCubicBilinearAlgorithm strassen = SchemeIO.readBilinear(new File(SECTION2, "2x2x2-r7-strassen-db11bcc.json"));
		NonCubicBilinearAlgorithm winograd = SchemeIO.readBilinear(new File(SECTION2, "2x2x2-r7-winograd_1971-511df05.json"));
		int[][][] sStr = seedFrom(strassen);
		int[][][] sWin = seedFrom(winograd);

		System.out.println();
		System.out.println("=== Exact symbolic enumeration over the full GL₂(ℚ)³ orbit (all rank-7 ⟨2,2,2⟩) ===");

		// Completeness self-check: increasing the direction bound must NOT grow the set.
		int[] bounds = { 2, 3, 4 };
		java.util.Set<String> prev = null;
		for (int b : bounds) {
			java.util.Set<String> ms = Ternary2x2x2Orbit.certifiedQMultisets(sStr, b);
			int[] axes = Ternary2x2x2Orbit.perAxisPatternCounts(sStr, b);
			System.out.printf("  dirBound=%d → per-axis patterns n/m/p = %d/%d/%d, distinct canonical multisets = %d%n",
					b, axes[0], axes[1], axes[2], ms.size());
			if (prev != null) {
				assertThat(ms).as("multiset set must be stable as direction bound grows").isEqualTo(prev);
			}
			prev = ms;
		}

		// Seed-independence: Strassen and Winograd are one GL-orbit → identical Q-sets.
		java.util.Set<String> qStr = Ternary2x2x2Orbit.certifiedQMultisets(sStr, 4);
		java.util.Set<String> qWin = Ternary2x2x2Orbit.certifiedQMultisets(sWin, 4);
		System.out.printf("%n  seed-independence: Strassen-seed Q-set size=%d, Winograd-seed Q-set size=%d, equal=%b%n",
				qStr.size(), qWin.size(), qStr.equals(qWin));
		assertThat(qWin).as("Q-orbit multiset set is seed-independent (de Groote)").isEqualTo(qStr);

		// Compare to the ternary-reachable set (the 6 from the orbit sweep).
		Ternary2x2x2Orbit.Result ternary = Ternary2x2x2Orbit.sweep(catalogSeeds(), -1, 1);
		java.util.Set<String> ternarySet = ternary.representatives.keySet();
		System.out.printf("%n  ℚ-orbit distinct canonical multisets : %d%n", qStr.size());
		System.out.printf("  ternary-reachable canonical multisets: %d%n", ternarySet.size());
		System.out.printf("  ternary ⊆ ℚ ? %b%n", qStr.containsAll(ternarySet));
		// Pin the paper's headline numbers: 40 over ℚ, ternary set contained in it.
		assertThat(qStr).hasSize(40);
		assertThat(qStr).as("ternary-reachable multisets ⊆ ℚ-orbit multisets")
				.containsAll(ternarySet);

		System.out.println();
		System.out.println("All distinct canonical recombination multisets (ℚ-orbit, at (9,8)³ reference split):");
		List<String> sorted = new ArrayList<>(qStr);
		java.util.Collections.sort(sorted);
		int i = 0;
		for (String k : sorted) {
			i++;
			String tern = ternarySet.contains(k) ? "ternary" : "ℚ-only ";
			System.out.printf("  #%-2d [%s]  %s%n", i, tern, Ternary2x2x2Orbit.pretty(k));
		}
		assertThat(qStr).isNotEmpty();
	}

	@Test
	@Tag("slow")
	void debug_strassen_vs_winograd_flips() throws Exception {
		var sota = eu.solven.matmul.recombination.Recombination.catalogResolver(
				eu.solven.matmul.algebra.Algebra.nonCommutative(eu.solven.matmul.algebra.Field.R));
		int[] alloc = { 9, 8 };
		for (String name : new String[] { "2x2x2-r7-strassen-db11bcc.json", "2x2x2-r7-winograd_1971-511df05.json" }) {
			NonCubicBilinearAlgorithm alg = SchemeIO.readBilinear(new File(SECTION2, name));
			var supports = eu.solven.matmul.recombination.AnalyticalMaskSearch.SchemeSupports.extract(alg);
			System.out.printf("%n%s at (9,8)³ → ⟨17,17,17⟩, all 8 masks:%n", name);
			for (int mask = 0; mask < 8; mask++) {
				int[] aA = ((mask & 1) != 0) ? new int[] { 8, 9 } : alloc;
				int[] aB = ((mask & 2) != 0) ? new int[] { 8, 9 } : alloc;
				int[] aC = ((mask & 4) != 0) ? new int[] { 8, 9 } : alloc;
				int[][] shapes = eu.solven.matmul.recombination.AnalyticalMaskSearch.shapesAt(supports, aA, aB, aC);
				long cost = eu.solven.matmul.recombination.AnalyticalMaskSearch.costOf(shapes, sota);
				System.out.printf("  mask=%d cost=%d  %s%n", mask, cost,
						Ternary2x2x2Orbit.pretty(sortedKey(shapes)));
			}
		}
	}

	private static String sortedKey(int[][] shapes) {
		String[] s = new String[shapes.length];
		for (int i = 0; i < shapes.length; i++) s[i] = shapes[i][0] + "," + shapes[i][1] + "," + shapes[i][2];
		java.util.Arrays.sort(s);
		return String.join("|", s);
	}

	@Test
	@Tag("slow")
	void catalog_multisets() throws Exception {
		Map<String, String> cat = catalogMultisets();
		System.out.println();
		System.out.println("=== Catalog ⟨2,2,2⟩ rank-7 multisets at (9,8)³ ===");
		Map<String, String> multisetToName = new LinkedHashMap<>();
		for (var e : cat.entrySet()) {
			System.out.printf("  %-55s  %s%n", e.getKey(), e.getValue());
			multisetToName.putIfAbsent(e.getValue(), e.getKey());
		}
		System.out.printf("  → %d catalog schemes realise %d distinct multisets%n",
				cat.size(), multisetToName.size());
		assertThat(cat).isNotEmpty();
	}

	/**
	 * Short capped brute-force run: measures node throughput and reports the
	 * distinct multisets discovered so far, labelling catalog matches. Not
	 * exhaustive at this cap — the goal here is to confirm the machinery and
	 * gauge feasibility before a full background run.
	 */
	@Test
	@Tag("slow")
	void bruteforce_capped() throws Exception {
		Map<String, String> cat = catalogMultisets();
		Map<String, String> multisetToCatalog = new LinkedHashMap<>();
		for (var e : cat.entrySet()) multisetToCatalog.putIfAbsent(e.getValue(), e.getKey());

		Ternary2x2x2MultisetEnumerator.Config cfg = new Ternary2x2x2MultisetEnumerator.Config();
		cfg.maxMillis = 60_000; // 60s probe
		cfg.progressEveryNodes = 20_000_000L;

		System.out.println();
		System.out.println("=== Brute-force ternary ⟨2,2,2⟩ rank-7 (capped 60s) ===");
		Ternary2x2x2MultisetEnumerator.Result res = Ternary2x2x2MultisetEnumerator.enumerate(cfg);

		System.out.printf("nodes=%,d  solutions=%,d  distinct-multisets=%d  exhaustive=%b  %,d ms%n",
				res.nodes, res.solutions, res.representatives.size(), res.exhaustive, res.elapsedMillis);
		System.out.println();
		System.out.println("Distinct multisets discovered:");
		int i = 0;
		for (var e : res.solutionCounts.entrySet()) {
			i++;
			String label = multisetToCatalog.getOrDefault(e.getKey(), "(not in catalog)");
			System.out.printf("  #%-2d  ×%-10d  %-60s  [%s]%n", i, e.getValue(), e.getKey(), label);
		}
		assertThat(res.solutions).isGreaterThan(0);
	}
}
