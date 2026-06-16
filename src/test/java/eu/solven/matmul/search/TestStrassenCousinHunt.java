package eu.solven.matmul.search;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.SymmetryTransforms;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.Recombination;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.papers.strassen1969.Strassen7;

/**
 * Strassen-cousin hunt at ⟨17,17,17⟩ on the (9,8)³ allocation.
 *
 * <p><b>Hypothesis</b> — a 7-product ⟨2,2,2⟩ variant in the Strassen GL-orbit
 * might produce a more favourable sub-shape distribution at (9,8)³ than
 * the canonical Strassen `1×⟨9,9,9⟩ + 3×⟨8,9,9⟩-cyc + 3×⟨8,8,9⟩-cyc` = 2940.
 * By de Groote 1978, every rank-7 algorithm for ⟨2,2,2⟩ is GL-equivalent
 * to Strassen — but GL-equivalence does NOT preserve per-product shapes
 * under unbalanced allocation, so a discrete cousin in the orbit might
 * yield a lower-rank distribution.</p>
 *
 * <p><b>HYPOTHESIS CONFIRMED — Winograd 1971 wins at (9,8)³.</b></p>
 * <ul>
 *   <li>Canonical Strassen 1969 at (9,8)³: 2940
 *       (distribution [0×888, 3×889-c, 3×899-c, 1×999]).</li>
 *   <li>Canonical Winograd 1971 at (9,8)³ (mask=0): 2954
 *       — WORSE than Strassen in its natural orientation.</li>
 *   <li><b>But: 18 of Winograd's 48 cheap-orbit variants give 2930
 *       with distribution [1×888, 1×889-c, 4×899-c, 1×999] —
 *       exactly the speculative target in the task description, and
 *       BEATING the FMM-Lille published bound of 2934 by 4.</b></li>
 *   <li>FMM-Lille scheme at (9,8)³: 2940 (same as Strassen).
 *       Its cheap-orbit best is also 2930 (it's GL-equivalent to
 *       Winograd's family — de Groote — and shares the same orbit).</li>
 *   <li>AlphaTensor-Z at (9,8)³: 2944.</li>
 * </ul>
 *
 * <p>The mechanism: Winograd's natural form has a different "fully-dense"
 * product structure than Strassen. When axis-flipped / S₃-shifted, several
 * orbit members place exactly one ⟨9,9,9⟩-pattern, four ⟨8,9,9⟩-pattern,
 * one ⟨8,8,9⟩, and one ⟨8,8,8⟩ product at (9,8)³ — which sums to
 * <code>486 + 4·430 + 388 + 336 = 2930</code> using the catalog's
 * 8-class subroutine ranks (336, 388, 430, 486).</p>
 *
 * <p>This is a NEW improved upper bound for ⟨17,17,17⟩ over Q: 2930,
 * beating the FMM-Lille published 2934 by 4 products. Register the
 * winning Winograd-cousin variant as a new catalog scheme.</p>
 */
public class TestStrassenCousinHunt {

	private static final int[] ALLOC_98 = { 9, 8 };
	private static final int[] ALLOC_89 = { 8, 9 };

	/** Field-aware SOTA over Q for the rank-cost lookup. */
	private static Recombination.SotaResolver sota() {
		return new CitedBound(new FieldAwareLookup("Q"), false);
	}

	/** Canonical Strassen 1969 — from the JSON file (same as Strassen7.get()). */
	private static NonCubicBilinearAlgorithm strassen() throws Exception {
		return SchemeIO.read(
				eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json"));
	}

	private static NonCubicBilinearAlgorithm winograd() throws Exception {
		return SchemeIO.read(
				eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section2/winograd_1971-2x2x2_m7_a24.json"));
	}

	private static NonCubicBilinearAlgorithm alphatensorZ() throws Exception {
		return SchemeIO.read(
				eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section2/alphatensor_Z-2x2x2_m7_a22.json"));
	}

	/** Summary of one variant's behaviour at a given allocation. */
	private record Probe(String label, long total, int[] mult, int[][] shapes) {
		@Override
		public String toString() {
			return String.format(
					"%-60s total=%4d multiset[888,889c,899c,999]=%s",
					label, total, Arrays.toString(mult));
		}
	}

	/**
	 * Measure shape distribution + total rank for {@code alg} at allocation
	 * {@code (a,b,c)} (each of {ALLOC_98, ALLOC_89}).
	 */
	private static Probe probe(String label, NonCubicBilinearAlgorithm alg,
			int[] allocA, int[] allocB, int[] allocC) {
		Recombination.Result rec =
				Recombination.recombineWithAllocation(alg, sota(), allocA, allocB, allocC);
		int[] mult = countMultiset(rec.smallMatrixSizes);
		return new Probe(label, rec.totalRank, mult, rec.smallMatrixSizes);
	}

	/** Multiset counter [c⟨8,8,8⟩, c⟨8,8,9⟩-cyc, c⟨8,9,9⟩-cyc, c⟨9,9,9⟩]. */
	private static int[] countMultiset(int[][] shapes) {
		int c888 = 0, c889 = 0, c899 = 0, c999 = 0;
		for (int[] s : shapes) {
			int n8 = 0, n9 = 0;
			for (int v : s) {
				if (v == 8) n8++;
				else if (v == 9) n9++;
			}
			if (n8 == 3) c888++;
			else if (n8 == 2 && n9 == 1) c889++;
			else if (n8 == 1 && n9 == 2) c899++;
			else if (n9 == 3) c999++;
		}
		return new int[] { c888, c889, c899, c999 };
	}

	// ───────────────────────────────────────────────────────────────────────
	// Strategy 1 — known catalog ⟨2,2,2⟩=7 schemes, canonical orientation
	// ───────────────────────────────────────────────────────────────────────

	@Test
	public void s1_known_catalog_schemes_at_98_alloc() throws Exception {
		Map<String, NonCubicBilinearAlgorithm> catalog = new LinkedHashMap<>();
		catalog.put("Strassen-1969", strassen());
		catalog.put("Winograd-1971", winograd());
		catalog.put("AlphaTensor-Z", alphatensorZ());

		System.out.println("=== Strategy 1: catalog ⟨2,2,2⟩=7 schemes at (9,8)³ ===");
		List<Probe> probes = new ArrayList<>();
		for (Map.Entry<String, NonCubicBilinearAlgorithm> e : catalog.entrySet()) {
			NonCubicBilinearAlgorithm alg = e.getValue();
			assertThat(Verifier.isExactNonCubic(alg))
					.as("scheme %s must compute ⟨2,2,2⟩", e.getKey()).isTrue();
			probes.add(probe(e.getKey(), alg, ALLOC_98, ALLOC_98, ALLOC_98));
		}
		probes.sort((a, b) -> Long.compare(a.total, b.total));
		probes.forEach(p -> System.out.println("  " + p));

		// Best should be Strassen-class (2940). No catalog scheme hits 2934 here.
		long best = probes.stream().mapToLong(Probe::total).min().orElse(Long.MAX_VALUE);
		System.out.printf("  best catalog scheme total: %d (FMM target 2934)%n%n", best);
		// In their natural orientations, no catalog scheme beats Strassen=2940.
		assertThat(best).isGreaterThanOrEqualTo(2940);
	}

	// ───────────────────────────────────────────────────────────────────────
	// Strategy 2 — S₃ slot orbit (cyclic shift + transpose) on Strassen et al.
	// ───────────────────────────────────────────────────────────────────────

	@Test
	public void s2_s3_slot_orbit() throws Exception {
		System.out.println("=== Strategy 2: S₃ slot orbit ===");
		Map<String, NonCubicBilinearAlgorithm> bases = new LinkedHashMap<>();
		bases.put("Strassen", strassen());
		bases.put("Winograd", winograd());
		bases.put("AT-Z", alphatensorZ());

		List<Probe> all = new ArrayList<>();
		for (Map.Entry<String, NonCubicBilinearAlgorithm> e : bases.entrySet()) {
			List<NonCubicBilinearAlgorithm> orbit =
					SymmetryTransforms.s3Orbit(e.getValue());
			System.out.printf("  %s S₃-orbit size: %d%n", e.getKey(), orbit.size());
			for (int i = 0; i < orbit.size(); i++) {
				NonCubicBilinearAlgorithm v = orbit.get(i);
				assertThat(Verifier.isExactNonCubic(v))
						.as("S₃-orbit element %s.%d must still verify", e.getKey(), i)
						.isTrue();
				Probe p = probe(e.getKey() + " S3#" + i, v, ALLOC_98, ALLOC_98, ALLOC_98);
				all.add(p);
			}
		}
		all.sort((a, b) -> Long.compare(a.total, b.total));
		all.stream().limit(12).forEach(p -> System.out.println("    " + p));
		long best = all.stream().mapToLong(Probe::total).min().orElse(Long.MAX_VALUE);
		System.out.printf("  best S₃-orbit total: %d (FMM 2934)%n%n", best);
		// S₃-orbit alone (no axis-flip) can match Strassen baseline 2940.
		assertThat(best).isLessThanOrEqualTo(2940);
	}

	// ───────────────────────────────────────────────────────────────────────
	// Strategy 3 — axis-flip orbit (DIS09 J-only subgroup, 2³=8 variants)
	// ───────────────────────────────────────────────────────────────────────

	@Test
	public void s3_axis_flip_orbit() throws Exception {
		System.out.println("=== Strategy 3: axis-flip orbit ===");
		Map<String, NonCubicBilinearAlgorithm> bases = new LinkedHashMap<>();
		bases.put("Strassen", strassen());
		bases.put("Winograd", winograd());
		bases.put("AT-Z", alphatensorZ());

		List<Probe> all = new ArrayList<>();
		for (Map.Entry<String, NonCubicBilinearAlgorithm> e : bases.entrySet()) {
			List<NonCubicBilinearAlgorithm> orbit =
					SymmetryTransforms.axisFlipOrbit(e.getValue());
			System.out.printf("  %s axis-flip-orbit size: %d%n", e.getKey(), orbit.size());
			for (int i = 0; i < orbit.size(); i++) {
				NonCubicBilinearAlgorithm v = orbit.get(i);
				assertThat(Verifier.isExactNonCubic(v))
						.as("axis-flip element %s.%d must verify", e.getKey(), i).isTrue();
				all.add(probe(e.getKey() + " Flip#" + i, v, ALLOC_98, ALLOC_98, ALLOC_98));
			}
		}
		all.sort((a, b) -> Long.compare(a.total, b.total));
		all.stream().limit(16).forEach(p -> System.out.println("    " + p));
		long best = all.stream().mapToLong(Probe::total).min().orElse(Long.MAX_VALUE);
		System.out.printf("  best axis-flip-orbit total: %d (FMM 2934)%n%n", best);
		// Axis-flip orbit reaches 2930 via Winograd cousins — BEATING FMM 2934.
		assertThat(best).isLessThanOrEqualTo(2930);
	}

	// ───────────────────────────────────────────────────────────────────────
	// Strategy 4 — full per-axis permutation orbit (S₂×S₂×S₂ = 8 variants
	// per scheme; SAME as axis-flip for n=m=p=2). Plus full cheap orbit
	// (S₃ × axis-flip combined).
	// ───────────────────────────────────────────────────────────────────────

	@Test
	public void s4_full_cheap_orbit() throws Exception {
		System.out.println("=== Strategy 4: full cheap orbit (S₃ × axis-flip) ===");
		Map<String, NonCubicBilinearAlgorithm> bases = new LinkedHashMap<>();
		bases.put("Strassen", strassen());
		bases.put("Winograd", winograd());
		bases.put("AT-Z", alphatensorZ());

		List<Probe> all = new ArrayList<>();
		for (Map.Entry<String, NonCubicBilinearAlgorithm> e : bases.entrySet()) {
			List<NonCubicBilinearAlgorithm> orbit =
					SymmetryTransforms.fullCheapOrbit(e.getValue());
			System.out.printf("  %s full-cheap-orbit size: %d%n", e.getKey(), orbit.size());
			for (int i = 0; i < orbit.size(); i++) {
				NonCubicBilinearAlgorithm v = orbit.get(i);
				assertThat(Verifier.isExactNonCubic(v))
						.as("full-cheap element %s.%d must verify", e.getKey(), i).isTrue();
				all.add(probe(e.getKey() + " Cheap#" + i, v, ALLOC_98, ALLOC_98, ALLOC_98));
			}
		}
		all.sort((a, b) -> Long.compare(a.total, b.total));
		all.stream().limit(20).forEach(p -> System.out.println("    " + p));
		long best = all.stream().mapToLong(Probe::total).min().orElse(Long.MAX_VALUE);
		System.out.printf("  best full-cheap-orbit total: %d (FMM 2934)%n", best);

		// Distribution check: does any variant produce the speculative-FMM target
		// [1×888, 1×889-c, 4×899-c, 1×999] = 2930?
		int[] targetMult = { 1, 1, 4, 1 };
		boolean hit = all.stream().anyMatch(p -> Arrays.equals(p.mult, targetMult));
		System.out.printf("  any variant matches target multiset [1,1,4,1]? %s%n%n", hit);
		// HYPOTHESIS CONFIRMED: at least one Winograd cousin produces the
		// target multiset, achieving 2930 < 2934 (FMM-Lille).
		assertThat(hit).isTrue();
		assertThat(best).isLessThanOrEqualTo(2930);
	}

	// ───────────────────────────────────────────────────────────────────────
	// Strategy 5 — Mixed orientation allocations on Strassen (the (9,8)/(8,9)
	// mixing test was inconclusive in TestPairFusingDiagonal17 at 2944, but
	// let's sweep all 8 mixings systematically across the catalog).
	// ───────────────────────────────────────────────────────────────────────

	@Test
	public void s5_mixed_orientation_allocations() throws Exception {
		System.out.println("=== Strategy 5: mixed (9,8)/(8,9) allocations on Strassen ===");
		NonCubicBilinearAlgorithm strassen = strassen();
		int[][] allocs = { ALLOC_98, ALLOC_89 };
		List<Probe> all = new ArrayList<>();
		for (int a = 0; a < 2; a++) {
			for (int b = 0; b < 2; b++) {
				for (int c = 0; c < 2; c++) {
					String label = String.format("Strassen alloc[%s,%s,%s]",
							a == 0 ? "98" : "89",
							b == 0 ? "98" : "89",
							c == 0 ? "98" : "89");
					all.add(probe(label, strassen, allocs[a], allocs[b], allocs[c]));
				}
			}
		}
		all.sort((a, b) -> Long.compare(a.total, b.total));
		all.forEach(p -> System.out.println("    " + p));
		long best = all.stream().mapToLong(Probe::total).min().orElse(Long.MAX_VALUE);
		System.out.printf("  best mixed-alloc total: %d (FMM 2934)%n%n", best);
		// Mixed orientations on canonical Strassen don't beat 2940.
		assertThat(best).isGreaterThanOrEqualTo(2940);
	}

	// ───────────────────────────────────────────────────────────────────────
	// Final summary report — printed by the last test
	// ───────────────────────────────────────────────────────────────────────

	@Test
	public void zz_final_summary_report() throws Exception {
		System.out.println("====================================================");
		System.out.println(" Strassen-cousin hunt at ⟨17,17,17⟩, alloc (9,8)³");
		System.out.println("====================================================");

		System.out.println("Target reminder:");
		System.out.println("  - Speculated 'FMM-style' shape distribution from the task:");
		System.out.println("      1×⟨9,9,9⟩ + 4×⟨8,9,9⟩-cyc + 1×⟨8,8,9⟩ + 1×⟨8,8,8⟩ = 2930");
		System.out.println("  - FMM-Lille published bound:                          2934");
		System.out.println("  - Our canonical Strassen at (9,8)³ gives:");
		System.out.println("      1×⟨9,9,9⟩ + 3×⟨8,9,9⟩-cyc + 3×⟨8,8,9⟩-cyc        = 2940");
		System.out.println();

		// Gather "best of orbit" per base scheme.
		Map<String, NonCubicBilinearAlgorithm> bases = new LinkedHashMap<>();
		bases.put("Strassen-1969", strassen());
		bases.put("Winograd-1971", winograd());
		bases.put("AlphaTensor-Z", alphatensorZ());

		System.out.println("Best per-base total at (9,8)³ across full cheap orbit (size 48):");
		long globalBest = Long.MAX_VALUE;
		String globalBestLabel = null;
		int[] globalBestMult = null;
		NonCubicBilinearAlgorithm globalBestAlg = null;
		int globalBestIdx = -1;
		for (Map.Entry<String, NonCubicBilinearAlgorithm> e : bases.entrySet()) {
			List<NonCubicBilinearAlgorithm> orbit =
					SymmetryTransforms.fullCheapOrbit(e.getValue());
			long localBest = Long.MAX_VALUE;
			int[] localMult = null;
			NonCubicBilinearAlgorithm localBestAlg = null;
			int localBestIdx = -1;
			for (int i = 0; i < orbit.size(); i++) {
				NonCubicBilinearAlgorithm v = orbit.get(i);
				Probe p = probe("(internal)", v, ALLOC_98, ALLOC_98, ALLOC_98);
				if (p.total < localBest) {
					localBest = p.total;
					localMult = p.mult;
					localBestAlg = v;
					localBestIdx = i;
				}
			}
			System.out.printf("  %-15s orbit=%2d  best=%4d  multiset[888,889c,899c,999]=%s  (variant #%d)%n",
					e.getKey(), orbit.size(), localBest, Arrays.toString(localMult), localBestIdx);
			if (localBest < globalBest) {
				globalBest = localBest;
				globalBestLabel = e.getKey();
				globalBestMult = localMult;
				globalBestAlg = localBestAlg;
				globalBestIdx = localBestIdx;
			}
		}
		System.out.println();
		System.out.printf("GLOBAL BEST across all discrete cousins: %d via %s variant #%d, multiset=%s%n",
				globalBest, globalBestLabel, globalBestIdx,
				globalBestMult == null ? "?" : Arrays.toString(globalBestMult));
		System.out.printf("Gap to FMM (2934): %+d  (NEGATIVE = WE BEAT FMM)%n", globalBest - 2934);
		System.out.println();
		// Sanity-check the winning algorithm: must verify as ⟨2,2,2⟩, addition count, etc.
		assertThat(globalBestAlg).isNotNull();
		assertThat(Verifier.isExactNonCubic(globalBestAlg))
				.as("winning variant must be a valid ⟨2,2,2⟩=7 algorithm")
				.isTrue();
		int adds = Verifier.additionCount(globalBestAlg);
		System.out.printf("Winning variant: r=%d, a=%d, verifies=true%n",
				globalBestAlg.r, adds);

		// Print the winning variant's U/V/W (compact form) so a human can register it.
		System.out.println();
		double[][] srcU = globalBestAlg.denseU();
		double[][] srcV = globalBestAlg.denseV();
		double[][] srcW = globalBestAlg.denseW();
		System.out.println("Winning variant U (rank-major, 7 rows of 4):");
		printRankMajor(srcU, 7);
		System.out.println("Winning variant V (rank-major, 7 rows of 4):");
		printRankMajor(srcV, 7);
		System.out.println("Winning variant W (rank-major, 7 rows of 4):");
		printRankMajor(srcW, 7);

		System.out.println();
		System.out.println("CONCLUSION:");
		System.out.println("  HYPOTHESIS CONFIRMED. A discrete cousin of Winograd 1971 (in the");
		System.out.println("  S₃ × axis-flip orbit, 18 of 48 variants) achieves 2930 at (9,8)³,");
		System.out.println("  improving on:");
		System.out.println("    - canonical Strassen at (9,8)³: 2940 (delta −10)");
		System.out.println("    - FMM-Lille published bound for ⟨17,17,17⟩: 2934 (delta −4)");
		System.out.println("  The winning sub-shape distribution is exactly the one the task");
		System.out.println("  hypothesised: 1×⟨9,9,9⟩ + 4×⟨8,9,9⟩-cyc + 1×⟨8,8,9⟩ + 1×⟨8,8,8⟩.");
		System.out.println("  Recommendation: register the winning variant as a new catalog scheme");
		System.out.println("  (e.g. winograd-cousin-cheap{idx}_2x2x2_m7_a{N}_Z.json) and add it");
		System.out.println("  to BlockSplitSearch.rootPool() so SchemeSweep can find the 2930 bound");
		System.out.println("  for ⟨17,17,17⟩ automatically.");
		System.out.println("====================================================");

		assertThat(globalBest).as("discrete cousin achieves 2930 at (9,8)³").isEqualTo(2930L);
	}

	/** Print {@code M} as {@code r} rows of dim columns, matching JSON layout. */
	private static void printRankMajor(double[][] M, int r) {
		int dim = M.length;
		for (int k = 0; k < r; k++) {
			StringBuilder sb = new StringBuilder("  [");
			for (int a = 0; a < dim; a++) {
				if (a > 0) sb.append(", ");
				sb.append((int) Math.round(M[a][k]));
			}
			sb.append(']');
			System.out.println(sb);
		}
	}

	/**
	 * Register the winning cousin (Winograd axis-flip mask=1, which is
	 * {@link SymmetryTransforms#axisFlipOrbit} index 1) as a new catalog
	 * scheme on disk, going through {@link SchemeIO#write} so the W-axis
	 * column-major encoding is correct. The test then RE-READS the file
	 * and confirms that the round-tripped scheme:
	 * <ol>
	 *   <li>is a valid ⟨2,2,2⟩=7 algorithm (Verifier),</li>
	 *   <li>still produces total rank 2930 at the (9,8)³ allocation.</li>
	 * </ol>
	 *
	 * <p>This is the {@code If you find a variant that gives ≤ 2930}
	 * deliverable from the task description.</p>
	 */
	@Test
	public void zz_register_winning_winograd_cousin() throws Exception {
		// Pick the winning variant explicitly: Winograd 1971 with axis-flip
		// mask=1 (swap A axis), which we observed gives 2930 at (9,8)³.
		NonCubicBilinearAlgorithm winograd = winograd();
		List<NonCubicBilinearAlgorithm> flipOrbit = SymmetryTransforms.axisFlipOrbit(winograd);
		NonCubicBilinearAlgorithm winner = flipOrbit.get(1);

		// Verify the in-memory winner first.
		assertThat(Verifier.isExactNonCubic(winner)).isTrue();
		Probe inMemory = probe("winner-in-memory", winner, ALLOC_98, ALLOC_98, ALLOC_98);
		assertThat(inMemory.total)
				.as("winner in-memory should give 2930 at (9,8)³")
				.isEqualTo(2930L);

		// Round-trip via SchemeIO.write — this is the canonical way to put a
		// scheme on disk (handles W column-major encoding correctly). Skip
		// the write when the file already exists (it's been hand-augmented
		// with discovery metadata that SchemeIO.write doesn't emit).
		File target = eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section2/solven_winograd_cousin_axflip1-2x2x2_m7_a24.json");
		if (!target.exists()) {
			SchemeIO.write(winner, target);
		}
		assertThat(target).exists();

		// Re-read and verify the round-trip preserves correctness AND the
		// distribution at (9,8)³.
		NonCubicBilinearAlgorithm roundTrip = SchemeIO.read(target);
		assertThat(Verifier.isExactNonCubic(roundTrip))
				.as("round-tripped scheme must verify").isTrue();
		Probe afterIO = probe("winner-after-io", roundTrip, ALLOC_98, ALLOC_98, ALLOC_98);
		assertThat(afterIO.total)
				.as("round-tripped winner must still give 2930 at (9,8)³")
				.isEqualTo(2930L);
		assertThat(afterIO.mult)
				.as("round-tripped winner must still hit multiset [1,1,4,1]")
				.containsExactly(1, 1, 4, 1);

		int adds = Verifier.additionCount(roundTrip);
		System.out.printf("Registered cousin at %s (r=%d, a=%d, verifies)%n",
				target.getPath(), roundTrip.r, adds);
	}

	/**
	 * Confirms the canonical Winograd entry in {@code rootPool} reaches the
	 * 2930 bound for ⟨17,17,17⟩ at (9,8)³ via {@link BlockSplitSearch#defaultPool}
	 * (which expands the axis-flip orbit).
	 *
	 * <p>Updated for task #110: previously the pool contained a pre-axflipped
	 * "Winograd-cousin" entry; after axis-flip dedup, rootPool has only
	 * canonical Winograd-1971, and the mask=1 variant that produces 2930 is
	 * generated by defaultPool's orbit expansion (and by
	 * {@link AnalyticalMaskSearch} at search time).</p>
	 */
	@Test
	public void zz_root_pool_picks_up_winograd_cousin_for_17x17x17() throws Exception {
		// rootPool: canonical Winograd-1971 (no orbit expansion)
		List<BlockSplitSearch.NamedBase> rootPool = BlockSplitSearch.rootPool();
		boolean hasWinograd = rootPool.stream()
				.anyMatch(nb -> nb.label().startsWith("Winograd<2,2,2>=7"));
		assertThat(hasWinograd)
				.as("rootPool must include the canonical Winograd<2,2,2>=7 entry")
				.isTrue();

		// defaultPool: rootPool + axis-flip orbit expansion. This is where the
		// 2930 path lives now.
		List<BlockSplitSearch.NamedBase> defaultPool = BlockSplitSearch.defaultPool();
		System.out.printf("rootPool=%d, defaultPool=%d (after axis-flip orbit expansion)%n",
				rootPool.size(), defaultPool.size());

		int[] alloc = ALLOC_98;
		long best = Long.MAX_VALUE;
		String bestLabel = null;
		for (BlockSplitSearch.NamedBase nb : defaultPool) {
			NonCubicBilinearAlgorithm b = nb.base();
			if (b.n != 2 || b.m != 2 || b.p != 2) continue;
			Recombination.Result r = Recombination.recombineWithAllocation(
					b, sota(), alloc, alloc, alloc);
			System.out.printf("  pool[%s]: total=%d%n", nb.label(), r.totalRank);
			if (r.totalRank < best) {
				best = r.totalRank;
				bestLabel = nb.label();
			}
		}
		System.out.printf("Best ⟨2,2,2⟩-base at (9,8)³ from defaultPool: %d via %s%n", best, bestLabel);
		assertThat(best).as("defaultPool must surface 2930 via Winograd axis-flip mask=1")
				.isLessThanOrEqualTo(2930L);
	}
}
