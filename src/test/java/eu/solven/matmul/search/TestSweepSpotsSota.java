package eu.solven.matmul.search;

import eu.solven.matmul.recombination.BlockSplitSearch;

import static org.assertj.core.api.Assertions.assertThat;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Optional;

import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;

import eu.solven.matmul.verifiers.Verifier;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.papers.dis2009.PanTrilinearAggregation;
import tools.jackson.databind.JsonNode;

/**
 * Fast, parameter-specific <em>regression guards</em> on the search/materialise
 * pipeline: with the real catalog, the recursive materialiser must still spot the
 * known SOTA (or better) for a spread of shapes that each exercise a distinct
 * mechanism. This is the cheap counterpart to the full {@code SchemeSweep} /
 * {@code VerifyAllSchemes} runs — it builds the lookup once and probes a handful
 * of shapes, so it runs in seconds and fails loudly when an engine silently
 * regresses (e.g. the 2026-06-10 empty-{@code extendedPool} bug, or the
 * bud-ordering bug that hid ⟨8,9,9⟩=430).
 *
 * <p>Each bound is an <em>upper</em> bound (≤): the pipeline may legitimately find
 * something better, but never worse. Keep the shape list small + low-dim so the
 * suite stays fast.</p>
 */
public class TestSweepSpotsSota {

	private static FieldAwareLookup lookup;
	private static RecursiveMaterialiser mat;

	@BeforeAll
	static void setUp() {
		lookup = new FieldAwareLookup("Q");
		List<BlockSplitSearch.NamedBase> pool = BlockSplitSearch.defaultPool();
		RecursiveClosureSota sota = new RecursiveClosureSota(lookup, pool, true, true);
		// dry-run (no disk writes); composes from the real catalog.
		mat = new RecursiveMaterialiser(lookup, pool, sota, null, false, true);
	}

	@ParameterizedTest(name = "⟨{0},{1},{2}⟩ ≤ {3}")
	@CsvSource({
			// shape         SOTA   mechanism exercised
			"2, 2, 2,    7",   // Strassen — direct disk hit
			"3, 3, 3,   23",   // Laderman — direct disk hit
			"4, 4, 4,   49",   // Kronecker ⟨2,2,2⟩⊗⟨2,2,2⟩
			"3, 7, 8,  126",   // concat ⟨3,7,4⟩ +ₚ ⟨3,7,4⟩ (orientation-aware)
			"7, 7, 7,  250",   // Strassen-recursion recombination
			"6, 8, 9,  296",   // serendipitous bud-product (bud-ordering fix)
			// ⟨8,9,9⟩=430: FILL-mode disk-presence (the derived 430 is on disk, persisted via
			// MaterialiseSerendipitousWins after the 2026-06-23 σ-base-selection fix). The
			// COMPUTE path is guarded separately by compute_pipeline_reaches_8x9x9_430.
			"8, 9, 9,  430",   // serendipitous bud-product (σ-aware V-bud base)
			"4, 8, 12, 272",   // serendipitous
			"8, 8, 12, 504",   // serendipitous
			"4, 20, 14, 736",  // serendipitous ⟨2,4,7⟩-base band-20 win (fmm-react 2026-07-06, was 755)
			"9, 9, 21, 1058",  // serendipitous ⟨3,3,7⟩⊗⟨3,3,3⟩ band-21 win (fmm-react 2026-07-06)
	})
	public void materialise_spots_sota(int n, int m, int p, int sota) {
		Optional<RecursiveMaterialiser.Result> r = mat.materialise(n, m, p);
		assertThat(r).as("⟨%d,%d,%d⟩ should resolve", n, m, p).isPresent();
		assertThat(r.get().alg().r)
				.as("⟨%d,%d,%d⟩ rank must be ≤ SOTA %d (regression if higher)", n, m, p, sota)
				.isLessThanOrEqualTo(sota);
		assertThat(Verifier.passesRandomMatmulSpotCheck(r.get().alg()))
				.as("⟨%d,%d,%d⟩ result must verify", n, m, p).isTrue();
	}

	/**
	 * COMPUTE-path guard for the 2026-06-23 σ-base-selection fix (distinct from the
	 * fill-mode {@link #materialise_spots_sota}, which only checks the DISK best). An
	 * IMPROVE-mode materialiser actually composes ⟨8,9,9⟩ and must reach the
	 * serendipitous SOTA 430 = (⟨4,3,3⟩=29−3)⊗⟨2,3,3⟩+⟨6,3,3⟩=40. That requires
	 * {@code trySerendipitous} to feed the size-3 V-bud base (budScore 4) to
	 * {@code SerendipitousSearch.bestFor}, not just the budScore-MAX sibling (11, a
	 * U-bud with σ_V=0 here). A regression to the count-based picker returns 432.
	 */
	@Test
	public void compute_pipeline_reaches_8x9x9_430() {
		List<BlockSplitSearch.NamedBase> pool = BlockSplitSearch.defaultPool();
		RecursiveClosureSota sota = new RecursiveClosureSota(lookup, pool, true, true);
		// improveExisting=true → composes instead of returning the dense 432 import; no disk write.
		RecursiveMaterialiser improver =
				new RecursiveMaterialiser(lookup, pool, sota, null, false, true, true);
		// Restrict to the serendipitous strategy: it's the one under test, and skipping the
		// recombination B&B keeps the guard fast (~seconds, not ~25s).
		improver.setStrategies(java.util.Set.of(RecursiveMaterialiser.STRAT_SERENDIPITOUS));
		Optional<RecursiveMaterialiser.Result> r = improver.materialise(8, 9, 9);
		assertThat(r).as("⟨8,9,9⟩ should resolve").isPresent();
		assertThat(r.get().alg().r)
				.as("compose() must reach the serendipitous SOTA 430 via the σ-selected size-3 "
						+ "V-bud base (regression to the budScore picker → 432)")
				.isLessThanOrEqualTo(430);
		assertThat(Verifier.passesRandomMatmulSpotCheck(r.get().alg()))
				.as("⟨8,9,9⟩=430 result must verify").isTrue();
	}

	/**
	 * Regression guard for the 2026-06-23 projection-parent orientation-pinning fix
	 * (project_projection_parent_orientation_not_pinned). ⟨19,19,20⟩ = Project(⟨20,19,20⟩);
	 * its parent ⟨19,20,20⟩ has two equal 20-axes, and pinning it as an ORIENTED
	 * {@code 20x19x20@hash} let replay re-{@code orientAs} ambiguously to a worse-projecting
	 * axis → predict 4154 / build 4237 → fatal {@code assertRebuildNotWorse}. The fix pins the
	 * NATIVE {@code 19x20x20@hash} + an exact-perm {@code Transpose} (so the lineage is
	 * {@code Project(Transpose(19x20x20, "ABC->CAB"), …)}), making predict==build. Materialise
	 * must NOT throw and must reach master's 4154. Projection-only strategy keeps it fast.
	 */
	@Test
	public void projection_parent_orientation_pinned_19x19x20() {
		List<BlockSplitSearch.NamedBase> pool = BlockSplitSearch.defaultPool();
		RecursiveClosureSota sota = new RecursiveClosureSota(lookup, pool, true, true);
		// derive-best (last arg) so a TIE with an already-on-disk 4154 still returns the composed
		// result (improve-mode would return empty on a non-strict-improvement, masking the build).
		RecursiveMaterialiser improver =
				new RecursiveMaterialiser(lookup, pool, sota, null, false, false, true, true);
		improver.setStrategies(java.util.Set.of(RecursiveMaterialiser.STRAT_PROJECTION));
		Optional<RecursiveMaterialiser.Result> r = improver.materialise(19, 19, 20);
		assertThat(r).as("⟨19,19,20⟩ must resolve (no projection-divergence throw)").isPresent();
		assertThat(r.get().alg().r)
				.as("⟨19,19,20⟩ projection must build at its predicted rank (≤4154 = master), not "
						+ "diverge to 4237 — the parent orientation must be pinned bit-exactly")
				.isLessThanOrEqualTo(4154);
		assertThat(Verifier.passesRandomMatmulSpotCheck(r.get().alg()))
				.as("⟨19,19,20⟩=4154 result must verify").isTrue();
	}

	/**
	 * Fast catalog invariant guarding the 2026-06-22 dis09-cube phantom class: no
	 * <em>derived</em> ⟨n,n,n⟩ file may duplicate the Pan-TA formula rank
	 * {@code cubicBound(n)} under a NON-formula lineage. Such a file is the
	 * {@code CORRUPT_RANK} phantom that was purged: it stamped the right rank
	 * (4340/5566/7000/8658 == its {@code known/…dis09_Q…} twin) but a bogus
	 * {@code Project(⟨n,n,n+1⟩)} lineage that only replays to ~4378. With the
	 * honest twin already pricing {@code findRank} at 4340, the duplicate did not
	 * change the score — it silently mis-led {@code resolveSubScheme} into BUILDING
	 * the ⟨n,n,n⟩ block as the worse projection, so the whole ⟨2k,2k,2k+2⟩ family
	 * diverged (⟨20,20,22⟩ evaluated 4950 / built 4988). The honest cube lives in
	 * {@code known/…dis09_Q…} as a {@code DIS09Lemma4(n)} atom; any derived cube at
	 * the same rank with a non-formula lineage is the phantom. Pure file-scan — no
	 * {@code materialise} (which is unbounded and minutes-to-hours on these shapes)
	 * — so it stays in milliseconds, per the fast-guard rule.
	 */
	@Test
	public void no_phantom_dis09_cube_duplicates() throws Exception {
		Path derived = Path.of("src/main/resources/schemes/derived");
		List<String> offenders = new ArrayList<>();
		try (var paths = Files.walk(derived)) {
			for (Path p : (Iterable<Path>) paths.filter(Files::isRegularFile)
					.filter(f -> f.getFileName().toString().matches("^(\\d+)x\\1x\\1-.*\\.json"))::iterator) {
				JsonNode d = SchemeIO.parseJson(p.toFile());
				if (!d.has("m") || !d.has("n")) continue;
				int n = d.get("n").get(0).asInt();
				int m = d.get("m").asInt();
				String lineage = d.has("lineage_str") ? d.get("lineage_str").asText() : "";
				if (m == PanTrilinearAggregation.cubicBound(n) && !lineage.contains("DIS09Lemma4")) {
					offenders.add(derived.relativize(p) + "  (m=" + m + "==cubicBound(" + n
							+ "), lineage='" + lineage + "')");
				}
			}
		}
		assertThat(offenders)
				.as("derived ⟨n,n,n⟩ cube(s) claim the Pan-TA rank cubicBound(n) via a non-formula "
						+ "lineage — the purged CORRUPT_RANK phantom is back (the buildable cube is the "
						+ "known/…dis09_Q… DIS09Lemma4(n) atom; a derived twin at the same rank is a phantom)")
				.isEmpty();
	}

	/**
	 * ⟨17,17,17⟩ has been a long-contested, hard shape: the plain search only
	 * reaches the <b>2940</b> floor, while the catalog holds sub-2940 results
	 * (FMM 2934, LRP/derived 2930) as <b>maxDim&gt;16 lineage-only stubs</b>.
	 * {@code materialise()} deliberately skips those stubs (it returns 2940 here),
	 * so we guard via {@code findRank}, which is stub-inclusive — checking a stub
	 * in a unit test is fine (per the user). The bound 2934 sits below the 2940
	 * search floor, so it fails loudly if the sub-2940 import/derivation is ever
	 * lost (e.g. a folder reorg or over-eager cleanup dropping the stub).
	 */
	@Test
	public void retains_hard_won_17x17x17_below_search_floor() {
		assertThat(lookup.findRank(17, 17, 17))
				.as("⟨17,17,17⟩ must retain a sub-2940 result (2930/2934-class stub); "
						+ "2940 would mean the hard-won import/derivation was lost")
				.isLessThanOrEqualTo(2934);
	}

	/**
	 * Disk-presence guards for the fmm-gap 2026-07-07 dim-7-outer-base wins
	 * (`Recombination(base=⟨3,4,7⟩:63 DPS, deficient block)`). The default pool
	 * caps non-cubic bases at maxBaseDim=5, so the COMPUTE pipeline cannot
	 * re-derive these — they exist only as exact-verified on-disk stubs reached
	 * via {@code SchemeSweep --base=3x4x7}. Losing the stubs (folder reorg,
	 * over-eager purge) would silently regress ⟨12,16,27⟩ to 2988+ (FMM: 2984)
	 * and ⟨11,16,28⟩ to 2925.
	 */
	@Test
	public void retains_fmm_gap_dim7_base_wins() {
		assertThat(lookup.findRank(12, 16, 27))
				.as("⟨12,16,27⟩ must retain the ⟨3,4,7⟩-outer-base 2964 stub (beats FMM 2984)")
				.isLessThanOrEqualTo(2964);
		assertThat(lookup.findRank(11, 16, 28))
				.as("⟨11,16,28⟩ must retain the ⟨3,4,7⟩-outer-base 2894 stub (ties FMM)")
				.isLessThanOrEqualTo(2894);
		// fmm-gap 2026-07-07 (2nd run): the ⟨28,29,31⟩ chain. ⟨3,25,28⟩=1520 is the
		// ⟨2,5,7⟩-outer deficient-A recomb (3=2+1) — another dim-7 base the default
		// pool cannot re-derive; ⟨28,29,31⟩=13091 concat-cascades from it.
		assertThat(lookup.findRank(3, 25, 28))
				.as("⟨3,25,28⟩ must retain the ⟨2,5,7⟩-outer-base 1520 stub (ties FMM)")
				.isLessThanOrEqualTo(1520);
		assertThat(lookup.findRank(28, 29, 31))
				.as("⟨28,29,31⟩ must retain the 13091 chain result (ties FMM)")
				.isLessThanOrEqualTo(13091);
		// fmm-gap 2026-07-07 (3rd run): ⟨3,4,6⟩=54 outer (dim-6, also above the
		// maxBaseDim=5 pool cap) with BOTH A and C deficient (14=[5,5,4], 29=[5×5,4]).
		assertThat(lookup.findRank(14, 28, 29))
				.as("⟨14,28,29⟩ must retain the ⟨3,4,6⟩-outer-base 6494 stub (beats FMM 6498)")
				.isLessThanOrEqualTo(6494);
		// fmm-gap 2026-07-07 (4th run): ⟨2,5,6⟩=47 outer, deficient A-split 3=2+1 with
		// 17 of 47 products isolating the width-1 block (thin-A ⟨3,·,·⟩ band opener).
		assertThat(lookup.findRank(3, 20, 30))
				.as("⟨3,20,30⟩ must retain the ⟨2,5,6⟩-outer-base 1300 stub (ties FMM)")
				.isLessThanOrEqualTo(1300);
		// fmm-gap 2026-07-07 (5th run): ⟨2,4,4⟩=26 outer with the two-axis-uneven alloc
		// [3,4 | 3,4,3,4 | 6,6,6,6] — in-pool base, but the full thorough pool starves
		// the per-base alloc budget; only a --baseFilter-concentrated run finds it.
		assertThat(lookup.findRank(7, 14, 24))
				.as("⟨7,14,24⟩ must retain the ⟨2,4,4⟩-outer-base 1514 stub (ties FMM)")
				.isLessThanOrEqualTo(1514);
		// fmm-gap 2026-07-08: serendipitous ⟨5,7,7⟩:176 (Kauers–Wood, bud-bases
		// import — size-5 U-bud + size-4 V-bud) ⊗ ⟨4,4,4⟩:48. Only reachable once
		// the serendipitous search resolves STUB fusion targets (the ⟨4,4,20⟩=230
		// ConcatCols stub) — see TestSerendipitousStubInner for the mechanism guard.
		assertThat(lookup.findRank(20, 28, 28))
				.as("⟨20,28,28⟩ must retain the serendipitous 8434 stub (beats FMM 8438)")
				.isLessThanOrEqualTo(8434);
		assertThat(lookup.findRank(21, 25, 28))
				.as("⟨21,25,28⟩ must retain the serendipitous 8118 collateral")
				.isLessThanOrEqualTo(8118);
		// fmm-gap 2026-07-08 (2nd run): ⟨2,3,3⟩:15 grid [6,8|9,9,9|9,9,9] — only
		// reachable since the AllocationOptimizer coordinate-descent seeding (the
		// stagnation cap starved single-axis-unbalanced allocs on ≥27-wide axes;
		// see TestAllocationOrdering / TestAllocationOptimizer for the mechanism).
		assertThat(lookup.findRank(14, 27, 27))
				.as("⟨14,27,27⟩ must retain the [6,8|9×3|9×3] grid 5862 stub (ties FMM)")
				.isLessThanOrEqualTo(5862);
		assertThat(lookup.findRank(17, 22, 29))
				.as("⟨17,22,29⟩ must retain the seeding-unlocked 6125 collateral")
				.isLessThanOrEqualTo(6125);
		// fmm-gap 2026-07-08 (3rd run): concat ⟨7,14,2⟩:152 + ⟨7,14,30⟩:1865 — only
		// buildable since diskBest resolves stub-only leaves by replay (the composed
		// build threw "missing sub-algorithm" on the ⟨7,14,30⟩ stub; see
		// TestDiskBestStubLeaf for the mechanism guard).
		assertThat(lookup.findRank(7, 14, 32))
				.as("⟨7,14,32⟩ must retain the stub-leaf concat 2017 (ties FMM)")
				.isLessThanOrEqualTo(2017);
		// fmm-gap 2026-07-09: ⟨7,7,5⟩-oriented Kauers–Wood base ⊗ ⟨3,4,6⟩ — only
		// reachable once budBasesAt offers ALL S₃ orientations on dims-repeat
		// shapes (see TestSerendipitousStubInner.ambiguous_orientation_*).
		assertThat(lookup.findRank(21, 28, 30))
				.as("⟨21,28,30⟩ must retain the orientation-unmasked 9473 (ties FMM index)")
				.isLessThanOrEqualTo(9473);
		assertThat(lookup.findRank(4, 18, 30))
				.as("⟨4,18,30⟩ must retain the 1394 collateral (beats FMM's artifact)")
				.isLessThanOrEqualTo(1394);
		// fmm-gap 2026-07-09 (index-recipe decode): most "absorbing-class" rows were
		// plain recombinations over CAP-EXCLUDED dim-6/7 bases (leafCount = base
		// rank!) — ⟨3,5,6⟩:68, ⟨5,5,7⟩:127, ⟨5,6,7⟩:150 via --base, plus FMM's
		// support-rich ⟨2,5,5⟩:40 rep (bud-bases import) for ⟨3,29,29⟩.
		assertThat(lookup.findRank(11, 20, 23))
				.as("⟨11,20,23⟩ must retain the ⟨3,5,6⟩-base 3035 (BEATS FMM 3039)")
				.isLessThanOrEqualTo(3035);
		assertThat(lookup.findRank(20, 20, 27))
				.as("⟨20,20,27⟩ must retain the ⟨5,5,7⟩-base 5986 (BEATS FMM 6006)")
				.isLessThanOrEqualTo(5986);
		assertThat(lookup.findRank(20, 23, 28))
				.as("⟨20,23,28⟩ must retain the ⟨5,6,7⟩-base 7100 (ties FMM)")
				.isLessThanOrEqualTo(7100);
		assertThat(lookup.findRank(3, 29, 29))
				.as("⟨3,29,29⟩ must retain the support-rep ⟨2,5,5⟩ recomb 1840 (ties FMM)")
				.isLessThanOrEqualTo(1840);
		// fmm-gap 2026-07-09 (leaf-level Pan pairing): two same-shape cubic leaves
		// of a ⟨2,2,2⟩:7 recombination fused via PanPairProduct — leaves are
		// formally independent bilinear problems, so pairing is base-agnostic
		// (RecombinationWithPairN; the KMW-2026 "merged recursive calls" device).
		assertThat(lookup.findRank(20, 23, 23))
				.as("⟨20,23,23⟩ must retain the leaf-paired 5906 (ties FMM; was the −39 gap)")
				.isLessThanOrEqualTo(5906);
		assertThat(lookup.findRank(19, 19, 22))
				.as("⟨19,19,22⟩ must retain the leaf-paired 4536 (ties FMM)")
				.isLessThanOrEqualTo(4536);
		assertThat(lookup.findRank(9, 11, 22))
				.as("⟨9,11,22⟩ must retain the support-rep ⟨3,4,4⟩ recomb 1374 (ties FMM)")
				.isLessThanOrEqualTo(1374);
	}

	/**
	 * Disk-presence guard for the six HK task-#9 shapes (issue #7, Adriano 2026,
	 * {@code schemes/known/*-adriano_2026-*}): the only {@code g = gcd(n, p/2) ≥ 6}
	 * shapes in the swept range, at the EXACT Hopcroft–Kerr formula
	 * {@code ⌈(3pn+max(p,n))/2⌉}. Our own emitter ({@code HopcroftKerr2bcAsymmetric},
	 * {@code schemes/constructed/}) provably cannot reach the formula there (+1..+3),
	 * and FMM-Lille lists the rank index-only — so nothing in the COMPUTE pipeline can
	 * re-derive these; losing the imported files would silently regress all six back
	 * to +1..+3 with no other symptom. Assert ≤ (a genuine improvement must never
	 * break this).
	 */
	/**
	 * Disk-presence guard for the issue-#8 family (2026-09-30): the 22 FMM-Lille
	 * "recipe" ranks Marcos Adriano executed as explicit schemes, RE-DERIVED here as
	 * recombination stubs — outer HK ⟨2,4,4⟩=26 over unequal blocks for 14 (the
	 * Perminov-ZT / AlphaTensor-Z root reps registered in {@code rootPool} for this),
	 * Strassen over the new Perminov ⟨7,12,16⟩=876 / ⟨7,15,16⟩=1081 pieces for 8.
	 * Thirteen of the fourteen land strictly BELOW the recipe (⟨10,19,31⟩ 3492 vs
	 * FMM 3532, ⟨10,23,23⟩ 3146 vs 3183, …). Bounds are OUR ranks (≤): losing the
	 * stubs, the ⟨2,4,4⟩ root reps, or the 876/1081 pieces regresses silently by
	 * +2…+40 with no other symptom.
	 */
	@Test
	public void retains_issue8_lille_recipe_ranks() {
		int[][] rows = {
				{ 6, 14, 25, 1322 }, { 10, 19, 31, 3492 }, { 10, 22, 25, 3288 }, { 10, 22, 29, 3772 },
				{ 10, 23, 23, 3146 }, { 10, 23, 26, 3534 }, { 10, 23, 27, 3648 }, { 10, 23, 30, 4036 },
				{ 10, 23, 31, 4150 }, { 10, 25, 26, 3850 }, { 10, 26, 26, 4014 }, { 10, 26, 27, 4144 },
				{ 10, 26, 29, 4482 }, { 10, 27, 27, 4276 },
				{ 13, 23, 31, 5396 }, { 13, 23, 32, 5552 }, { 13, 24, 31, 5540 }, { 13, 25, 32, 6008 },
				{ 13, 29, 32, 6910 }, { 13, 30, 31, 6913 }, { 13, 31, 32, 7322 }, { 15, 31, 32, 8185 } };
		for (int[] r : rows) {
			assertThat(lookup.findRank(r[0], r[1], r[2]))
					.as("⟨%d,%d,%d⟩ must retain the issue-#8 re-derivation at %d (FMM recipe rank or below)",
							r[0], r[1], r[2], r[3])
					.isLessThanOrEqualTo(r[3]);
		}
	}

	/**
	 * COMPUTE-path guard for the ⟨2,4,4⟩=26 root registration (issue #8): the
	 * {@code defaultPool()} (= rootPool + axis-flips) must carry the ⟨2,4,4⟩ reps and
	 * an IMPROVE-mode materialiser restricted to them must compose ⟨6,14,25⟩ at the
	 * FMM recipe rank 1322 from disk leaves (⟨3,4,7⟩, ⟨3,3,6⟩, …) — the AlphaTensor-Z
	 * rep's support does it; the hk71 rep alone gives 1324 (the old catalog value).
	 * A pool that silently drops the reps (or keeps only one) regresses to ≥ 1324.
	 * The pool is filtered to the ⟨2,4,4⟩ entries so the guard stays fast (seconds).
	 */
	@Test
	public void compute_pipeline_reaches_6x14x25_1322_via_244_root() {
		List<BlockSplitSearch.NamedBase> pool = BlockSplitSearch.defaultPool().stream()
				.filter(nb -> nb.label().contains("2,4,4")).toList();
		assertThat(pool).as("defaultPool must carry the ⟨2,4,4⟩ root reps").isNotEmpty();
		eu.solven.matmul.recombination.Recombination.SotaResolver diskSota = (a, b, c) -> {
			if (a == 0 || b == 0 || c == 0) return 0;
			if (a == 1) return b * c;
			if (b == 1) return a * c;
			if (c == 1) return a * b;
			return lookup.findRank(a, b, c);
		};
		// improveExisting=true → composes instead of returning the on-disk stub; deriveBest=true
		// → do NOT prune at the disk incumbent (the 1322 stub IS on disk, so a strict-improvement
		// bound would prune the very derivation under test); no write (writeRoot=null).
		RecursiveMaterialiser improver =
				new RecursiveMaterialiser(lookup, pool, diskSota, null, false, false, true, true);
		improver.setStrategies(java.util.Set.of(RecursiveMaterialiser.STRAT_RECOMBINATION));
		Optional<RecursiveMaterialiser.Result> r = improver.materialise(6, 14, 25);
		assertThat(r).as("⟨6,14,25⟩ should resolve").isPresent();
		assertThat(r.get().alg().r)
				.as("⟨2,4,4⟩-root recombination must reach the FMM recipe rank 1322 (hk71 rep alone: 1324)")
				.isLessThanOrEqualTo(1322);
		assertThat(Verifier.passesRandomMatmulSpotCheck(r.get().alg()))
				.as("⟨6,14,25⟩ result must verify").isTrue();
	}

	/**
	 * COMPUTE-path guard for DEGENERATE serendipitous products (2026-09-30): most of
	 * Perminov's serendipitous 17–32 band is {@code base ⊗ˢ ⟨1,b,c⟩} — the second
	 * factor has a unit axis, so the base's buds on that axis fuse into ⟨k,b,c⟩ blocks
	 * (⟨13,20,21⟩ = ⟨13,4,3⟩:123 ⊗ˢ ⟨1,5,7⟩ = 3165, catalog was 3291). Prediction priced
	 * the unit-axis inner (findRank → naive) and the replayer resolved it, but the
	 * materialiser's BUILD resolver returned empty for ⟨1,5,7⟩ (no catalog file — it is
	 * the naive scheme), so every such candidate was dropped as "unbuildable" and 621
	 * formats sat above Perminov's digest. Needs the bud-rich base
	 * ({@code bud-bases/section13/3x4x13-r123-perminov_serbase_*}) on disk too — the
	 * importer used to drop it under its (shape, rank) idempotence key.
	 * {@code deriveBest} so the on-disk 3165 stub does not prune the derivation.
	 */
	@Test
	public void compute_pipeline_reaches_13x20x21_3165_via_unit_axis_inner() {
		List<BlockSplitSearch.NamedBase> pool = BlockSplitSearch.defaultPool();
		eu.solven.matmul.recombination.Recombination.SotaResolver diskSota = (a, b, c) -> {
			if (a == 0 || b == 0 || c == 0) return 0;
			if (a == 1) return b * c;
			if (b == 1) return a * c;
			if (c == 1) return a * b;
			return lookup.findRank(a, b, c);
		};
		RecursiveMaterialiser improver =
				new RecursiveMaterialiser(lookup, pool, diskSota, null, false, false, true, true);
		improver.setStrategies(java.util.Set.of(RecursiveMaterialiser.STRAT_SERENDIPITOUS));
		Optional<RecursiveMaterialiser.Result> r = improver.materialise(13, 20, 21);
		assertThat(r).as("⟨13,20,21⟩ should resolve via base ⊗ˢ ⟨1,5,7⟩").isPresent();
		assertThat(r.get().alg().r)
				.as("degenerate serendipitous product must reach Perminov's 3165 (unit-axis inner dropped → 3291)")
				.isLessThanOrEqualTo(3165);
		assertThat(Verifier.passesRandomMatmulSpotCheck(r.get().alg()))
				.as("⟨13,20,21⟩ result must verify").isTrue();
	}

	/**
	 * Disk-presence guard for the Khoruzhii–Serafin–Gelß–Pokutta 2026 LITA cubes
	 * (REFERENCES [82]; {@code known/section{N}/{N}x{N}x{N}-r{R}-khoruzhii_2026-*},
	 * imported 2026-09-30 via {@code ImportKhoruzhiiLita}): the ⟨N,N,N⟩ ranks FMM-Lille's
	 * 2026-09 index cites, every one below our June-LITA {@code TA_lita} stubs by
	 * 42 (13³) … 1086 (31³). Our emitter cannot re-derive them (it ports the older
	 * generator), and every ⟨≤N⟩ family member's projection stub hangs off the cube —
	 * losing a cube file regresses the whole family silently. Assert ≤ (a better
	 * cube must never break this).
	 */
	@Test
	public void retains_kgp_lita_cubes() {
		int[][] rows = { { 13, 1379 }, { 14, 1594 }, { 15, 1977 }, { 16, 2237 }, { 17, 2723 }, { 18, 3032 },
				{ 19, 3633 }, { 20, 3995 }, { 21, 4723 }, { 22, 5142 }, { 23, 6009 }, { 24, 6489 }, { 25, 7507 },
				{ 26, 8052 }, { 27, 9233 }, { 28, 9847 }, { 29, 11203 }, { 30, 11890 }, { 31, 13433 }, { 32, 14197 } };
		for (int[] r : rows) {
			assertThat(lookup.findRank(r[0], r[0], r[0]))
					.as("⟨%d,%d,%d⟩ must retain the KGP-2026 LITA cube at %d (FMM index rank)", r[0], r[0], r[0], r[1])
					.isLessThanOrEqualTo(r[1]);
		}
	}

	@Test
	public void retains_issue7_hk_task9_formula_schemes() {
		int[][] shapes = { { 2, 12, 18 }, { 2, 14, 21 }, { 2, 16, 24 }, { 2, 18, 27 }, { 2, 20, 30 }, { 2, 24, 30 } };
		for (int[] s : shapes) {
			int p = s[1], n = s[2];
			int formula = (int) Math.ceil((3.0 * p * n + Math.max(p, n)) / 2.0);
			assertThat(lookup.findRank(2, p, n))
					.as("⟨2,%d,%d⟩ must retain the issue-#7 scheme at the HK formula %d "
							+ "(the emitter's constructed/ file sits above it)", p, n, formula)
					.isLessThanOrEqualTo(formula);
		}
	}

	/**
	 * The extended template pool must see the whole catalog tree. A
	 * {@code listFiles("section*")} on the schemes root (pre-2026-06-10 bug)
	 * silently returned an empty pool after the known/derived/curated split,
	 * crippling the search; this guards against that regression.
	 */
	/**
	 * Recombination pool must keep CONTENT-distinct bases at the same (shape, rank), not
	 * dedup to one. Two different ⟨2,4,4⟩=26 schemes (hk71 vs alphatensor_Z) tile a target
	 * DIFFERENTLY — ⟨5,20,26⟩ reaches 1700 via the alphatensor_Z one but only 1716 via hk71.
	 * The 2026-06-23 fix changed `extendedPool`/`buildPool` dedup from `shape:r` to
	 * `shape:r:contentHash`; a regression to shape:r dedup loses the better base and reopens
	 * the residual large-unbalanced master-regressions. Guards ≥2 distinct ⟨2,4,4⟩=26 schemes.
	 */
	@Test
	public void pool_keeps_content_distinct_244_bases() {
		List<BlockSplitSearch.NamedBase> pool = BlockSplitSearch.buildPool(RecombinationPoolConfig.includeDerived(), "Q");
		long distinct244 = pool.stream()
				.map(BlockSplitSearch.NamedBase::base)
				.filter(b -> b.n == 2 && b.m == 4 && b.p == 4 && b.r == 26)
				.map(b -> eu.solven.matmul.catalog.SchemeIO.contentHash(b))
				.distinct().count();
		assertThat(distinct244)
				.as("pool must keep ≥2 content-distinct ⟨2,4,4⟩=26 bases (they recombine "
						+ "differently); shape:r dedup would collapse to 1 and lose the better base")
				.isGreaterThanOrEqualTo(2);
	}

	@Test
	public void extended_pool_is_not_empty() {
		assertThat(BlockSplitSearch.extendedPool(8))
				.as("extendedPool(8) must load catalog leaves from known/derived/curated, not be seed-only")
				.hasSizeGreaterThan(50);
	}

	/**
	 * ⟨5,32,32⟩ = 3320 (= FMM-Lille) must be reachable through the SchemeSweep
	 * evaluate path — {@code buildPool(includeDerived)} + {@code findBestStrategy}
	 * — via the HK ⟨2,4,4⟩=26 recombination: allocA=[3,2] (n: 5=3+2),
	 * allocB=allocC=[8,8,8,8] (each 32=4·8) → 16×⟨3,8,8⟩=145 + 10×⟨2,8,8⟩=100.
	 *
	 * <p>The committed catalog held 3446 because the DEFAULT {@code rootPool}
	 * omitted ⟨2,4,4⟩ as an outer base; until 2026-09-30 the ⟨2,4,4⟩ base lived
	 * only in the derived-inclusive (extended) pool — since issue #8 two reps are
	 * {@code rootPool} entries too (labelled {@code HK<2,4,4>=26 (… rep)}), so
	 * the winner may carry either the extended-pool spelling ({@code 2x4x4}) or
	 * the root-pool one ({@code 2,4,4}). This guards the mechanism: a regression
	 * that drops ⟨2,4,4⟩ from both pools, breaks 4-way ({@code [8,8,8,8]})
	 * allocations, or loses the ⟨3,8,8⟩/⟨2,8,8⟩ leaves would push this back to
	 * 3446 and fail. SOTA-or-better (≤), so a future improvement never breaks it.</p>
	 */
	@Test
	public void includeDerived_sweep_finds_5x32x32_3320_via_2x4x4() {
		List<BlockSplitSearch.NamedBase> pool = BlockSplitSearch.buildPool(RecombinationPoolConfig.includeDerived());
		CitedBound sota = new CitedBound(lookup);
		// bound just above 3320 so the recombination B&B prunes hard and the test
		// stays fast, while still letting the 3320 route through.
		Optional<BlockSplitSearch.NonCubicStrategy> best = BlockSplitSearch.findBestStrategy(
				5, 32, 32, pool, sota, false,
				RecombinationPoolConfig.UNBOUNDED_IMBALANCE, RecombinationPoolConfig.UNBOUNDED_COMBINATIONS, 0, 3446L);
		assertThat(best).as("⟨5,32,32⟩ must resolve via the includeDerived pool").isPresent();
		assertThat(best.get().rank())
				.as("⟨5,32,32⟩ must reach FMM's 3320 or better (regression → 3446 = ⟨2,4,4⟩ base lost)")
				.isLessThanOrEqualTo(3320L);
		assertThat(best.get().recombination())
				.as("the 3320 route is a recombination, not concat/kronecker").isNotNull();
		assertThat(best.get().label())
				.as("the winning outer base must be ⟨2,4,4⟩ (extended-pool '2x4x4' or root-pool '2,4,4' spelling)")
				.containsAnyOf("2x4x4", "2,4,4");
	}
}
