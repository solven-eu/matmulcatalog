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
	 * omits ⟨2,4,4⟩ as an outer base; the ⟨2,4,4⟩ base lives only in the
	 * derived-inclusive (extended) pool. This guards the mechanism: a regression
	 * that drops ⟨2,4,4⟩ from the extended pool, breaks 4-way ({@code [8,8,8,8]})
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
				.as("the winning outer base must be ⟨2,4,4⟩").contains("2x4x4");
	}
}
