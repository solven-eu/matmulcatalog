package eu.solven.matmul.search;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;
import java.util.List;
import java.util.Optional;

import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.catalog.SerendipitousBudProduct;
import eu.solven.matmul.catalog.SerendipitousSearch;

/**
 * Regression guards for the stub-blind serendipitous search (fmm-gap 2026-07-08,
 * ⟨20,28,28⟩). The silent failure: {@code SerendipitousSearch} priced candidates
 * with a {@code findWithSource}-based oracle (−1 for any shape whose best entry
 * is a lineage-only stub) and built through {@code findWithSource} too — so a
 * candidate whose enlarged fusion target existed only as a stub (⟨4,4,20⟩=230,
 * a ConcatCols stub) was silently dropped, and FMM's published
 * {@code (⟨5,7,7⟩:176 − 5) ⊗ ⟨4,4,4⟩:48 + ⟨4,4,20⟩:230 = 8438} construction
 * looked unreachable. Fix: predict via {@code findRank} (stub-inclusive) and
 * build via a stub-capable {@link SerendipitousBudProduct.InnerResolver}.
 */
public class TestSerendipitousStubInner {

	private static final String FMM_BASE =
			"src/main/resources/schemes/bud-bases/section7/fmm-lille_5x7x7_r176_a3315.json";

	private static FieldAwareLookup lookup;
	private static NonCubicBilinearAlgorithm base;

	@BeforeAll
	public static void setUp() throws Exception {
		lookup = new FieldAwareLookup("Q");
		base = SchemeIO.read(new File(FMM_BASE));
	}

	/** Stub-capable resolver — same construction as RecursiveMaterialiser's. */
	private static SerendipitousBudProduct.InnerResolver replayingResolver() {
		LineageReplayer replayer = LineageReplayer.withDefaultPool(lookup);
		return (n, m, p) -> {
			Optional<NonCubicBilinearAlgorithm> direct =
					lookup.findWithSource(n, m, p).map(FieldAwareLookup.WithSource::alg);
			if (direct.isPresent() && lookup.findRank(n, m, p) >= direct.get().r) {
				return direct;
			}
			return lookup.findFile(n, m, p)
					.map(f -> {
						try {
							return replayer.replayFromFile(f.toFile()).orientAs(n, m, p).orElse(null);
						} catch (RuntimeException e) {
							return null;
						}
					})
					.or(() -> direct);
		};
	}

	/**
	 * The FMM Kauers–Wood ⟨5,7,7⟩:176 base carries a size-5 U-bud (⟨1,1,5⟩ class)
	 * plus a size-4 V-bud, pricing ⟨20,28,28⟩ = 8448 − 10 − 4 = 8434 &lt; FMM's
	 * 8438. With a stub-capable resolver the search must BUILD it, not just
	 * price it — the enlarged ⟨4,4,20⟩ inner is a stub.
	 */
	@Test
	public void bestFor_builds_through_stub_inner() {
		Optional<SerendipitousSearch.Hit> hit = SerendipitousSearch.bestFor(
				20, 28, 28, List.of(base), lookup, 8440, replayingResolver());
		assertThat(hit).as("serendipitous ⟨20,28,28⟩ from the fmm ⟨5,7,7⟩ base").isPresent();
		assertThat(hit.get().scheme().r).isLessThanOrEqualTo(8434);
	}

	/**
	 * Phase-1 pricing must be stub-inclusive: with the historical
	 * findWithSource oracle the ⟨4,4,20⟩-fusing candidate predicted −1 and the
	 * whole search came back EMPTY even when handed the winning base directly.
	 * The default-resolver overload must at least not crash on the stub-priced
	 * candidate (it falls through when the build fails) — the WIN path is the
	 * resolver overload above.
	 */
	@Test
	public void default_resolver_does_not_crash_on_stub_priced_candidate() {
		// Must terminate without throwing; result may be empty or a worse build.
		SerendipitousSearch.bestFor(20, 28, 28, List.of(base), lookup, 8440);
	}

	/**
	 * Orientation-ambiguity guard (fmm-gap 2026-07-09, ⟨21,28,30⟩): a ⟨5,7,7⟩
	 * file reaches ⟨7,7,5⟩ by TWO distinct S₃ orientations with different bud
	 * profiles — budBasesAt must offer BOTH (the σ-paying one is the second;
	 * single-orientation offering priced 9477 instead of FMM's 9473), and the
	 * persisted pin must carry the EXACT axisMap so replay is deterministic.
	 */
	@Test
	public void ambiguous_orientation_offers_both_bud_profiles() {
		NonCubicBilinearAlgorithm fmm = base; // native ⟨5,7,7⟩
		long best = Long.MAX_VALUE;
		java.util.Set<Long> costs = new java.util.TreeSet<>();
		for (int[] perm : NonCubicBilinearAlgorithm.ORIENT_PERMS) {
			NonCubicBilinearAlgorithm o = fmm.orientByPerm(perm);
			if (o.n != 7 || o.m != 7 || o.p != 5) continue;
			long c = SerendipitousBudProduct.serendipitousCost(o, lookup, 3, 4, 6);
			costs.add(c);
			best = Math.min(best, c);
		}
		assertThat(costs).as("two orientations reach ⟨7,7,5⟩ with DIFFERENT σ").hasSizeGreaterThan(1);
		assertThat(best).as("the σ-best orientation prices FMM's ⟨21,28,30⟩ recipe").isLessThanOrEqualTo(9473);
	}

	/** The persisted ⟨21,28,30⟩=9473 stub carries an explicit axisMap pin and
	 *  replays to exactly the claimed rank (dim-repeat-safe orientation). */
	@Test
	public void ambiguous_orientation_stub_replays_exactly() {
		File stub = new File(
				"src/main/resources/schemes/derived/section30/21x28x30-r9473-derived-f794828.json");
		assertThat(stub).exists();
		NonCubicBilinearAlgorithm replayed =
				LineageReplayer.withDefaultPool(lookup).replayFromFile(stub);
		assertThat(replayed.r).isEqualTo(9473);
	}

	/**
	 * The persisted ⟨20,28,28⟩=8434 stub replays through its stub inner
	 * (LineageReplayer.resolveInner) to exactly the claimed rank.
	 */
	@Test
	public void serendip_stub_with_stub_inner_replays() {
		File stub = new File(
				"src/main/resources/schemes/derived/section28/20x28x28-r8434-derived-5d8d9ca.json");
		assertThat(stub).exists();
		NonCubicBilinearAlgorithm replayed =
				LineageReplayer.withDefaultPool(lookup).replayFromFile(stub);
		assertThat(replayed.r).isEqualTo(8434);
		assertThat(new int[] { replayed.n, replayed.m, replayed.p }).containsExactly(20, 28, 28);
	}
}
