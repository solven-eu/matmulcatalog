package eu.solven.matmul.search;

import eu.solven.matmul.recombination.BlockSplitSearch;

import static org.assertj.core.api.Assertions.assertThat;

import java.util.List;
import java.util.Optional;
import java.util.Set;

import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.Lineage;
import eu.solven.matmul.recombination.Recombination;

/**
 * Regression guard for rectangular Pan trilinear-aggregation (TA) <em>within</em>
 * decomposition. TA is not a top-level strategy: the two off-diagonal cross-blocks
 * of a symmetric ⟨1,2,2⟩ peel of ⟨N,N+s,N+s⟩ are cyclic rotations and fuse at
 * {@code fusedRank=N·s·N+N·N+N·s+s·N} (= {@link
 * eu.solven.matmul.recombination.PairedSubProducts#pairCost}) instead of the two leaves'
 * summed rank. The canonical case is ⟨26,29,29⟩ = ⟨26,26,26⟩ 8658 + TA(⟨26,3,26⟩,
 * ⟨26,26,3⟩) 2860 + ⟨26,3,3⟩ 175 = <b>11693</b> (FMM parity; our prior best 11808).
 *
 * <p>Two silent failures this guards:</p>
 * <ol>
 *   <li><b>Scoring</b>: the default unbalanced/no-peel recombination path used the
 *       pairing-BLIND allocation optimizer, so the TA saving was invisible
 *       (findBestStrategy returned 11808/concat, never 11693). The pairing-aware
 *       mask sweep is now run alongside it.</li>
 *   <li><b>Materialisation</b>: a recombination scored with pairing was BUILT by
 *       gluing the cross-pair independently (→ un-fused 11841, rejected). The
 *       ⟨1,2,2⟩ peel now routes through {@code buildPeeledViaTa}; the corner must
 *       resolve to the NON-COMMUTATIVE ⟨26,3,3⟩=175, not the commutative 159.</li>
 * </ol>
 */
public class TestTaPeelDecomposition {

	private static FieldAwareLookup lookup;
	private static List<BlockSplitSearch.NamedBase> peelPool;
	private static Recombination.SotaResolver sota;

	@BeforeAll
	static void setUp() {
		lookup = new FieldAwareLookup("Q");
		// The symmetric-peel carrier: the naïve ⟨1,2,2⟩ scheme (rank 4).
		peelPool = List.of(new BlockSplitSearch.NamedBase(
				"base<1x2x2>=4", NonCubicBilinearAlgorithm.naive(1, 2, 2)));
		sota = (a, b, c) -> lookup.findRank(a, b, c);
	}

	@Test
	public void scoring_considers_TA_for_26x29x29_via_1x2x2() {
		// Scoring guard: the recombination scorer must SEE the TA cross-fusion (≤11693),
		// not just the un-paired 11841 / concat 11839 the optimizer alone reports.
		Optional<BlockSplitSearch.NonCubicStrategy> picked = BlockSplitSearch.findBestStrategy(
				26, 29, 29, peelPool, sota, false,
				RecombinationPoolConfig.UNBOUNDED_IMBALANCE, Integer.MAX_VALUE, 0, Long.MAX_VALUE);
		assertThat(picked).as("a ⟨1,2,2⟩-peel recombination strategy must be found").isPresent();
		assertThat(picked.get().rank())
				.as("⟨26,29,29⟩ via ⟨1,2,2⟩ must be priced ≤ 11693 (TA cross-fusion), got %s",
						picked.get())
				.isLessThanOrEqualTo(11693L);
	}

	@Test
	public void materialise_builds_verified_TA_peel_26x29x29() {
		RecursiveClosureSota closure = new RecursiveClosureSota(lookup, peelPool, true, true);
		// dry-run (writeRoot=null), balancedOnly=false so the unbalanced [26,3] peel is
		// explored, improveExisting=false, deriveBest=true.
		RecursiveMaterialiser mat = new RecursiveMaterialiser(
				lookup, peelPool, closure, null, false, false, false, true);
		mat.setStrategies(Set.of(RecursiveMaterialiser.STRAT_RECOMBINATION));

		Optional<RecursiveMaterialiser.Result> r = mat.materialise(26, 29, 29);
		assertThat(r).as("⟨26,29,29⟩ must materialise via the ⟨1,2,2⟩ TA peel").isPresent();
		assertThat(r.get().alg().r)
				.as("⟨26,29,29⟩ TA peel rank must be ≤ 11693 (regression if higher)")
				.isLessThanOrEqualTo(11693);
		assertThat(r.get().lineage())
				.as("the construction must be recorded as a generic TA-fused recombination node")
				.isInstanceOf(Lineage.RecombinationTaN.class);
		assertThat(Verifier.passesRandomMatmulSpotCheck(r.get().alg()))
				.as("the built ⟨26,29,29⟩ scheme must verify").isTrue();
	}
}
