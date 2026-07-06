package eu.solven.matmul.recombination;

import static org.assertj.core.api.Assertions.assertThat;

import java.util.List;
import java.util.Optional;
import java.util.Set;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.recombination.BlockSplitSearch.CandidateKind;
import eu.solven.matmul.recombination.BlockSplitSearch.NonCubicStrategy;
import eu.solven.matmul.search.CitedBound;
import eu.solven.matmul.search.RecursiveMaterialiser;

/**
 * Guards the {@link CandidateKind} election filter of
 * {@link BlockSplitSearch#findBestStrategy}: each {@code --strategies} token maps
 * 1:1 to the kinds it may elect. Historically the {@code recombination} token ran
 * the whole generic upward search and could silently return a Kronecker or concat
 * pick — the same token-dishonesty family as the old {@code kron} alias.
 */
public class TestStrategyKinds {

	private static final FieldAwareLookup LOOKUP = new FieldAwareLookup("Q");
	private static final CitedBound SOTA = new CitedBound(LOOKUP);
	/** One Strassen base — enough for ⟨4,4,4⟩ to admit kron, concat AND recombination. */
	private static final List<BlockSplitSearch.NamedBase> POOL = List.of(
			new BlockSplitSearch.NamedBase("2x2x2",
					LOOKUP.find(2, 2, 2).orElseThrow()));

	private static Optional<NonCubicStrategy> search(Set<CandidateKind> kinds) {
		return BlockSplitSearch.findBestStrategy(4, 4, 4, POOL, SOTA, false,
				Integer.MAX_VALUE, Integer.MAX_VALUE, 0, Long.MAX_VALUE, kinds);
	}

	@Test
	public void kron_kind_elects_only_kronecker() {
		NonCubicStrategy s = search(Set.of(CandidateKind.KRONECKER)).orElseThrow();
		assertThat(s.kronecker()).isNotNull();
		assertThat(BlockSplitSearch.kindOf(s)).isEqualTo(CandidateKind.KRONECKER);
		assertThat(s.rank()).isEqualTo(49);   // Strassen ⊗ Strassen
	}

	@Test
	public void concat_kind_elects_only_concat() {
		NonCubicStrategy s = search(Set.of(CandidateKind.CONCAT)).orElseThrow();
		assertThat(s.concat()).isNotNull();
		assertThat(BlockSplitSearch.kindOf(s)).isEqualTo(CandidateKind.CONCAT);
	}

	@Test
	public void recombination_kind_no_longer_leaks_kron_or_concat() {
		// ⟨4,4,4⟩ over the ⟨2,2,2⟩ base: the balanced [2,2]³ recombination ties the
		// Kron product at 49 — but with only RECOMBINATION electable the returned
		// strategy must BE a recombination, never the (equal-rank) Kronecker pick.
		NonCubicStrategy s = search(Set.of(CandidateKind.RECOMBINATION)).orElseThrow();
		assertThat(s.recombination()).isNotNull();
		assertThat(s.kronecker()).isNull();
		assertThat(s.concat()).isNull();
		assertThat(BlockSplitSearch.kindOf(s)).isEqualTo(CandidateKind.RECOMBINATION);
	}

	@Test
	public void all_kinds_matches_the_legacy_overload() {
		Optional<NonCubicStrategy> legacy = BlockSplitSearch.findBestStrategy(
				4, 4, 4, POOL, SOTA, false, Integer.MAX_VALUE, Integer.MAX_VALUE, 0, Long.MAX_VALUE);
		Optional<NonCubicStrategy> all = search(BlockSplitSearch.ALL_KINDS);
		assertThat(all.map(NonCubicStrategy::rank)).isEqualTo(legacy.map(NonCubicStrategy::rank));
	}

	@Test
	public void token_to_kind_mapping_is_one_to_one() {
		assertThat(RecursiveMaterialiser.kindsFor(Set.of(RecursiveMaterialiser.STRAT_KRONECKER)))
				.containsExactly(CandidateKind.KRONECKER);
		assertThat(RecursiveMaterialiser.kindsFor(Set.of(RecursiveMaterialiser.STRAT_CONCAT)))
				.containsExactly(CandidateKind.CONCAT);
		// Pan-TA pair fusion is NOT a kind: it is a saving within recombination.
		assertThat(RecursiveMaterialiser.kindsFor(Set.of(RecursiveMaterialiser.STRAT_RECOMBINATION)))
				.containsExactlyInAnyOrder(CandidateKind.RECOMBINATION, CandidateKind.METHOD);
		// Serendipitous / projection are not upward-search kinds at all.
		assertThat(RecursiveMaterialiser.kindsFor(Set.of(
				RecursiveMaterialiser.STRAT_SERENDIPITOUS, RecursiveMaterialiser.STRAT_PROJECTION)))
				.isEmpty();
	}
}
