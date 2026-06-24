package eu.solven.matmul.search;

import eu.solven.matmul.recombination.BlockSplitSearch;

import static org.assertj.core.api.Assertions.assertThat;

import java.util.List;
import java.util.Optional;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.recombination.Recombination.SotaResolver;
import eu.solven.matmul.recombination.BlockSplitSearch.NamedBase;
import eu.solven.matmul.recombination.BlockSplitSearch.NonCubicStrategy;

/**
 * The search must not rely on UNVERIFIED (theoretical / formula-only) bounds.
 * Canonical trigger: HK1971 predicts ⟨2,10,15⟩=233, one below the buildable
 * Perminov 234, but we have no construction for it. The predictor must still
 * EMIT the bound (flagged unverified, for display), yet {@code findBestStrategy}
 * must never ELECT it.
 */
public class TestUnverifiedBoundsIgnored {

	private static SotaResolver sota(FieldAwareLookup lk) {
		return (p, q, r) -> {
			if (p == 0 || q == 0 || r == 0) return 0;
			if (p == 1) return q * r;
			if (q == 1) return p * r;
			if (r == 1) return p * q;
			int v = lk.findRank(p, q, r);
			return v >= Integer.MAX_VALUE / 100 ? p * q * r : v;
		};
	}

	@Test
	public void hk1971PredictionIsEmittedButFlaggedUnverified() {
		FieldAwareLookup lk = new FieldAwareLookup("R");
		List<ConstructiveMethod.Prediction> preds = MethodCatalog.predictAll(2, 10, 15, sota(lk));
		Optional<ConstructiveMethod.Prediction> hk =
				preds.stream().filter(p -> "HK1971".equals(p.label())).findFirst();
		assertThat(hk).as("HK1971 still predicts (for display)").isPresent();
		assertThat(hk.get().predictedRank()).isEqualTo(233);
		assertThat(hk.get().verified()).as("but the HK ⟨2,m,n⟩ formula bound is unverified").isFalse();
	}

	@Test
	public void searchDoesNotElectTheUnverifiedHkBound() {
		FieldAwareLookup lk = new FieldAwareLookup("R");
		List<NamedBase> pool = BlockSplitSearch.buildPool(PoolConfig.simple());
		Optional<NonCubicStrategy> best =
				BlockSplitSearch.findBestStrategy(2, 10, 15, pool, sota(lk), false);
		assertThat(best).isPresent();
		// The elected strategy must be a real one, not the unverified HK formula.
		assertThat(best.get().label()).as("must not elect HK1971").isNotEqualTo("HK1971");
		assertThat(best.get().rank()).as("must not be the unbuildable 233").isGreaterThanOrEqualTo(234L);
	}
}
