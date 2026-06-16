package eu.solven.matmul.search;

import java.util.Optional;

import eu.solven.matmul.catalog.Recombination.SotaResolver;
import eu.solven.matmul.papers.sedoglavic2017.SedoglavicProp1;

/**
 * {@link ConstructiveMethod} wrapper around
 * {@link eu.solven.matmul.papers.sedoglavic2017.SedoglavicProp1#predict}.
 * Applies only to cubic targets {@code ⟨n,n,n⟩} with n ≥ 3.
 *
 * <p>For (u,v) splits with u &gt; v: predicts
 * {@code ⟨u,u,u⟩ + 3·⟨u,u,v⟩ + 3·⟨v,v,u⟩}. For u = v = k: uses the Pan TA
 * pair-cost extension. Returns the minimum over all valid splits.</p>
 *
 * <p>Wired into {@link BlockSplitSearch#findBestStrategy} via
 * {@link MethodCatalog#all()}.</p>
 */
public final class SedoglavicProp1Method implements ConstructiveMethod {

	@Override
	public String name() { return "SedoglavicProp1"; }

	@Override
	public String paperRef() { return "Sedoglavic 2017 hal-01572046v2"; }

	@Override
	public Optional<Prediction> predict(int n, int m, int p, SotaResolver sota) {
		if (n != m || m != p) return Optional.empty(); // cubic only
		if (n < 3) return Optional.empty();
		return SedoglavicProp1.predict(n, sota).map(sp -> {
			String label = sp.usesPanPairCost()
					? "Sedoglavic-doubling[k=" + sp.u() + "]"
					: "Sedoglavic-prop1[u=" + sp.u() + ",v=" + sp.v() + "]";
			return new Prediction(sp.predictedRank(), label, sp.lineageCompact(), this);
		});
	}
}
