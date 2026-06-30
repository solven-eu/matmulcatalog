package eu.solven.matmul.search;

import java.util.Optional;

import eu.solven.matmul.recombination.Recombination.SotaResolver;
import eu.solven.matmul.papers.dis2009.PanTrilinearAggregation;

/**
 * {@link ConstructiveMethod} wrapper around
 * {@link PanTrilinearAggregation#cubicBound}. Applies to cubic targets
 * {@code ⟨n,n,n⟩} with n ≥ 2.
 *
 * <p>Closed-form Islam 2009 / DIS09 §3 odd-n recipe:</p>
 * <pre>
 *   n odd:  R(⟨n,n,n⟩) ≤ (n³ + 15n² + 14n − 6) / 3
 *   n even: R(⟨n,n,n⟩) ≤ (n³ + 12n² + 11n) / 3
 * </pre>
 *
 * <p>Non-commutative; valid over any field of characteristic 0 (the
 * coefficients are rational, denominators bounded by {@code n/2 + 1}).</p>
 *
 * <p>Useful as a baseline at cubic shapes where no better composition
 * exists. The Sedoglavic Prop 1 bound (the Strassen recombination multiset
 * at a {@code [u,v]} split) usually beats this one at the same shape
 * (e.g. ⟨11⟩³: Sedoglavic 873 vs Pan 1098) — but that bound is reached by
 * the recombination path, not by a registered method.</p>
 */
public final class PanTrilinearAggregationMethod implements ConstructiveMethod {

	@Override
	public String name() { return "PanTrilinearAggregation"; }

	@Override
	public String paperRef() { return "Pan 1980 / Islam 2009 / DIS09 §3 odd-n closed form"; }

	@Override
	public Optional<Prediction> predict(int n, int m, int p, SotaResolver sota) {
		if (n != m || m != p) return Optional.empty(); // cubic only
		if (n < 2) return Optional.empty();
		long r = PanTrilinearAggregation.cubicBound(n);
		String parityTag = (n % 2 == 0) ? "even" : "odd";
		String lineage = String.format("PanTrilinearAggregation(n=%d, %s)=%d", n, parityTag, r);
		return Optional.of(new Prediction(r, "PanTA", lineage, this));
	}
}
