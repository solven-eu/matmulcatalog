package eu.solven.matmul.papers.hopcroftkerr1971;

import java.util.Optional;

import eu.solven.matmul.recombination.Recombination.SotaResolver;
import eu.solven.matmul.search.ConstructiveMethod;
import eu.solven.matmul.search.ConstructiveMethod.Prediction;

/**
 * {@link ConstructiveMethod} wrapper around
 * {@link HopcroftKerrBound#forShape}. Applies to the {@code ⟨a,b,c⟩}
 * family where at least one axis equals 2 — the canonical Hopcroft-Kerr
 * 1971 setup.
 *
 * <p>Closed-form predictor: {@code R(⟨x,2,y⟩) ≤ ⌈(3xy + max(x,y))/2⌉}.
 * Returns the constructive recipe in lineage_compact form for the
 * winning axis orientation.</p>
 *
 * <p>Non-commutative; valid over any algebra.</p>
 */
public final class HopcroftKerr1971Method implements ConstructiveMethod {

	@Override
	public String name() { return "HopcroftKerr1971"; }

	@Override
	public String paperRef() { return "Hopcroft-Kerr 1971 (Cornell TR 69-44 / SIAM 1971)"; }

	@Override
	public Optional<Prediction> predict(int n, int m, int p, SotaResolver sota) {
		long r = HopcroftKerrBound.forShape(n, m, p);
		if (r < 0) return Optional.empty();
		String lineage = String.format("HopcroftKerr1971(%d,%d,%d)=%d", n, m, p, r);
		// UNVERIFIED: this is the HK closed-form bound only. We have not ported a
		// construction that covers every parity/range (e.g. even-p asymmetric
		// ⟨2,10,15⟩=233 is predicted but not constructible here — the buildable
		// catalog value is Perminov 234). The search must not rely on it; where a
		// real HK scheme exists it is on disk and reached via findRank instead.
		return Optional.of(new Prediction(r, "HK1971", lineage, this, false));
	}
}
