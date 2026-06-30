package eu.solven.matmul.papers.waksman1970;

import java.util.Optional;

import eu.solven.matmul.recombination.Recombination.SotaResolver;
import eu.solven.matmul.search.ConstructiveMethod;
import eu.solven.matmul.search.MethodCatalog;
import eu.solven.matmul.search.ConstructiveMethod.Prediction;

/**
 * {@link ConstructiveMethod} wrapper around
 * {@link WaksmanBound#forShape}. Applies to any shape {@code ⟨a,b,c⟩}
 * with a, b, c ≥ 1.
 *
 * <p><strong>COMMUTATIVE ONLY.</strong> Waksman 1970's identity exploits
 * scalar commutativity; the resulting scheme does NOT lift to recursive
 * matmul over a non-commutative ring (matrix scalars). Therefore this
 * method is <strong>not registered in {@link MethodCatalog#all()}</strong>
 * (which targets the NC search) — it lives here for use by a future
 * commutative-aware lookup ({@code #65}) and to keep the constructor
 * vocabulary uniform across papers.</p>
 *
 * <p>Closed-form bound:</p>
 * <pre>
 *   b even: R(⟨a,b,c⟩) ≤ b·(a·c + a + c − 1) / 2
 *   b odd:  R(⟨a,b,c⟩) ≤ (b−1)·(a·c + a + c − 1)/2 + a·c
 * </pre>
 */
public final class Waksman1970Method implements ConstructiveMethod {

	@Override
	public String name() { return "Waksman1970"; }

	@Override
	public String paperRef() { return "Waksman 1970"; }

	@Override
	public Optional<Prediction> predict(int n, int m, int p, SotaResolver sota) {
		if (n < 1 || m < 1 || p < 1) return Optional.empty();
		long r = WaksmanBound.forShape(n, m, p);
		String lineage = String.format("Waksman1970(%d,%d,%d)=%d [commutative]", n, m, p, r);
		return Optional.of(new Prediction(r, "Waksman1970-cmt", lineage, this));
	}

	/** True for every Waksman prediction — used by commutative-aware callers. */
	public boolean isCommutative() { return true; }
}
