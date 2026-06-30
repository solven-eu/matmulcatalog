package eu.solven.matmul.search;

import java.util.Optional;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.papers.khoruzhii2026.LitaTaConstruction;
import eu.solven.matmul.papers.khoruzhii2026.LitaTrilinearAggregation;
import eu.solven.matmul.recombination.Recombination.AlgorithmLookup;
import eu.solven.matmul.recombination.Recombination.SotaResolver;

/**
 * {@link ConstructiveMethod} for the Khoruzhii–Gelß–Pokutta 2026 <em>Local
 * Improvements to Trilinear Aggregation</em> (LITA) cubic construction. Applies
 * to cubic {@code ⟨n,n,n⟩} with {@code n ≥ 19} (the construction's domain).
 *
 * <p>Closed-form ranks (Q-rational, non-commutative; see {@link
 * LitaTrilinearAggregation}):</p>
 * <pre>
 *   n even: (4n³ + 45n² + 116n + 84)/12
 *   n odd:  (4n³ + 57n² + 14n − 15)/12 − ⌊3(n−1)/8⌋
 * </pre>
 *
 * <p>BUILDABLE: materialises explicit factors via {@link LitaTaConstruction#build(int)}.
 * The emitted lineage {@code TA_lita(n=N)} is replayable by {@link LineageReplayer}
 * (so a stub referencing it can be expanded on demand). Odd N are exact-verifiable;
 * even N are dense and verify via the spot-check path — but both are realisable
 * schemes at the predicted rank, so predictions are {@code verified}.</p>
 *
 * <p>This is an ALTERNATIVE to {@link PanTrilinearAggregationMethod} (TA_dis); the
 * search adds both as candidates and keeps the min. LITA wins for odd {@code n≥19}
 * and large even {@code n} (see {@code TrilinearAggregations}).</p>
 */
public final class LitaTrilinearAggregationMethod implements ConstructiveMethod {

	@Override
	public String name() {
		return "LitaTrilinearAggregation";
	}

	@Override
	public String paperRef() {
		return "Khoruzhii–Gelß–Pokutta 2026 (LITA); github.com/khoruzhii/lita";
	}

	@Override
	public Optional<Prediction> predict(int n, int m, int p, SotaResolver sota) {
		if (n != m || m != p) {
			return Optional.empty(); // cubic only
		}
		if (n < LitaTrilinearAggregation.MIN_N) {
			return Optional.empty();
		}
		long r = LitaTrilinearAggregation.cubicRank(n);
		// Replayable parametric lineage — LineageReplayer.resolveParametric maps
		// TA_lita(n=N) → LitaTaConstruction.build(N).
		String lineage = "TA_lita(n=" + n + ")";
		return Optional.of(new Prediction(r, "TA_lita", lineage, this));
	}

	@Override
	public Optional<NonCubicBilinearAlgorithm> construct(int n, int m, int p, AlgorithmLookup atoms) {
		if (n != m || m != p || n < LitaTrilinearAggregation.MIN_N) {
			return Optional.empty();
		}
		return Optional.of(LitaTaConstruction.build(n));
	}
}
