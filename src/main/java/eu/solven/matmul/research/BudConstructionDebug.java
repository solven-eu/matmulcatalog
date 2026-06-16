package eu.solven.matmul.research;

import java.util.EnumSet;
import java.util.List;
import java.util.Set;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.catalog.SerendipitousBudProduct;
import eu.solven.matmul.catalog.SerendipitousBudProduct.BudType;
import lombok.extern.slf4j.Slf4j;

/** Bisect which bud type makes productViaBuds emit an invalid scheme. */
@Slf4j
public final class BudConstructionDebug {
	private BudConstructionDebug() {}

	public static void main(String[] args) throws Exception {
		FieldAwareLookup lookup = new FieldAwareLookup(Field.Q);
		// Bud-richest ⟨3,4,4⟩ (the failing base ⟨6,16,20⟩ = ⟨3,4,4⟩ ⊗ ⟨2,4,5⟩).
		NonCubicBilinearAlgorithm base = richest(lookup, 3, 4, 4);
		int[] inner = { 2, 4, 5 };
		var dec = SerendipitousBudProduct.findBuds(base);
		long u = dec.buds().stream().filter(b -> b.type() == BudType.U).count();
		long v = dec.buds().stream().filter(b -> b.type() == BudType.V).count();
		long w = dec.buds().stream().filter(b -> b.type() == BudType.W).count();
		log.info("base ⟨3,4,4⟩ r={} buds: U={} V={} W={} trivial={}", base.r, u, v, w, dec.trivial().length);

		for (Set<BudType> allow : List.of(
				EnumSet.noneOf(BudType.class),
				EnumSet.of(BudType.U), EnumSet.of(BudType.V), EnumSet.of(BudType.W),
				EnumSet.allOf(BudType.class))) {
			try {
				NonCubicBilinearAlgorithm built = SerendipitousBudProduct.productViaBudsTyped(
						base, lookup, inner[0], inner[1], inner[2], allow);
				boolean ok = Verifier.passesRandomMatmulSpotCheck(built);
				log.info("  allow={} → r={} verify={}", allow, built.r, ok);
			} catch (Exception e) {
				log.info("  allow={} → ERROR {}", allow, e.toString());
			}
		}
	}

	private static NonCubicBilinearAlgorithm richest(FieldAwareLookup lookup, int n, int m, int p)
			throws Exception {
		NonCubicBilinearAlgorithm best = null;
		int bestScore = -1;
		for (var path : lookup.findFiles(n, m, p)) {
			try {
				NonCubicBilinearAlgorithm a = SchemeIO.read(path.toFile());
				int sc = BudBaseFactory.budScore(a);
				if (sc > bestScore) {
					bestScore = sc;
					best = a;
				}
			} catch (Exception ignored) {
				// skip
			}
		}
		return best;
	}
}
