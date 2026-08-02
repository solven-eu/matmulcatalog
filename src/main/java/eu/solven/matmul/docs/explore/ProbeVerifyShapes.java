package eu.solven.matmul.docs.explore;

import java.util.List;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.recombination.BlockSplitSearch;
import eu.solven.matmul.search.RecursiveClosureSota;
import eu.solven.matmul.search.RecursiveMaterialiser;
import eu.solven.matmul.verifiers.Verifier;

/**
 * Throwaway: materialise each requested shape from the catalog (replaying the
 * persisted lineage) and run an exact / sampled matmul verification. Args:
 * shapes as {@code NxMxP}. Prints rank and PASS/FAIL per shape.
 */
public final class ProbeVerifyShapes {
	private ProbeVerifyShapes() {}

	public static void main(String[] args) {
		FieldAwareLookup lookup = new FieldAwareLookup("Q");
		List<BlockSplitSearch.NamedBase> pool = BlockSplitSearch.defaultPool();
		RecursiveClosureSota sota = new RecursiveClosureSota(lookup, pool, true, true);
		RecursiveMaterialiser mat = new RecursiveMaterialiser(lookup, pool, sota, null, false, true);
		for (String shape : args) {
			String[] pp = shape.split("x");
			int n = Integer.parseInt(pp[0]), m = Integer.parseInt(pp[1]), p = Integer.parseInt(pp[2]);
			NonCubicBilinearAlgorithm alg =
					mat.materialise(n, m, p).map(RecursiveMaterialiser.Result::alg).orElse(null);
			if (alg == null) {
				System.out.printf("⟨%d,%d,%d⟩ : NOT RESOLVABLE%n", n, m, p);
				continue;
			}
			boolean ok = Verifier.passesRandomMatmulSpotCheck(alg);
			System.out.printf("⟨%d,%d,%d⟩ rank=%d  %s%n", n, m, p, alg.r, ok ? "PASS" : "*** FAIL ***");
		}
	}
}
