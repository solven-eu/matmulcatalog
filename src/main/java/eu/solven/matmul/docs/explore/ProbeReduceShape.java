package eu.solven.matmul.docs.explore;

import java.util.List;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.recombination.BlockSplitSearch;
import eu.solven.matmul.search.RecursiveClosureSota;
import eu.solven.matmul.search.RecursiveMaterialiser;
import eu.solven.matmul.search.flip.FlipScheme;

/**
 * Throwaway: resolve each requested shape to its catalog-best explicit scheme
 * (materialising the lineage) and check whether the merge-aware {@code reduce()}
 * pass — drop-zeros + two-slot ± fusion to fixpoint — lands below the stored
 * rank. A drop below the catalog rank means the stored SOTA is a loose
 * survivor-count upper bound, not the true rank.
 *
 * <p>Args: shapes as {@code NxMxP}, e.g. {@code 9x17x17 9x17x18 9x17x19}.</p>
 */
public final class ProbeReduceShape {
	private ProbeReduceShape() {}

	public static void main(String[] args) {
		FieldAwareLookup lookup = new FieldAwareLookup("Q");
		List<BlockSplitSearch.NamedBase> pool = BlockSplitSearch.defaultPool();
		RecursiveClosureSota sota = new RecursiveClosureSota(lookup, pool, true, true);
		RecursiveMaterialiser mat = new RecursiveMaterialiser(lookup, pool, sota, null, false, true);
		for (String shape : args) {
			String[] parts = shape.split("x");
			int n = Integer.parseInt(parts[0]), m = Integer.parseInt(parts[1]), p = Integer.parseInt(parts[2]);
			NonCubicBilinearAlgorithm alg = mat.materialise(n, m, p).map(RecursiveMaterialiser.Result::alg).orElse(null);
			if (alg == null) {
				System.out.printf("⟨%d,%d,%d⟩ : not resolvable%n", n, m, p);
				continue;
			}
			int stored = alg.r;
			FlipScheme fs = FlipScheme.of(alg);
			fs.reduce();
			int reduced = fs.rank();
			System.out.printf("⟨%d,%d,%d⟩ stored=%d  reduced=%d  delta=%d%s%n",
					n, m, p, stored, reduced, stored - reduced,
					reduced < stored ? "   *** REDUCE BEATS STORED ***" : "");
		}
	}
}
