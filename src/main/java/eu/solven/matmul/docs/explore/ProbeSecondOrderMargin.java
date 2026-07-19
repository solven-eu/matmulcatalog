package eu.solven.matmul.docs.explore;

import java.io.File;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.catalog.Compose;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.search.flip.FlipScheme;

/**
 * Throwaway probe (research question, NOT a committed metric): does projecting a
 * base one index down ever create <em>proportional</em> survivors — i.e. is the
 * "second-order margin" the projection metric ignores actually nonzero in the
 * catalog?
 *
 * <p>{@code Compose.project} returns the DCE'd projection ({@code .r} =
 * {@code ProjectionSearch.survivorCount} — first-order margin only). Wrapping it
 * in {@code FlipScheme.of(..).reduce()} additionally runs {@code mergeOnce} to
 * fixpoint (drop-zero + two-slot ± merges). The gap {@code dce − reduced} is the
 * second-order margin the current metric never counts.</p>
 *
 * <p>Args: one or more scheme JSON paths. For each, tries every single-index drop
 * on each axis and reports the largest second-order margin seen.</p>
 */
public final class ProbeSecondOrderMargin {
	private ProbeSecondOrderMargin() {}

	public static void main(String[] args) throws Exception {
		for (String path : args) {
			NonCubicBilinearAlgorithm base;
			try {
				base = SchemeIO.read(new File(path));
			} catch (Exception stubOrUnreadable) {
				continue; // lineage stub / no explicit matrices — skip silently
			}
			System.out.printf("%n=== %s  ⟨%d,%d,%d⟩ r=%d ===%n",
					new File(path).getName(), base.n, base.m, base.p, base.r);
			int bestGap = 0;
			String bestWhere = "none";
			int childrenProbed = 0;
			// axis 0 = n, 1 = m, 2 = p
			for (int axis = 0; axis < 3; axis++) {
				int dim = axis == 0 ? base.n : axis == 1 ? base.m : base.p;
				if (dim <= 1) continue;
				for (int drop = 0; drop < dim; drop++) {
					int[] keepN = keepAllBut(base.n, axis == 0 ? drop : -1);
					int[] keepM = keepAllBut(base.m, axis == 1 ? drop : -1);
					int[] keepP = keepAllBut(base.p, axis == 2 ? drop : -1);
					NonCubicBilinearAlgorithm proj = Compose.project(base, keepN, keepM, keepP);
					int dce = proj.r;                       // first-order (survivorCount)
					FlipScheme fs = FlipScheme.of(proj);
					fs.reduce();                            // + second-order merges
					int reduced = fs.rank();
					int gap = dce - reduced;
					childrenProbed++;
					if (gap > bestGap) {
						bestGap = gap;
						bestWhere = String.format("axis=%s drop=%d → ⟨%d,%d,%d⟩ dce=%d reduced=%d",
								"nmp".charAt(axis), drop, proj.n, proj.m, proj.p, dce, reduced);
					}
				}
			}
			System.out.printf("  probed %d one-drop children; MAX second-order margin = %d  (%s)%n",
					childrenProbed, bestGap, bestWhere);
		}
	}

	private static int[] keepAllBut(int dim, int drop) {
		int keep = drop < 0 ? dim : dim - 1;
		int[] out = new int[keep];
		int w = 0;
		for (int i = 0; i < dim; i++) {
			if (i == drop) continue;
			out[w++] = i;
		}
		return out;
	}
}
