package eu.solven.matmul.docs.explore;

import java.io.File;

import eu.solven.matmul.FactorMatrix;
import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.catalog.Compose;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.search.flip.FlipScheme;

/**
 * Throwaway: dump the surviving products of a single-index projection so the
 * second-order merge can be read off by hand. Args: {@code file axis dropIndex}
 * where axis ∈ {n,m,p}. Prints each surviving product as its U, V, W column
 * vectors (over the projected shape), then the merge-aware rank.
 */
public final class ProbeProjectionMerge {
	private ProbeProjectionMerge() {}

	public static void main(String[] args) throws Exception {
		NonCubicBilinearAlgorithm base = SchemeIO.read(new File(args[0]));
		int axis = "nmp".indexOf(args[1].charAt(0));
		int drop = Integer.parseInt(args[2]);
		int[] keepN = keepAllBut(base.n, axis == 0 ? drop : -1);
		int[] keepM = keepAllBut(base.m, axis == 1 ? drop : -1);
		int[] keepP = keepAllBut(base.p, axis == 2 ? drop : -1);
		System.out.printf("== PARENT ⟨%d,%d,%d⟩ r=%d  (W rows = a·%d+c) ==%n",
				base.n, base.m, base.p, base.r, base.p);
		{
			FactorMatrix u = base.u(), v = base.v(), w = base.w();
			for (int l = 0; l < base.r; l++) {
				System.out.printf("  p%-2d  U=%s  V=%s  W=%s%n", l, col(u, l), col(v, l), col(w, l));
			}
		}
		NonCubicBilinearAlgorithm proj = Compose.project(base, keepN, keepM, keepP);
		System.out.printf("%nparent ⟨%d,%d,%d⟩ r=%d  →  drop %s[%d]  →  ⟨%d,%d,%d⟩ dce=%d%n",
				base.n, base.m, base.p, base.r, args[1], drop, proj.n, proj.m, proj.p, proj.r);
		FactorMatrix u = proj.u(), v = proj.v(), w = proj.w();
		System.out.printf("  layout: U rows=(a·%d+b) [%dx%d], V rows=(b·%d+c) [%dx%d], W rows=(a·%d+c) [%dx%d]%n",
				proj.m, proj.n, proj.m, proj.p, proj.m, proj.p, proj.p, proj.n, proj.p);
		for (int l = 0; l < proj.r; l++) {
			System.out.printf("  P%-2d  U=%s  V=%s  W=%s%n", l, col(u, l), col(v, l), col(w, l));
		}
		FlipScheme fs = FlipScheme.of(proj);
		fs.reduce();
		System.out.printf("  merge-aware rank = %d  (second-order margin = %d)%n",
				fs.rank(), proj.r - fs.rank());
	}

	private static String col(FactorMatrix f, int c) {
		StringBuilder sb = new StringBuilder("[");
		for (int r = 0; r < f.rows(); r++) {
			if (r > 0) sb.append(' ');
			double x = f.get(r, c);
			sb.append(x == Math.rint(x) ? Integer.toString((int) x) : Double.toString(x));
		}
		return sb.append(']').toString();
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
