package eu.solven.matmul.commutative;

import eu.solven.matmul.papers.waksman1970.Waksman1970;

import eu.solven.matmul.NonBilinearAlgorithm;

public class DebugWaksman {
	public static void main(String[] a) {
		int n = 2;
		NonBilinearAlgorithm alg = Waksman1970.build(n);
		int m = alg.m, p = alg.p, r = alg.r;

		double sumSq = 0;
		System.out.println("=== A·B constraint violations ===");
		for (int i = 0; i < n; i++) for (int l = 0; l < p; l++) {
			int il = i*p+l;
			for (int al = 0; al < n; al++) for (int be = 0; be < m; be++)
			for (int ga = 0; ga < m; ga++) for (int de = 0; de < p; de++) {
				int ab = al*m+be, gd = ga*p+de;
				double approx = 0;
				for (int k = 0; k < r; k++) approx += alg.W[il][k] * (alg.Ua[ab][k]*alg.Vb[gd][k] + alg.Va[ab][k]*alg.Ub[gd][k]);
				double target = (al==i && be==ga && de==l) ? 1.0 : 0.0;
				double diff = target - approx;
				if (Math.abs(diff) > 1e-9) System.out.printf("  c[%d,%d] coef of A[%d,%d]·B[%d,%d]: got %.4f want %.4f%n", i,l,al,be,ga,de,approx,target);
				sumSq += diff*diff;
			}
		}
		System.out.println("=== A·A constraint violations ===");
		for (int i = 0; i < n; i++) for (int l = 0; l < p; l++) {
			int il = i*p+l;
			for (int ab = 0; ab < n*m; ab++) for (int ab2 = ab; ab2 < n*m; ab2++) {
				double approx = 0;
				for (int k = 0; k < r; k++) {
					if (ab == ab2) approx += alg.W[il][k]*alg.Ua[ab][k]*alg.Va[ab][k];
					else approx += alg.W[il][k]*(alg.Ua[ab][k]*alg.Va[ab2][k] + alg.Ua[ab2][k]*alg.Va[ab][k]);
				}
				if (Math.abs(approx) > 1e-9) System.out.printf("  c[%d,%d] coef of A[%d]·A[%d]: %.4f%n", i,l,ab,ab2,approx);
				sumSq += approx*approx;
			}
		}
		System.out.println("=== B·B constraint violations ===");
		for (int i = 0; i < n; i++) for (int l = 0; l < p; l++) {
			int il = i*p+l;
			for (int gd = 0; gd < m*p; gd++) for (int gd2 = gd; gd2 < m*p; gd2++) {
				double approx = 0;
				for (int k = 0; k < r; k++) {
					if (gd == gd2) approx += alg.W[il][k]*alg.Ub[gd][k]*alg.Vb[gd][k];
					else approx += alg.W[il][k]*(alg.Ub[gd][k]*alg.Vb[gd2][k] + alg.Ub[gd2][k]*alg.Vb[gd][k]);
				}
				if (Math.abs(approx) > 1e-9) System.out.printf("  c[%d,%d] coef of B[%d]·B[%d]: %.4f%n", i,l,gd,gd2,approx);
				sumSq += approx*approx;
			}
		}
		System.out.printf("Total sumSq = %.6f  residual = %.6f%n", sumSq, Math.sqrt(sumSq));
	}
}
