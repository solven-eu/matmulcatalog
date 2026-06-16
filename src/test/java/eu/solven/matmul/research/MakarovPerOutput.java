package eu.solven.matmul.research;

import eu.solven.matmul.NonBilinearAlgorithm;
import eu.solven.matmul.papers.makarov1986.Makarov22;

/**
 * Per-output residual diagnostic for Makarov22. Prints each c_{i,j}'s
 * sum-of-squared-error so we can see which output formulas are wrong.
 */
public final class MakarovPerOutput {
	private MakarovPerOutput() {}

	public static void main(String[] args) {
		NonBilinearAlgorithm alg = Makarov22.buildDefault();
		int n = alg.n, m = alg.m, p = alg.p, r = alg.r;

		System.out.println("Per-output sumSq (target − computed coefficient of each monomial):");
		double totalSq = 0;
		for (int i = 0; i < n; i++) {
			for (int l = 0; l < p; l++) {
				int ilIdx = i * p + l;
				double outSq = 0;
				int wrongCount = 0;
				StringBuilder wrongList = new StringBuilder();
				// All A·B cross-terms (commutative-symmetrised).
				for (int alpha = 0; alpha < n; alpha++) {
					for (int beta = 0; beta < m; beta++) {
						int abIdx = alpha * m + beta;
						for (int gamma = 0; gamma < m; gamma++) {
							for (int delta = 0; delta < p; delta++) {
								int gdIdx = gamma * p + delta;
								double approx = 0;
								for (int k = 0; k < r; k++) {
									approx += alg.W[ilIdx][k] *
											(alg.Ua[abIdx][k] * alg.Vb[gdIdx][k]
													+ alg.Va[abIdx][k] * alg.Ub[gdIdx][k]);
								}
								double target = (alpha == i && beta == gamma && delta == l) ? 1.0 : 0.0;
								double diff = target - approx;
								if (Math.abs(diff) > 1e-9) {
									wrongCount++;
									if (wrongList.length() < 200) {
										wrongList.append(String.format("a%d%d·b%d%d:%+g ",
												alpha + 1, beta + 1, gamma + 1, delta + 1, -diff));
									}
								}
								outSq += diff * diff;
							}
						}
					}
				}
				// A·A cross-terms (must be 0).
				for (int alpha = 0; alpha < n; alpha++) {
					for (int beta = 0; beta < m; beta++) {
						int abIdx = alpha * m + beta;
						for (int alphaP = 0; alphaP < n; alphaP++) {
							for (int betaP = 0; betaP < m; betaP++) {
								if ((alpha * m + beta) > (alphaP * m + betaP)) continue;
								int abIdx2 = alphaP * m + betaP;
								double approx = 0;
								for (int k = 0; k < r; k++) {
									approx += alg.W[ilIdx][k] *
											(alg.Ua[abIdx][k] * alg.Va[abIdx2][k]
													+ alg.Ua[abIdx2][k] * alg.Va[abIdx][k]);
								}
								if (Math.abs(approx) > 1e-9) {
									wrongCount++;
									if (wrongList.length() < 200) {
										wrongList.append(String.format("a%d%d·a%d%d:%+g ",
												alpha + 1, beta + 1, alphaP + 1, betaP + 1, approx));
									}
								}
								outSq += approx * approx;
							}
						}
					}
				}
				// B·B cross-terms (must be 0).
				for (int gamma = 0; gamma < m; gamma++) {
					for (int delta = 0; delta < p; delta++) {
						int gdIdx = gamma * p + delta;
						for (int gammaP = 0; gammaP < m; gammaP++) {
							for (int deltaP = 0; deltaP < p; deltaP++) {
								if ((gamma * p + delta) > (gammaP * p + deltaP)) continue;
								int gdIdx2 = gammaP * p + deltaP;
								double approx = 0;
								for (int k = 0; k < r; k++) {
									approx += alg.W[ilIdx][k] *
											(alg.Ub[gdIdx][k] * alg.Vb[gdIdx2][k]
													+ alg.Ub[gdIdx2][k] * alg.Vb[gdIdx][k]);
								}
								if (Math.abs(approx) > 1e-9) {
									wrongCount++;
									if (wrongList.length() < 200) {
										wrongList.append(String.format("b%d%d·b%d%d:%+g ",
												gamma + 1, delta + 1, gammaP + 1, deltaP + 1, approx));
									}
								}
								outSq += approx * approx;
							}
						}
					}
				}
				System.out.printf("c%d%d  sumSq=%6.2f  wrongTerms=%d   %s%n",
						i + 1, l + 1, outSq, wrongCount, wrongList);
				totalSq += outSq;
			}
		}
		System.out.printf("Total sumSq = %.2f  (residual = sqrt = %.4f)%n",
				totalSq, Math.sqrt(totalSq));
	}
}
