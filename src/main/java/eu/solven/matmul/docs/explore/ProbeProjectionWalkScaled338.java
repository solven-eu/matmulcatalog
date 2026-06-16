package eu.solven.matmul.docs.explore;

import java.util.ArrayList;
import java.util.List;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.ProjectionSearch;
import eu.solven.matmul.search.flip.FlipGraphWalk;
import eu.solven.matmul.search.flip.FlipObjectives;
import eu.solven.matmul.search.flip.FlipScheme;
import lombok.extern.slf4j.Slf4j;

/**
 * SCALED-LIFT margin walk (user "Go", 2026-06-12): unlocks RATIONAL bases for
 * the integer flip engine. Any rational scheme lifts to an integer
 * decomposition of {@code c·T}: clear each product's denominators per factor
 * (u←sU·u, v←sV·v, w←sW·w), then equalize the per-product scale
 * {@code P_l = sU·sV·sW} to a common {@code c = max P_l} by multiplying w by
 * {@code c/P_l}. Flips preserve the sum (so every vertex still decomposes
 * c·T), and the projection objective is SUPPORT-based — scale-invariant — so
 * the walk is meaningful verbatim. Un-scale (divide W by c) at the end.
 *
 * <p>Basin: the rank-55 rational ⟨3,3,8⟩ schemes — the catalog's ACTUAL
 * projection parents for ⟨3,3,7⟩=49 (proj 49, the tie). Success: projected
 * rank ≤ 48 → materialize + verify → catalog win.</p>
 */
@Slf4j
public class ProbeProjectionWalkScaled338 {

	public static void main(String[] args) throws Exception {
		FieldAwareLookup q = new FieldAwareLookup(Field.Q);
		List<NonCubicBilinearAlgorithm> seeds = new ArrayList<>();
		q.findByHash(3, 3, 8, "b766126").ifPresent(ws -> seeds.add(ws.alg()));

		var objective = FlipObjectives.projectedTo(3, 3, 7, 1);
		long bestOverall = Long.MAX_VALUE;
		for (NonCubicBilinearAlgorithm seed : seeds) {
			Scaled lift = lift(seed);
			log.info("lifted r={} c={} maxAbs={} scaledBrent={} seedProj={}",
					seed.r, lift.c, maxAbs(lift.alg), scaledBrentOk(lift),
					ProjectionSearch.projectedRank(lift.alg, 3, 3, 7, 1));
			for (long rng = 1; rng <= 6; rng++) {
				FlipScheme fs = FlipScheme.of(lift.alg);
				FlipGraphWalk.Config cfg = new FlipGraphWalk.Config(400_000, rng, 16, 0.02,
						5_000, 2, 100_000, 0.02, false);
				FlipGraphWalk.Result res = FlipGraphWalk.walk(fs, objective, cfg);
				bestOverall = Math.min(bestOverall, res.bestCost());
				log.info("scaled r{} rng={}: bestProj={} (seed {}) rank={} restarts={}",
						seed.r, rng, res.bestCost(), res.seedCost(), res.best().rank(),
						res.restarts());
				if (res.bestCost() <= 48) {
					NonCubicBilinearAlgorithm unscaled = unscale(res.best().toAlgorithm(), lift.c);
					boolean spot = Verifier.passesRandomMatmulSpotCheck(unscaled);
					var hit = ProjectionSearch.bestFor(3, 3, 7, List.of(unscaled), 49, 1);
					log.info("WIN candidate: walked base r={} unscaled spotCheck={} → projection {}",
							unscaled.r, spot, hit.map(h -> h.rank() + " (spot "
									+ Verifier.passesRandomMatmulSpotCheck(h.scheme()) + ")")
									.orElse("FAILED TO MATERIALIZE"));
					return;
				}
			}
		}
		log.info("no improvement: best projected {} vs catalog 49 (rank-55 rational basin, "
				+ "6 walks/seed)", bestOverall);
	}

	record Scaled(NonCubicBilinearAlgorithm alg, long c) {}

	static Scaled lift(NonCubicBilinearAlgorithm a) {
		double[][] u = a.denseU();
		double[][] v = a.denseV();
		double[][] w = a.denseW();
		int r = a.r;
		long[] pu = new long[r];
		long[] pv = new long[r];
		long[] pw = new long[r];
		long c = 1;
		for (int l = 0; l < r; l++) {
			pu[l] = clearScale(u, l);
			pv[l] = clearScale(v, l);
			pw[l] = clearScale(w, l);
			c = Math.max(c, pu[l] * pv[l] * pw[l]);
		}
		double[][] u2 = new double[u.length][r];
		double[][] v2 = new double[v.length][r];
		double[][] w2 = new double[w.length][r];
		for (int l = 0; l < r; l++) {
			long fix = c / (pu[l] * pv[l] * pw[l]);
			if (fix * pu[l] * pv[l] * pw[l] != c) {
				throw new IllegalStateException("non-dyadic scale mix at product " + l);
			}
			for (int x = 0; x < u.length; x++) {
				u2[x][l] = Math.rint(u[x][l] * pu[l]);
			}
			for (int x = 0; x < v.length; x++) {
				v2[x][l] = Math.rint(v[x][l] * pv[l]);
			}
			for (int x = 0; x < w.length; x++) {
				w2[x][l] = Math.rint(w[x][l] * pw[l] * fix);
			}
		}
		return new Scaled(new NonCubicBilinearAlgorithm(a.n, a.m, a.p, u2, v2, w2), c);
	}

	/** Smallest power of 2 making column l of the factor integer. */
	private static long clearScale(double[][] mat, int l) {
		long s = 1;
		while (s <= 64) {
			boolean ok = true;
			for (double[] row : mat) {
				double x = row[l] * s;
				if (Math.abs(x - Math.rint(x)) > 1e-9) {
					ok = false;
					break;
				}
			}
			if (ok) {
				return s;
			}
			s *= 2;
		}
		throw new IllegalStateException("denominator not dyadic ≤64 in product " + l);
	}

	static NonCubicBilinearAlgorithm unscale(NonCubicBilinearAlgorithm a, long c) {
		double[][] w = a.denseW();
		for (double[] row : w) {
			for (int l = 0; l < row.length; l++) {
				row[l] /= c;
			}
		}
		return new NonCubicBilinearAlgorithm(a.n, a.m, a.p, a.denseU(), a.denseV(), w);
	}

	/** Exact check that the lifted integer factors sum to c·T_matmul. */
	static boolean scaledBrentOk(Scaled s) {
		NonCubicBilinearAlgorithm a = s.alg();
		double[][] u = a.denseU();
		double[][] v = a.denseV();
		double[][] w = a.denseW();
		for (int i = 0; i < a.n; i++) {
			for (int j = 0; j < a.m; j++) {
				for (int j2 = 0; j2 < a.m; j2++) {
					for (int k = 0; k < a.p; k++) {
						for (int i2 = 0; i2 < a.n; i2++) {
							for (int k2 = 0; k2 < a.p; k2++) {
								long sum = 0;
								for (int l = 0; l < a.r; l++) {
									sum += Math.round(u[i * a.m + j][l])
											* Math.round(v[j2 * a.p + k][l])
											* Math.round(w[i2 * a.p + k2][l]);
								}
								long expected = (i == i2 && j == j2 && k == k2) ? s.c() : 0;
								if (sum != expected) {
									return false;
								}
							}
						}
					}
				}
			}
		}
		return true;
	}

	private static long maxAbs(NonCubicBilinearAlgorithm a) {
		long max = 0;
		for (double[][] mat : new double[][][] { a.denseU(), a.denseV(), a.denseW() }) {
			for (double[] row : mat) {
				for (double x : row) {
					max = Math.max(max, Math.round(Math.abs(x)));
				}
			}
		}
		return max;
	}
}
