package eu.solven.matmul.docs.explore;

import java.util.ArrayList;
import java.util.List;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.verifiers.Verifier;
import eu.solven.matmul.papers.dis2009.PanTrilinearAggregation;
import lombok.extern.slf4j.Slf4j;

/**
 * Prototype: does a STRUCTURED LINEAR projection of the PanTA ⟨n,n,n⟩ cube reach
 * FMM's ⟨n−1,n,n⟩ rank, where our exhaustive COORDINATE-drop projection cannot?
 *
 * <p>FMM's {@code 27x28x28=10413} is {@code projection [[1,0],[0]]} of the
 * {@code 28x28x28=10556} PanTA cube. Our {@link eu.solven.matmul.catalog.ProjectionSearch}
 * only quotients a raw n-axis coordinate e_d (best margin μ=114 → 10442). The
 * hypothesis: FMM quotients an AGGREGATION direction {@code v = e_a ± e_d} of the
 * Pan construction, killing every product whose A-form OR C-form is parallel to v
 * — i.e. every aggregated {@code (A_a ± A_d)} / {@code (C_a ± C_d)} product — which
 * a coordinate drop leaves alive. This probe searches that {0,±1} two-coordinate
 * direction family.</p>
 *
 * <p>Quotienting direction v on the n-axis is the symmetric projection
 * (E embeds v⊥, Q projects onto v⊥ along v, Q·E=I) — a proven matmul-preserving
 * operator, so the surviving-product count IS the rank. The winner is rebuilt and
 * spot-checked as ⟨n−1,n,n⟩ to prove correctness by construction.</p>
 *
 * <pre>mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.explore.ProbeStructuredProjection -Dexec.args="28"</pre>
 */
@Slf4j
public final class ProbeStructuredProjection {
	private ProbeStructuredProjection() {}

	private static final double EPS = 1e-9;

	public static void main(String[] args) {
		int N = args.length > 0 ? Integer.parseInt(args[0]) : 28;
		NonCubicBilinearAlgorithm cube = PanTrilinearAggregation.build(N);
		int n = cube.n, m = cube.m, p = cube.p, r = cube.r;
		log.info("PanTA ⟨{},{},{}⟩ r={} exact={}", n, m, p, r,
				Verifier.passesRandomMatmulSpotCheck(cube));
		if (n > 31) throw new IllegalStateException("n-row bitmask assumes n≤31");
		double[][] U = cube.denseU(), V = cube.denseV(), W = cube.denseW();

		// Per-product n-row support masks (which A-rows / C-rows the product touches).
		int[] uMask = new int[r], wMask = new int[r];
		for (int k = 0; k < r; k++) {
			int um = 0, wm = 0;
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < m; j++) if (Math.abs(U[i * m + j][k]) > EPS) { um |= 1 << i; break; }
				for (int l = 0; l < p; l++) if (Math.abs(W[i * p + l][k]) > EPS) { wm |= 1 << i; break; }
			}
			uMask[k] = um; wMask[k] = wm;
		}

		// Baseline — best single coordinate quotient v=e_d (reproduces μ=114 → 10442).
		long bestCoord = Long.MAX_VALUE; int bestD = -1;
		for (int d = 0; d < n; d++) {
			int notd = ~(1 << d);
			long surv = 0;
			for (int k = 0; k < r; k++) if ((uMask[k] & notd) != 0 && (wMask[k] & notd) != 0) surv++;
			if (surv < bestCoord) { bestCoord = surv; bestD = d; }
		}
		log.info("baseline best COORDINATE quotient ⟨{},{},{}⟩ = {} (kill e_{})", n - 1, m, p, bestCoord, bestD);

		// Aggregation directions v = e_a + s·e_d, s ∈ {+1,−1}. A product dies iff its
		// A-form ∥ v (all rows ⊆ {a,d} with U[d,·]=s·U[a,·]) OR its C-form ∥ v.
		int[] best = {-1, -1, 0}; long bestDir = Long.MAX_VALUE;
		for (int a = 0; a < n; a++) {
			for (int d = a + 1; d < n; d++) {
				int notad = ~((1 << a) | (1 << d));
				for (int s : new int[] {1, -1}) {
					long surv = 0;
					for (int k = 0; k < r; k++) {
						boolean inDead = (uMask[k] & notad) == 0 && parallel(U, m, a, d, s, k);
						if (inDead) continue;
						boolean outDead = (wMask[k] & notad) == 0 && parallel(W, p, a, d, s, k);
						if (outDead) continue;
						surv++;
					}
					if (surv < bestDir) { bestDir = surv; best[0] = a; best[1] = d; best[2] = s; }
				}
			}
		}
		log.info("best AGGREGATION-direction quotient ⟨{},{},{}⟩ = {} (kill e_{} {} e_{})",
				n - 1, m, p, bestDir, best[0], best[2] > 0 ? "+" : "−", best[1]);
		log.info("SUMMARY ⟨{},{},{}⟩: coord-drop={}  best-direction={}  | FMM target=10413",
				n - 1, m, p, bestCoord, Math.min(bestCoord, bestDir));

		// Proof by construction: rebuild + verify the winning projection.
		if (bestDir < bestCoord) {
			NonCubicBilinearAlgorithm proj = buildProjection(cube, best[0], best[1], best[2], uMask, wMask);
			boolean ok = Verifier.passesRandomMatmulSpotCheck(proj, 20_000, 1e-9);
			log.info("VERIFIED winner ⟨{},{},{}⟩ r={} exact={}  (beats coord-drop {} by {})",
					proj.n, proj.m, proj.p, proj.r, ok, bestCoord, bestCoord - proj.r);
		} else {
			log.info("no aggregation direction beat the coordinate drop here.");
		}
	}

	/** True iff every column-vector of {@code F} (over the m/p axis) for product k is
	 *  ∥ (e_a + s·e_d), i.e. F[d,·] = s·F[a,·] (rows outside {a,d} already masked out). */
	private static boolean parallel(double[][] F, int q, int a, int d, int s, int k) {
		for (int j = 0; j < q; j++)
			if (Math.abs(F[d * q + j][k] - s * F[a * q + j][k]) > EPS) return false;
		return true;
	}

	/** Rebuild the projected ⟨n−1,m,p⟩ scheme that quotients direction v=e_a+s·e_d.
	 *  Basis of v⊥: {e_i : i∉{a,d}} ∪ {e_a − s·e_d} (the a-slot). E embeds it,
	 *  Q = projection onto v⊥ along v (so the a-slot output is (W_a − s·W_d)/2). */
	private static NonCubicBilinearAlgorithm buildProjection(
			NonCubicBilinearAlgorithm cube, int a, int d, int s, int[] uMask, int[] wMask) {
		int n = cube.n, m = cube.m, p = cube.p, r = cube.r, n2 = n - 1;
		int notad = ~((1 << a) | (1 << d));
		double[][] U = cube.denseU(), V = cube.denseV(), W = cube.denseW();
		// new n-row order = old rows except d; the slot that was 'a' carries direction e_a−s·e_d.
		int[] oldRow = new int[n2];
		for (int i = 0, w = 0; i < n; i++) if (i != d) oldRow[w++] = i;

		List<Integer> keep = new ArrayList<>();
		for (int k = 0; k < r; k++) {
			boolean inDead = (uMask[k] & notad) == 0 && parallel(U, m, a, d, s, k);
			if (inDead) continue;
			boolean outDead = (wMask[k] & notad) == 0 && parallel(W, p, a, d, s, k);
			if (outDead) continue;
			keep.add(k);
		}
		int r2 = keep.size();
		double[][] U2 = new double[n2 * m][r2], V2 = new double[m * p][r2], W2 = new double[n2 * p][r2];
		for (int kk = 0; kk < r2; kk++) {
			int k = keep.get(kk);
			for (int x = 0; x < m * p; x++) V2[x][kk] = V[x][k];
			for (int i2 = 0; i2 < n2; i2++) {
				int oi = oldRow[i2];
				boolean aSlot = oi == a;
				for (int j = 0; j < m; j++)
					U2[i2 * m + j][kk] = aSlot ? U[a * m + j][k] - s * U[d * m + j][k] : U[oi * m + j][k];
				for (int l = 0; l < p; l++)
					W2[i2 * p + l][kk] = aSlot ? (W[a * p + l][k] - s * W[d * p + l][k]) / 2.0 : W[oi * p + l][k];
			}
		}
		return new NonCubicBilinearAlgorithm(n2, m, p, U2, V2, W2);
	}
}
