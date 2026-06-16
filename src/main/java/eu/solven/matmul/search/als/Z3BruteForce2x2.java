package eu.solven.matmul.search.als;

import lombok.extern.slf4j.Slf4j;

import eu.solven.matmul.Verifier;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

/**
 * Exhaustive enumeration of Z/3-equivariant rank-7 decompositions of the
 * trilinear matmul tensor T(A,B,C) = trace(A·B·C) for n=2, over the alphabet
 * {-1, 0, +1}, with Strassen-shaped structure 1 + 2 (one fixed triple plus
 * two orbit generators).
 *
 * Why this matters: counting how many distinct canonical solutions exist tests
 * whether the symmetry quotient we've claimed in [[SOLVING_STRATEGIES]] §2.1
 * is actually the right one. De Groote 1978 proved Strassen is essentially
 * unique up to the full GL₂³ ⋊ S_3 gauge; over {-1, 0, +1} with only the
 * discrete part of that gauge quotiented out we expect a SMALL number of
 * canonical solutions — a large count would indicate we're missing a symmetry.
 *
 * Search space (per [[SOLVING_STRATEGIES]] §2.1):
 *   raw equivariant:          3^{28}  ≈  2.3 · 10^{13}
 *   after the canonicalization stack we apply here:  ~10^{7–8}
 *
 * Why T_trilin instead of T_matmul: the raw Z/3 cyclic action (u,v,w) → (v,w,u)
 * makes orbit contributions cyclic-symmetric in the tensor indices (a,b,c).
 * T_trilin IS cyclic-symmetric — T_matmul is not. A rank-r decomposition of
 * T_trilin converts to T_matmul by transposing the W factor.
 *
 * Canonicalization applied here:
 *   - u_fix: lex-min over (per-coordinate ±1 sign)  →  first nonzero entry positive
 *   - each orbit (u,v,w): cyclic-canonical  →  u lex-smallest of {u, v, w}
 *   - the two orbits: S_2-canonical  →  orbit1 ≤ orbit2 lex
 *
 * Canonicalization NOT applied (deliberately, to see the multiplicative effect):
 *   - per-term scaling beyond sign already absorbed (αβγ = 1 with α,β,γ ∈ {±1})
 *   - full monomial gauge (S_n ⋉ {±1}^n)^3 on the n² indices  →  factor 8³ = 512
 *     for n=2 — quotienting this is more involved (acts simultaneously on all
 *     vectors), so we leave it to a post-hoc dedupe of the printed solutions.
 *
 * Expected runtime: ~minutes on a single core thanks to the R2 diagonal prune.
 */
@Slf4j
public class Z3BruteForce2x2 {

	private static final int N = 2;
	private static final int N2 = 4;
	private static final int NUM_VECS = 81; // 3^4
	private static final int[][] VECS = new int[NUM_VECS][N2];
	static {
		// Index 0 → (-1,-1,-1,-1), index 80 → (+1,+1,+1,+1), lex order.
		for (int idx = 0; idx < NUM_VECS; idx++) {
			int x = idx;
			for (int i = N2 - 1; i >= 0; i--) {
				VECS[idx][i] = (x % 3) - 1;
				x /= 3;
			}
		}
	}

	/** All 27 triples (a,b,c) ∈ {-1,0,1}³ grouped by their product abc ∈ {-1, 0, +1}. */
	private static final int[][][] PER_INDEX_CANDS = buildPerIndexCands();

	private static int[][][] buildPerIndexCands() {
		// PER_INDEX_CANDS[p+1] = list of (a, b, c) triples with a*b*c == p.
		List<int[]> neg = new ArrayList<>();
		List<int[]> zero = new ArrayList<>();
		List<int[]> pos = new ArrayList<>();
		for (int a = -1; a <= 1; a++) {
			for (int b = -1; b <= 1; b++) {
				for (int c = -1; c <= 1; c++) {
					int prod = a * b * c;
					int[] triple = { a, b, c };
					if (prod == -1) neg.add(triple);
					else if (prod == 0) zero.add(triple);
					else pos.add(triple);
				}
			}
		}
		return new int[][][] { neg.toArray(new int[0][]), zero.toArray(new int[0][]), pos.toArray(new int[0][]) };
	}

	public static class Solution {
		public final int[] uFix;
		public final int[] u1, v1, w1;
		public final int[] u2, v2, w2;

		public Solution(int[] uFix, int[] u1, int[] v1, int[] w1, int[] u2, int[] v2, int[] w2) {
			this.uFix = uFix.clone();
			this.u1 = u1.clone(); this.v1 = v1.clone(); this.w1 = w1.clone();
			this.u2 = u2.clone(); this.v2 = v2.clone(); this.w2 = w2.clone();
		}

		@Override
		public String toString() {
			return String.format("u_fix=%s orbit1=(%s,%s,%s) orbit2=(%s,%s,%s)",
					fmt(uFix), fmt(u1), fmt(v1), fmt(w1), fmt(u2), fmt(v2), fmt(w2));
		}
	}

	public static List<Solution> enumerate(boolean printProgress) {
		int[][][] T = Verifier.intTrilinTensor(N);
		List<Solution> solutions = new ArrayList<>();

		long start = System.nanoTime();
		long checkedPairs = 0;
		long passedR2Diag = 0;
		long fullVerifications = 0;

		for (int uFixIdx = 0; uFixIdx < NUM_VECS; uFixIdx++) {
			int[] uFix = VECS[uFixIdx];
			if (!fixedCanonical(uFix)) continue;

			int[][][] R1 = subFixed(T, uFix);

			for (int u1Idx = 0; u1Idx < NUM_VECS; u1Idx++) {
				int[] u1 = VECS[u1Idx];
				for (int v1Idx = u1Idx; v1Idx < NUM_VECS; v1Idx++) {
					int[] v1 = VECS[v1Idx];
					for (int w1Idx = u1Idx; w1Idx < NUM_VECS; w1Idx++) {
						int[] w1 = VECS[w1Idx];

						int[][][] R2 = subOrbit(R1, u1, v1, w1);
						checkedPairs++;

						// R2 diagonal prune: R2[i,i,i] = 3·u2[i]·v2[i]·w2[i] with each factor
						// in {-1, 0, +1}, so R2[i,i,i] must be in {-3, 0, +3}. Exact int compare.
						int[] diagProd = new int[N2];
						boolean diagOk = true;
						for (int i = 0; i < N2; i++) {
							int d = R2[i][i][i];
							if (d == 3) diagProd[i] = 1;
							else if (d == 0) diagProd[i] = 0;
							else if (d == -3) diagProd[i] = -1;
							else { diagOk = false; break; }
						}
						if (!diagOk) continue;
						passedR2Diag++;

						// Per-index candidate enumeration. Each i has 4 (prod=±1) or 19 (prod=0)
						// (u2[i], v2[i], w2[i]) candidates. Iterate Cartesian product.
						int[][] c0 = PER_INDEX_CANDS[diagProd[0] + 1];
						int[][] c1 = PER_INDEX_CANDS[diagProd[1] + 1];
						int[][] c2 = PER_INDEX_CANDS[diagProd[2] + 1];
						int[][] c3 = PER_INDEX_CANDS[diagProd[3] + 1];

						int[] u2 = new int[N2];
						int[] v2 = new int[N2];
						int[] w2 = new int[N2];

						for (int[] t0 : c0) {
							u2[0] = t0[0]; v2[0] = t0[1]; w2[0] = t0[2];
							for (int[] t1 : c1) {
								u2[1] = t1[0]; v2[1] = t1[1]; w2[1] = t1[2];
								for (int[] t2 : c2) {
									u2[2] = t2[0]; v2[2] = t2[1]; w2[2] = t2[2];
									for (int[] t3 : c3) {
										u2[3] = t3[0]; v2[3] = t3[1]; w2[3] = t3[2];
										fullVerifications++;

										if (!orbitMatches(R2, u2, v2, w2)) continue;
										if (!orbitCyclicCanonical(u2, v2, w2)) continue;
										if (orbitSwapBreaks(u1, v1, w1, u2, v2, w2)) continue;

										solutions.add(new Solution(uFix, u1, v1, w1, u2, v2, w2));
									}
								}
							}
						}
					}
				}
			}
			if (printProgress) {
				long ms = (System.nanoTime() - start) / 1_000_000;
				log.info(String.format("u_fix=%s | pairs=%d, R2-diag-pass=%d, full-checks=%d, solutions=%d, %d ms%n",
						fmt(uFix), checkedPairs, passedR2Diag, fullVerifications, solutions.size(), ms));
			}
		}
		return solutions;
	}

	// ---- canonical filters ---------------------------------------------------

	/** First nonzero entry of u is positive (or u is all zeros). */
	static boolean fixedCanonical(int[] u) {
		for (int i = 0; i < N2; i++) {
			if (u[i] != 0) return u[i] > 0;
		}
		return true;
	}

	/** (u, v, w) is cyclic-canonical iff u is lex-smallest of {u, v, w}. */
	static boolean orbitCyclicCanonical(int[] u, int[] v, int[] w) {
		return lexCmp(u, v) <= 0 && lexCmp(u, w) <= 0;
	}

	/** Returns true iff orbit1 > orbit2 (i.e. should be skipped because we want orbit1 ≤ orbit2). */
	static boolean orbitSwapBreaks(int[] u1, int[] v1, int[] w1, int[] u2, int[] v2, int[] w2) {
		int c = lexCmp(u1, u2); if (c != 0) return c > 0;
		c = lexCmp(v1, v2); if (c != 0) return c > 0;
		return lexCmp(w1, w2) > 0;
	}

	static int lexCmp(int[] a, int[] b) {
		for (int i = 0; i < a.length; i++) {
			if (a[i] != b[i]) return Integer.compare(a[i], b[i]);
		}
		return 0;
	}

	// ---- tensor ops ----------------------------------------------------------

	static int[][][] subFixed(int[][][] T, int[] u) {
		int[][][] R = new int[N2][N2][N2];
		for (int a = 0; a < N2; a++)
			for (int b = 0; b < N2; b++)
				for (int c = 0; c < N2; c++)
					R[a][b][c] = T[a][b][c] - u[a] * u[b] * u[c];
		return R;
	}

	static int[][][] subOrbit(int[][][] R, int[] u, int[] v, int[] w) {
		int[][][] out = new int[N2][N2][N2];
		for (int a = 0; a < N2; a++)
			for (int b = 0; b < N2; b++)
				for (int c = 0; c < N2; c++)
					out[a][b][c] = R[a][b][c]
							- u[a] * v[b] * w[c]
							- v[a] * w[b] * u[c]
							- w[a] * u[b] * v[c];
		return out;
	}

	static boolean orbitMatches(int[][][] R, int[] u, int[] v, int[] w) {
		for (int a = 0; a < N2; a++) {
			for (int b = 0; b < N2; b++) {
				for (int c = 0; c < N2; c++) {
					int contrib = u[a] * v[b] * w[c] + v[a] * w[b] * u[c] + w[a] * u[b] * v[c];
					if (R[a][b][c] != contrib) return false;
				}
			}
		}
		return true;
	}

	// ---- in-loop canonicalization helpers ------------------------------------

	private static int[] negateInt(int[] u) {
		int[] out = new int[u.length];
		for (int i = 0; i < u.length; i++) out[i] = -u[i];
		return out;
	}

	private static int compareTriple(int[] u, int[] v, int[] w, int[] uR, int[] vR, int[] wR) {
		int c = lexCmp(u, uR); if (c != 0) return c;
		c = lexCmp(v, vR); if (c != 0) return c;
		return lexCmp(w, wR);
	}

	/**
	 * Conjugate a 4-vector u (representing 2×2 matrix M) by a diagonal monomial P:
	 *   new M[i,j] = s[i] · s[j] · M[σ_P(i), σ_P(j)].
	 * This is the action of the diagonal-gauge subgroup {(P, P, P)} on a single
	 * U/V/W slot. Preserves Z/3-equivariance and stabilizes u_fix = I.
	 */
	private static int[] conjMono(int[] u, int perm, int s0, int s1) {
		int p0 = perm == 0 ? 0 : 1;
		int p1 = perm == 0 ? 1 : 0;
		int[] s = { s0, s1 };
		int[] pP = { p0, p1 };
		int[] out = new int[N2];
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				out[i * N + j] = s[i] * s[j] * u[pP[i] * N + pP[j]];
			}
		}
		return out;
	}

	/**
	 * Returns true iff (u, v, w) is the lex-min representative of its orbit under
	 * per-term sign (αβγ=1, 4 elts) × cyclic-within-orbit (3 elts) — group size 12.
	 *
	 * Note: this does NOT include the diagonal monomial gauge. The gauge acts on
	 * the full decomposition (u_fix, orbit1, orbit2) simultaneously, so checking
	 * an orbit in isolation against the gauge would reject orbits whose canonical
	 * form depends on the paired orbit. The gauge stays in the post-hoc dedup.
	 */
	private static boolean isSignCyclicCanonical(int[] u, int[] v, int[] w) {
		final int[][] TERM_SIGNS = { { 1, 1, 1 }, { 1, -1, -1 }, { -1, 1, -1 }, { -1, -1, 1 } };
		for (int[] ts : TERM_SIGNS) {
			int[] uS = ts[0] == 1 ? u : negateInt(u);
			int[] vS = ts[1] == 1 ? v : negateInt(v);
			int[] wS = ts[2] == 1 ? w : negateInt(w);
			for (int rot = 0; rot < 3; rot++) {
				if (ts[0] == 1 && ts[1] == 1 && ts[2] == 1 && rot == 0) continue;
				int[] uR, vR, wR;
				if (rot == 0) { uR = uS; vR = vS; wR = wS; }
				else if (rot == 1) { uR = vS; vR = wS; wR = uS; }
				else { uR = wS; vR = uS; wR = vS; }
				if (compareTriple(uR, vR, wR, u, v, w) < 0) return false;
			}
		}
		return true;
	}

	/**
	 * (kept for reference — was over-restrictive when used in-loop because the
	 * diagonal gauge acts on the full decomposition, not on a single orbit.)
	 */
	@SuppressWarnings("unused")
	private static boolean isJointCanonical_DEPRECATED(int[] u, int[] v, int[] w) {
		final int[][] TERM_SIGNS = { { 1, 1, 1 }, { 1, -1, -1 }, { -1, 1, -1 }, { -1, -1, 1 } };
		final int[][] GAUGE_SIGNS = { { 1, 1 }, { 1, -1 }, { -1, 1 }, { -1, -1 } };
		for (int[] ts : TERM_SIGNS) {
			int[] uS = ts[0] == 1 ? u : negateInt(u);
			int[] vS = ts[1] == 1 ? v : negateInt(v);
			int[] wS = ts[2] == 1 ? w : negateInt(w);
			for (int rot = 0; rot < 3; rot++) {
				int[] uR, vR, wR;
				if (rot == 0) { uR = uS; vR = vS; wR = wS; }
				else if (rot == 1) { uR = vS; vR = wS; wR = uS; }
				else { uR = wS; vR = uS; wR = vS; }
				for (int gp = 0; gp < 2; gp++) {
					for (int[] gs : GAUGE_SIGNS) {
						if (ts[0] == 1 && ts[1] == 1 && ts[2] == 1
								&& rot == 0 && gp == 0 && gs[0] == 1 && gs[1] == 1) {
							continue; // skip identity
						}
						int[] uG = conjMono(uR, gp, gs[0], gs[1]);
						int[] vG = conjMono(vR, gp, gs[0], gs[1]);
						int[] wG = conjMono(wR, gp, gs[0], gs[1]);
						if (compareTriple(uG, vG, wG, u, v, w) < 0) return false;
					}
				}
			}
		}
		return true;
	}

	/**
	 * In-loop variant: applies isJointCanonical as a filter on both orbit1 and
	 * orbit2, in addition to all the filters used by {@link #enumerate}. Should
	 * produce ~1 solution directly (matching the post-hoc dedup output) and
	 * run substantially faster due to the inner-loop reduction.
	 */
	public static List<Solution> enumerateInLoopCanonical(boolean printProgress) {
		int[][][] T = Verifier.intTrilinTensor(N);
		List<Solution> solutions = new ArrayList<>();

		long start = System.nanoTime();

		for (int uFixIdx = 0; uFixIdx < NUM_VECS; uFixIdx++) {
			int[] uFix = VECS[uFixIdx];
			if (!fixedCanonical(uFix)) continue;
			int[][][] R1 = subFixed(T, uFix);
			for (int u1Idx = 0; u1Idx < NUM_VECS; u1Idx++) {
				solutions.addAll(enumerateForU1(uFix, R1, VECS[u1Idx], u1Idx));
			}
			if (printProgress) {
				long ms = (System.nanoTime() - start) / 1_000_000;
				log.info(String.format("[in-loop] u_fix=%s | solutions=%d, %d ms%n",
						fmt(uFix), solutions.size(), ms));
			}
		}
		return solutions;
	}

	/**
	 * Parallel variant: distributes (u_fix, u1) work units across a work-stealing
	 * thread pool. Each task computes its own R1, scans its v1/w1 sub-grid, and
	 * returns a local list of solutions; merging is a single sequential pass after
	 * all tasks complete.
	 *
	 * No shared mutable state in the hot path. Speedup is bounded by load
	 * balance: u_fix = (1, 0, 0, 1) carries almost all the work, so we split it
	 * along u1Idx to expose enough parallelism (~81 tasks just for that u_fix).
	 */
	public static List<Solution> enumerateInLoopCanonicalParallel(boolean printProgress) {
		final int[][][] T = Verifier.intTrilinTensor(N);

		// Pre-flatten canonical (uFixIdx, u1Idx) pairs into tasks.
		List<int[]> taskPairs = new ArrayList<>();
		for (int uFixIdx = 0; uFixIdx < NUM_VECS; uFixIdx++) {
			if (!fixedCanonical(VECS[uFixIdx])) continue;
			for (int u1Idx = 0; u1Idx < NUM_VECS; u1Idx++) {
				taskPairs.add(new int[] { uFixIdx, u1Idx });
			}
		}

		int parallelism = Runtime.getRuntime().availableProcessors();
		ExecutorService exec = Executors.newWorkStealingPool(parallelism);
		long start = System.nanoTime();
		try {
			List<Future<List<Solution>>> futures = new ArrayList<>(taskPairs.size());
			for (int[] pair : taskPairs) {
				final int[] uFix = VECS[pair[0]];
				final int u1Idx = pair[1];
				final int[] u1 = VECS[u1Idx];
				// Compute R1 inside the task — costs ~64 int ops, negligible vs the
				// per-task scan, and avoids sharing mutable state across tasks.
				futures.add(exec.submit(() -> enumerateForU1(uFix, subFixed(T, uFix), u1, u1Idx)));
			}
			List<Solution> solutions = new ArrayList<>();
			for (Future<List<Solution>> f : futures) {
				try {
					solutions.addAll(f.get());
				} catch (InterruptedException | ExecutionException e) {
					throw new RuntimeException(e);
				}
			}
			if (printProgress) {
				long ms = (System.nanoTime() - start) / 1_000_000;
				log.info(String.format("[parallel] %d tasks across %d threads | solutions=%d, %d ms%n",
						taskPairs.size(), parallelism, solutions.size(), ms));
			}
			return solutions;
		} finally {
			exec.shutdown();
		}
	}

	/**
	 * Per-(u_fix, u1) inner scan: iterates v1 ≥ u1, w1 ≥ u1 (cyclic-canonical),
	 * applies the R2-diagonal prune + per-index orbit2 enumeration with all the
	 * in-loop canonicalization filters. Pure function: reads its arguments and
	 * the static VECS / PER_INDEX_CANDS tables, returns a fresh list. Safe for
	 * concurrent invocation across distinct (u_fix, u1) pairs.
	 */
	private static List<Solution> enumerateForU1(int[] uFix, int[][][] R1, int[] u1, int u1Idx) {
		List<Solution> solutions = new ArrayList<>();
		for (int v1Idx = u1Idx; v1Idx < NUM_VECS; v1Idx++) {
			int[] v1 = VECS[v1Idx];
			for (int w1Idx = u1Idx; w1Idx < NUM_VECS; w1Idx++) {
				int[] w1 = VECS[w1Idx];

				if (!isSignCyclicCanonical(u1, v1, w1)) continue;

				int[][][] R2 = subOrbit(R1, u1, v1, w1);

				int[] diagProd = new int[N2];
				boolean diagOk = true;
				for (int i = 0; i < N2; i++) {
					int d = R2[i][i][i];
					if (d == 3) diagProd[i] = 1;
					else if (d == 0) diagProd[i] = 0;
					else if (d == -3) diagProd[i] = -1;
					else { diagOk = false; break; }
				}
				if (!diagOk) continue;

				int[][] c0 = PER_INDEX_CANDS[diagProd[0] + 1];
				int[][] c1 = PER_INDEX_CANDS[diagProd[1] + 1];
				int[][] c2 = PER_INDEX_CANDS[diagProd[2] + 1];
				int[][] c3 = PER_INDEX_CANDS[diagProd[3] + 1];

				int[] u2 = new int[N2];
				int[] v2 = new int[N2];
				int[] w2 = new int[N2];

				for (int[] t0 : c0) {
					u2[0] = t0[0]; v2[0] = t0[1]; w2[0] = t0[2];
					for (int[] t1 : c1) {
						u2[1] = t1[0]; v2[1] = t1[1]; w2[1] = t1[2];
						for (int[] t2 : c2) {
							u2[2] = t2[0]; v2[2] = t2[1]; w2[2] = t2[2];
							for (int[] t3 : c3) {
								u2[3] = t3[0]; v2[3] = t3[1]; w2[3] = t3[2];

								if (!orbitMatches(R2, u2, v2, w2)) continue;
								if (!orbitCyclicCanonical(u2, v2, w2)) continue;
								if (!isSignCyclicCanonical(u2, v2, w2)) continue;
								if (orbitSwapBreaks(u1, v1, w1, u2, v2, w2)) continue;

								solutions.add(new Solution(uFix, u1, v1, w1, u2, v2, w2));
							}
						}
					}
				}
			}
		}
		return solutions;
	}

	// ---- canonicalization over the monomial gauge × transpose × per-term signs ----

	private static final int[][] MONO_PERMS = { { 0, 1 }, { 1, 0 } };
	private static final int[][] SIGNS = { { 1, 1 }, { 1, -1 }, { -1, 1 }, { -1, -1 } };
	private static final int[][] TERM_SIGNS = { { 1, 1, 1 }, { 1, -1, -1 }, { -1, 1, -1 }, { -1, -1, 1 } };
	private static final int[][] GAUGES = buildGauges();

	private static int[][] buildGauges() {
		List<int[]> out = new ArrayList<>();
		for (int permP = 0; permP < 2; permP++) {
			for (int[] sP : SIGNS) {
				for (int permQ = 0; permQ < 2; permQ++) {
					for (int[] sQ : SIGNS) {
						for (int permR = 0; permR < 2; permR++) {
							for (int[] sR : SIGNS) {
								out.add(new int[] { permP, sP[0], sP[1], permQ, sQ[0], sQ[1], permR, sR[0], sR[1] });
							}
						}
					}
				}
			}
		}
		return out.toArray(new int[0][]);
	}

	/**
	 * Apply A → P A Q⁻¹ to a flattened 2×2 vector u (P, Q monomial).
	 * Derived: new u[i·n + j] = sign_P[i] · sign_Q[j] · u[σ_P(i)·n + σ_Q(j)].
	 * Valid for n=2 because monomial perms are self-inverse there.
	 */
	static int[] applyMono(int[] u, int permP, int sP0, int sP1, int permQ, int sQ0, int sQ1) {
		int[] pP = MONO_PERMS[permP], pQ = MONO_PERMS[permQ];
		int[] sP = { sP0, sP1 }, sQ = { sQ0, sQ1 };
		int[] out = new int[N2];
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				out[i * N + j] = sP[i] * sQ[j] * u[pP[i] * N + pQ[j]];
			}
		}
		return out;
	}

	/** Transpose action on a 4-vector representing a 2×2 matrix: swap entries 1 and 2. */
	static int[] transposeVec(int[] u) {
		return new int[] { u[0], u[2], u[1], u[3] };
	}

	static int[] negate(int[] u) {
		int[] out = new int[u.length];
		for (int i = 0; i < u.length; i++) out[i] = -u[i];
		return out;
	}

	/**
	 * Per-term sign normalization: pick (α, β, γ) ∈ {±1}³ with αβγ=1 to make
	 * (α·u, β·v, γ·w) lex-smallest. Returns the flattened normalized triple.
	 */
	static int[] normalizeTerm(int[] u, int[] v, int[] w) {
		int[] best = null;
		for (int[] s : TERM_SIGNS) {
			int[] uS = (s[0] == 1) ? u : negate(u);
			int[] vS = (s[1] == 1) ? v : negate(v);
			int[] wS = (s[2] == 1) ? w : negate(w);
			int[] flat = new int[3 * N2];
			System.arraycopy(uS, 0, flat, 0, N2);
			System.arraycopy(vS, 0, flat, N2, N2);
			System.arraycopy(wS, 0, flat, 2 * N2, N2);
			if (best == null || lexCmpArr(flat, best) < 0) best = flat;
		}
		return best;
	}

	static int lexCmpArr(int[] a, int[] b) {
		for (int i = 0; i < a.length; i++) {
			if (a[i] != b[i]) return Integer.compare(a[i], b[i]);
		}
		return 0;
	}

	/**
	 * Compute the canonical form of a solution: enumerate all 1024 elements of
	 * (monomial gauge)³ × transpose, transform the 7 rank-1 terms, sign-normalize
	 * each term, sort the multiset, and take the lex-min serialized form.
	 *
	 * Two raw solutions are equivalent iff they share the same canonical form.
	 */
	static String canonicalForm(Solution sol) {
		// Expand to 7 (u, v, w) rank-1 terms: 1 fixed + 3 cyclic of orbit1 + 3 cyclic of orbit2.
		int[][] U = { sol.uFix, sol.u1, sol.v1, sol.w1, sol.u2, sol.v2, sol.w2 };
		int[][] V = { sol.uFix, sol.v1, sol.w1, sol.u1, sol.v2, sol.w2, sol.u2 };
		int[][] W = { sol.uFix, sol.w1, sol.u1, sol.v1, sol.w2, sol.u2, sol.v2 };

		String best = null;
		for (int[] g : GAUGES) {
			int permP = g[0], sP0 = g[1], sP1 = g[2];
			int permQ = g[3], sQ0 = g[4], sQ1 = g[5];
			int permR = g[6], sR0 = g[7], sR1 = g[8];
			for (int tr = 0; tr < 2; tr++) {
				int[][] flats = new int[7][];
				for (int t = 0; t < 7; t++) {
					int[] u = applyMono(U[t], permP, sP0, sP1, permQ, sQ0, sQ1);
					int[] v = applyMono(V[t], permQ, sQ0, sQ1, permR, sR0, sR1);
					int[] w = applyMono(W[t], permR, sR0, sR1, permP, sP0, sP1);
					if (tr == 1) {
						// Transpose Z/2 action on a rank-1 term: (u, v, w) → (T(w), T(v), T(u)).
						int[] newU = transposeVec(w);
						int[] newV = transposeVec(v);
						int[] newW = transposeVec(u);
						u = newU; v = newV; w = newW;
					}
					flats[t] = normalizeTerm(u, v, w);
				}
				Arrays.sort(flats, Z3BruteForce2x2::lexCmpArr);
				StringBuilder sb = new StringBuilder();
				for (int[] f : flats) sb.append(Arrays.toString(f)).append('|');
				String key = sb.toString();
				if (best == null || key.compareTo(best) < 0) best = key;
			}
		}
		return best;
	}

	/** Group solutions by canonical form. */
	static Map<String, List<Solution>> dedupeByCanonicalForm(List<Solution> solutions) {
		Map<String, List<Solution>> groups = new HashMap<>();
		for (Solution s : solutions) {
			String canon = canonicalForm(s);
			groups.computeIfAbsent(canon, k -> new ArrayList<>()).add(s);
		}
		return groups;
	}

	static String fmt(int[] v) {
		StringBuilder sb = new StringBuilder("[");
		for (int i = 0; i < v.length; i++) {
			if (i > 0) sb.append(",");
			sb.append(v[i] >= 0 ? " " : "");
			sb.append(v[i]);
		}
		return sb.append("]").toString();
	}

	public static void main(String[] args) {
		int cores = Runtime.getRuntime().availableProcessors();
		log.info(String.format("Brute-force enumeration of Z/3-equivariant {-1,0,+1} rank-7 algorithms%n"
				+ "for the 2×2 trilinear matmul tensor, using %d cores.%n%n", cores));

		long t = System.nanoTime();
		List<Solution> solutions = enumerateInLoopCanonicalParallel(true);
		long ms = (System.nanoTime() - t) / 1_000_000;
		log.info(String.format("%nFound %d partial-canonical solutions in %d ms.%n", solutions.size(), ms));

		Map<String, List<Solution>> canonical = dedupeByCanonicalForm(solutions);
		log.info(String.format("After full-symmetry dedup: %d distinct canonical solution(s).%n%n",
				canonical.size()));

		int i = 0;
		for (Map.Entry<String, List<Solution>> e : canonical.entrySet()) {
			i++;
			log.info(String.format("Equivalence class %d (orbit size %d): %s%n",
					i, e.getValue().size(), e.getValue().get(0)));
		}
	}
}
