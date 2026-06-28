package eu.solven.matmul;

import java.math.BigInteger;

/**
 * Verifies a bilinear algorithm against the exact matrix-multiplication tensor.
 *
 *   T_exact[(a,b), (c,d), (i,j)] = δ(a=i) · δ(b=c) · δ(d=j)
 *   T_approx[a', b', c']         = Σ_k U[a'][k] · V[b'][k] · W[c'][k]
 *
 * Returns the Frobenius residual ‖T_exact − T_approx‖. Residual = 0 ⇔ the
 * algorithm is exact.
 */
public class Verifier {

	public static double[][][] matmulTensor(int n) {
		int n2 = n * n;
		double[][][] T = new double[n2][n2][n2];
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				for (int l = 0; l < n; l++) {
					T[i * n + l][l * n + j][i * n + j] = 1.0;
				}
			}
		}
		return T;
	}

	/**
	 * Trilinear-form matmul tensor: T[a,b,c] = 1 iff there exist i, j, k with
	 * a = i·n + j, b = j·n + k, c = k·n + i — i.e. the (a,b,c) entry of
	 * trace(A·B·C). Unlike {@link #matmulTensor}, this tensor IS cyclic-symmetric
	 * in (a,b,c) — trace(ABC) = trace(BCA) = trace(CAB) — which is the form
	 * needed for raw-Z/3-equivariant search.
	 *
	 * Any rank-r decomposition of {@code trilinTensor(n)} converts to a rank-r
	 * decomposition of {@code matmulTensor(n)} by transposing the W-factor
	 * (swapping the row/column flatten convention in W's matrix interpretation).
	 */
	public static double[][][] trilinTensor(int n) {
		int n2 = n * n;
		double[][][] T = new double[n2][n2][n2];
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				for (int k = 0; k < n; k++) {
					T[i * n + j][j * n + k][k * n + i] = 1.0;
				}
			}
		}
		return T;
	}

	/**
	 * Integer version of {@link #matmulTensor(int)} — every entry is exactly 0 or 1.
	 * Use this for discrete-alphabet search code (e.g. brute force over {-1, 0, +1})
	 * where the residual tensor stays in a small integer range and exact comparison
	 * is both correct and faster than floating-point with epsilon checks.
	 */
	public static int[][][] intMatmulTensor(int n) {
		int n2 = n * n;
		int[][][] T = new int[n2][n2][n2];
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				for (int l = 0; l < n; l++) {
					T[i * n + l][l * n + j][i * n + j] = 1;
				}
			}
		}
		return T;
	}

	/**
	 * Non-cubic version of {@link #intMatmulTensor(int)} for {@code ⟨n, m, p⟩}:
	 * compute {@code C = A·B} where A is {@code n×m}, B is {@code m×p}, C is
	 * {@code n×p}. Returns the integer tensor {@code T[a, b, c] = 1} iff there
	 * exist {@code (i, l, j)} with {@code a = i·m + l, b = l·p + j, c = i·p + j}.
	 * Shape: {@code int[n·m][m·p][n·p]}.
	 */
	public static int[][][] intMatmulTensor(int n, int m, int p) {
		int dimU = n * m, dimV = m * p, dimW = n * p;
		int[][][] T = new int[dimU][dimV][dimW];
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < p; j++) {
				for (int l = 0; l < m; l++) {
					T[i * m + l][l * p + j][i * p + j] = 1;
				}
			}
		}
		return T;
	}

	/**
	 * O(1) single-entry of the {@code ⟨n,m,p⟩} matmul tensor — the value
	 * {@link #intMatmulTensor(int, int, int)} would store at {@code T[a][b][c]},
	 * computed by decoding the flattened indices instead of materialising the
	 * (mostly-zero) {@code dimU×dimV×dimW} array. The dense tensor is
	 * {@code (n·m)·(m·p)·(n·p)} ints — 4.3&nbsp;GB at ⟨32,32,32⟩ — for only
	 * {@code n·m·p} ones, so the sampled verifiers (which read a handful of
	 * positions) MUST use this, not the array.
	 *
	 * <p>{@code a = i·m+l} (A entry), {@code b = l′·p+j} (B entry),
	 * {@code c = i″·p+j″} (C entry); the tensor is 1 iff {@code i==i″},
	 * {@code l==l′}, {@code j==j″} (the single contraction A[i,l]·B[l,j]→C[i,j]).</p>
	 */
	public static int matmulTensorEntry(int a, int b, int c, int n, int m, int p) {
		int i = a / m, l = a % m;        // a = i·m + l
		int lB = b / p, j = b % p;       // b = l·p + j
		int iC = c / p, jC = c % p;      // c = i·p + j
		return (i == iC && l == lB && j == jC) ? 1 : 0;
	}

	/** Integer version of {@link #trilinTensor(int)}. See {@link #intMatmulTensor(int)} for why. */
	public static int[][][] intTrilinTensor(int n) {
		int n2 = n * n;
		int[][][] T = new int[n2][n2][n2];
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				for (int k = 0; k < n; k++) {
					T[i * n + j][j * n + k][k * n + i] = 1;
				}
			}
		}
		return T;
	}

	public static double residual(BilinearAlgorithm alg) {
		return residualAgainst(alg, matmulTensor(alg.n));
	}

	/** Frobenius residual {@code ‖T − Σ_k u_k ⊗ v_k ⊗ w_k‖} against any target tensor. */
	public static double residualAgainst(BilinearAlgorithm alg, double[][][] target) {
		int n2 = alg.n * alg.n;
		double sumSq = 0.0;
		for (int a = 0; a < n2; a++) {
			for (int b = 0; b < n2; b++) {
				for (int c = 0; c < n2; c++) {
					double approx = 0.0;
					for (int k = 0; k < alg.r; k++) {
						approx += alg.U[a][k] * alg.V[b][k] * alg.W[c][k];
					}
					double diff = target[a][b][c] - approx;
					sumSq += diff * diff;
				}
			}
		}
		return Math.sqrt(sumSq);
	}

	public static boolean isExact(BilinearAlgorithm alg) {
		return residual(alg) < 1e-10;
	}

	/**
	 * Visit every {@code (a, b, c)} term in the matmul-tensor check loop,
	 * supplying the actual algorithm-computed coefficient and the
	 * expected integer target. Return {@code false} from {@link #accept}
	 * to short-circuit iteration (used by symbolic-diagnostic callers
	 * that cap their output).
	 */
	@FunctionalInterface
	public interface TermVisitor {
		/** @return true to continue, false to short-circuit. */
		boolean accept(int a, int b, int c, int expected, double actual);
	}

	/**
	 * Iterate every {@code (a, b, c)} term of the matmul-tensor check
	 * for a non-cubic algorithm. Computes {@code actual} via the
	 * rank-{@code r} bilinear expansion and looks up the integer
	 * target from {@link #intMatmulTensor}. Shared driver for
	 * {@link #residualNonCubic} (full sum) and
	 * {@link #symbolicDiscrepancies} (term-level diagnostics with
	 * early-exit).
	 */
	public static void forEachTerm(NonCubicBilinearAlgorithm alg, TermVisitor visitor) {
		double[][] srcU = alg.denseU(), srcV = alg.denseV(), srcW = alg.denseW();
		int dimU = alg.n * alg.m, dimV = alg.m * alg.p, dimW = alg.n * alg.p;
		int[][][] target = intMatmulTensor(alg.n, alg.m, alg.p);
		for (int a = 0; a < dimU; a++) {
			for (int b = 0; b < dimV; b++) {
				for (int c = 0; c < dimW; c++) {
					double actual = 0.0;
					for (int k = 0; k < alg.r; k++) {
						actual += srcU[a][k] * srcV[b][k] * srcW[c][k];
					}
					if (!visitor.accept(a, b, c, target[a][b][c], actual)) return;
				}
			}
		}
	}

	/** Non-cubic Frobenius residual against the {@code ⟨n,m,p⟩} matmul tensor. */
	public static double residualNonCubic(NonCubicBilinearAlgorithm alg) {
		double[] sumSq = { 0.0 };
		forEachTerm(alg, (a, b, c, expected, actual) -> {
			double diff = expected - actual;
			sumSq[0] += diff * diff;
			return true;
		});
		return Math.sqrt(sumSq[0]);
	}

	/**
	 * Exact non-cubic matmul check — SPARSE + symbolic-target. Replaces the old dense
	 * path ({@link #forEachTerm} → {@code denseU/V/W} + a dense {@code (nm)×(mp)×(np)}
	 * {@link #intMatmulTensor}), which OOM'd on large shapes (the ⟨30,32,32⟩ tensor alone
	 * is ~943M ints ≈ 3.8 GB; the heap dump pinned it to {@code Verifier.forEachTerm}).
	 *
	 * <p>Instead: accumulate each rank-1 term's contribution over the factors' SPARSE
	 * columns into a {@code (a,b,c) → coeff} map (memory ∝ distinct product terms, not the
	 * dense cube), and compare each accumulated coefficient against the matmul target,
	 * which is decoded in O(1) ({@link #matmulTargetOf}) — no dense tensor. Coefficients
	 * are summed exactly for integer / dyadic-rational schemes (the common case: ±1, ½, ¼
	 * are exact in IEEE-754, so no rounding); the per-term comparison is against the exact
	 * integer target (0/1), not a Frobenius double residual. Non-dyadic rationals (e.g.
	 * 1/3) remain ε-approximate here — full BigRational symbolic accumulation is a
	 * follow-up if a non-dyadic scheme ever needs certifying.</p>
	 */
	public static boolean isExactNonCubic(NonCubicBilinearAlgorithm alg) {
		return exactDiscrepancyNonCubic(alg, 1e-9) == null;
	}

	/** The matmul-tensor entry for flattened factor indices: {@code U}-index {@code a=i·m+j},
	 *  {@code V}-index {@code b=j'·p+l}, {@code W}-index {@code c=i'·p+l'}. The matmul tensor
	 *  is 1 iff {@code j==j' ∧ i==i' ∧ l==l'}, else 0 — decoded in O(1), no dense allocation. */
	private static int matmulTargetOf(int a, int b, int c, int m, int p) {
		int i = a / m, j = a % m;       // U row/col
		int j2 = b / p, l = b % p;      // V row/col
		int i2 = c / p, l2 = c % p;     // W row/col
		return (j == j2 && i == i2 && l == l2) ? 1 : 0;
	}

	/** A sparse factor column: parallel (row, value) arrays of the non-zeros. */
	private record Col(int[] rows, double[] vals) {}

	private static Col col(FactorMatrix f, int k) {
		java.util.ArrayList<Integer> rs = new java.util.ArrayList<>();
		java.util.ArrayList<Double> vs = new java.util.ArrayList<>();
		f.forEachInColumn(k, (row, val) -> {
			if (val != 0.0) { rs.add(row); vs.add(val); }
		});
		int[] rows = new int[rs.size()];
		double[] vals = new double[vs.size()];
		for (int i = 0; i < rows.length; i++) { rows[i] = rs.get(i); vals[i] = vs.get(i); }
		return new Col(rows, vals);
	}

	/** Pack a flattened-index triple into one long. Each index &lt; 2²¹ (covers nm,mp,np
	 *  up to ~2M — far beyond any catalogued shape, where max is 32·32=1024). */
	private static long packTerm(long a, long b, long c) {
		return (a << 42) | (b << 21) | c;
	}

	/**
	 * First exact discrepancy of a non-cubic matmul scheme, or {@code null} if exact.
	 * SPARSE + SYMBOLIC: every coefficient is scaled to an exact integer by the common
	 * denominator {@code D} (lcm of all denominators) and accumulated as {@link BigInteger}
	 * — no doubles in the decision, and no dense factor/tensor. Each matmul-target term
	 * must equal {@code D³}; every other accumulated term must cancel to exactly 0.
	 * O(Σ_k |supp Uₖ|·|supp Vₖ|·|supp Wₖ|) time, O(#distinct terms) memory. The {@code eps}
	 * argument is retained for signature compatibility and ignored (the check is exact).
	 */
	public static String exactDiscrepancyNonCubic(NonCubicBilinearAlgorithm alg, double eps) {
		final int n = alg.n, m = alg.m, p = alg.p, r = alg.r;
		FactorMatrix uMat = alg.u(), vMat = alg.v(), wMat = alg.w();
		// Pass 1 (sparse): common denominator D = lcm of every non-zero coefficient's denom.
		BigInteger D = BigInteger.ONE;
		for (FactorMatrix f : new FactorMatrix[] { uMat, vMat, wMat }) {
			for (int k = 0; k < r; k++) {
				for (double v : col(f, k).vals) {
					D = SymbolicVerifier.lcm(D, SymbolicVerifier.denominatorOf(v));
				}
			}
		}
		BigInteger D3 = D.multiply(D).multiply(D);
		// Pass 2 (sparse): accumulate each rank-1 term's outer product as scaled integers.
		java.util.HashMap<Long, BigInteger> acc = new java.util.HashMap<>();
		for (int k = 0; k < r; k++) {
			Col cu = col(uMat, k), cv = col(vMat, k), cw = col(wMat, k);
			BigInteger[] uN = scaleInts(cu.vals, D), vN = scaleInts(cv.vals, D), wN = scaleInts(cw.vals, D);
			for (int ai = 0; ai < cu.rows.length; ai++) {
				int a = cu.rows[ai];
				for (int bi = 0; bi < cv.rows.length; bi++) {
					BigInteger uv = uN[ai].multiply(vN[bi]);
					if (uv.signum() == 0) continue;
					int b = cv.rows[bi];
					for (int ci = 0; ci < cw.rows.length; ci++) {
						BigInteger term = uv.multiply(wN[ci]);
						if (term.signum() == 0) continue;
						acc.merge(packTerm(a, b, cw.rows[ci]), term, BigInteger::add);
					}
				}
			}
		}
		// Every matmul-target triple (i,j,l) must accumulate to exactly D³; remove as we go.
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				for (int l = 0; l < p; l++) {
					BigInteger got = acc.remove(packTerm(i * m + j, j * p + l, i * p + l));
					if (got == null) got = BigInteger.ZERO;
					if (!got.equals(D3)) {
						return String.format("⟨%d,%d,%d⟩ matmul term (i=%d,j=%d,l=%d) = %s, want D³=%s",
								n, m, p, i, j, l, got, D3);
					}
				}
			}
		}
		// Anything left must cancel to EXACTLY 0.
		for (java.util.Map.Entry<Long, BigInteger> e : acc.entrySet()) {
			if (e.getValue().signum() != 0) {
				long key = e.getKey();
				long a = (key >> 42) & 0x1FFFFF, b = (key >> 21) & 0x1FFFFF, c = key & 0x1FFFFF;
				return String.format("⟨%d,%d,%d⟩ spurious term (a=%d,b=%d,c=%d) = %s, want 0",
						n, m, p, a, b, c, e.getValue());
			}
		}
		return null;
	}

	private static BigInteger[] scaleInts(double[] vals, BigInteger D) {
		BigInteger[] out = new BigInteger[vals.length];
		for (int i = 0; i < vals.length; i++) {
			out[i] = SymbolicVerifier.numeratorScaled(vals[i], D);
		}
		return out;
	}

	// ── size/density-aware dispatcher ───────────────────────────────────────────

	/** Which concrete verifier {@link #verifyAuto} chose. */
	public enum VerifyStrategy {
		/** {@link #isExactNonCubic} — a real algebraic proof; memory ∝ distinct product terms. */
		EXACT_SYMBOLIC,
		/** {@link #passesRandomMatmulSpotCheck} — O(dim) memory, randomised (false-accept ≈ 0). */
		RANDOM_SPOT_CHECK
	}

	/** Outcome of {@link #verifyAuto}: the verdict, WHICH strategy produced it, and the
	 *  estimated exact-verifier work that drove the choice (so callers/logs can be honest
	 *  about whether they got a proof or a spot-check). */
	public record Verdict(boolean ok, VerifyStrategy strategy, long estimatedTerms) {
		/** True iff the verdict is a genuine algebraic proof (not the randomised path). */
		public boolean isProof() { return strategy == VerifyStrategy.EXACT_SYMBOLIC; }
	}

	/**
	 * Default ceiling on the exact verifier's generated-term count (≈ both its time in
	 * BigInteger multiplies AND a loose upper bound on its HashMap memory). Above this we
	 * fall back to the spot-check. 30M terms ≈ a few seconds and &lt;~2 GB worst case;
	 * the ⟨30,32,32⟩ AlphaEvolve-based scheme (dense, rank 14863) estimates far above it,
	 * which is exactly the case that OOM'd the exact path at 12 GB.
	 */
	public static final long DEFAULT_MAX_EXACT_TERMS = 30_000_000L;

	/**
	 * Verify a non-cubic scheme with the STRONGEST AFFORDABLE strategy, chosen by size and
	 * density — not a fixed verifier. Sparse / small schemes get the exact symbolic proof
	 * ({@link #isExactNonCubic}); dense / large schemes (where the exact term-map would blow
	 * up — the AlphaEvolve-based ⟨30,32,32⟩ being the canonical OOM) get the memory-light
	 * random matmul spot-check ({@link #passesRandomMatmulSpotCheck}). The specific verifiers
	 * stay first-class — this only routes between them.
	 *
	 * <p>The gate is {@link #estimateExactTerms}: the sum over products of
	 * {@code |supp Uₖ|·|supp Vₖ|·|supp Wₖ|}, which is what the exact accumulator iterates
	 * (and an upper bound on the distinct terms it stores). Cheap to compute (O(nnz)).</p>
	 */
	public static Verdict verifyAuto(NonCubicBilinearAlgorithm alg, long maxExactTerms) {
		long est = estimateExactTerms(alg, maxExactTerms);
		if (est <= maxExactTerms) {
			return new Verdict(isExactNonCubic(alg), VerifyStrategy.EXACT_SYMBOLIC, est);
		}
		return new Verdict(passesRandomMatmulSpotCheck(alg), VerifyStrategy.RANDOM_SPOT_CHECK, est);
	}

	/** {@link #verifyAuto(NonCubicBilinearAlgorithm, long)} with {@link #DEFAULT_MAX_EXACT_TERMS}. */
	public static Verdict verifyAuto(NonCubicBilinearAlgorithm alg) {
		return verifyAuto(alg, DEFAULT_MAX_EXACT_TERMS);
	}

	/**
	 * Estimate the exact verifier's generated-term count: {@code Σ_k nnz(Uₖ)·nnz(Vₖ)·nnz(Wₖ)}.
	 * Counts non-zeros per column over the sparse factors (O(total nnz), no dense). Short-circuits
	 * at {@code cap}: returns {@link Long#MAX_VALUE} as soon as the running sum exceeds it (and on
	 * any overflow), so a hugely-dense scheme is classified without summing the whole thing.
	 */
	public static long estimateExactTerms(NonCubicBilinearAlgorithm alg, long cap) {
		int r = alg.r;
		int[] nu = colNonZeros(alg.u(), r), nv = colNonZeros(alg.v(), r), nw = colNonZeros(alg.w(), r);
		long total = 0;
		for (int k = 0; k < r; k++) {
			long terms = (long) nu[k] * nv[k] * nw[k];
			total += terms;
			if (total < 0 || total > cap) return Long.MAX_VALUE; // overflow or clearly over budget
		}
		return total;
	}

	/** Non-zero count per column of a sparse factor (index k → nnz of column k). */
	private static int[] colNonZeros(FactorMatrix f, int r) {
		int[] nz = new int[r];
		for (int k = 0; k < r; k++) {
			int[] c = { 0 };
			f.forEachInColumn(k, (row, val) -> { if (val != 0.0) c[0]++; });
			nz[k] = c[0];
		}
		return nz;
	}

	/**
	 * <strong>Fast randomised verification</strong> for non-cubic bilinear
	 * algorithms: instead of computing the full Frobenius residual against
	 * the matmul tensor (O(n²m²p²·r) — ~10¹¹ ops for ⟨21,21,21⟩),
	 * execute the algorithm on {@code samples} random {@code A, B} pairs,
	 * compute the naive {@code C = A·B} for comparison, and check the
	 * Frobenius norm of {@code |C_algo − C_naive|} per sample.
	 *
	 * <p>Complexity: O(samples · r · (nm + mp + np)) — for ⟨21,21,21⟩=5258
	 * with samples=10, that's ~10⁸ ops vs ~10¹¹ for the tensor verifier
	 * (≈ 1000× faster). False acceptance probability is essentially zero
	 * for any non-trivial discrepancy because random A,B sample the full
	 * polynomial identity densely.</p>
	 *
	 * @param alg     scheme to verify
	 * @param samples number of random (A, B) trials (≥ 1; 5–10 is plenty)
	 * @param epsilon per-trial Frobenius tolerance; pick relative to
	 *                {@code 10·n·m·p · alg.r · max|coefficient|² · sqrt(n·p)}
	 *                or use the default in {@link #passesRandomMatmulSpotCheck}
	 * @return true iff every trial's residual is below {@code epsilon}
	 */
	public static boolean passesRandomMatmulSpotCheck(NonCubicBilinearAlgorithm alg,
			int samples, double epsilon) {
		java.util.Random rng = new java.util.Random(0xC0DE);
		int n = alg.n, m = alg.m, p = alg.p, r = alg.r;
		FactorMatrix U = alg.u(), V = alg.v(), W = alg.w();
		for (int trial = 0; trial < samples; trial++) {
			// Flattened A (i·m+j) and B (j·p+l) — same fill order as the factor-matrix
			// row indexing, so the RNG sequence (seed 0xC0DE) is unchanged.
			double[] aFlat = new double[n * m];
			double[] bFlat = new double[m * p];
			for (int x = 0; x < aFlat.length; x++) aFlat[x] = rng.nextGaussian();
			for (int x = 0; x < bFlat.length; x++) bFlat[x] = rng.nextGaussian();

			// α_k = Σ U[:,k]·A, β_k = Σ V[:,k]·B — allocation-free sparse dot products.
			double[] alpha = new double[r];
			double[] beta  = new double[r];
			for (int k = 0; k < r; k++) {
				alpha[k] = U.dotColumn(k, aFlat);
				beta[k]  = V.dotColumn(k, bFlat);
			}

			// C[i·p+l] = Σ_k W[:,k]·γ_k — sparse AXPY accumulation.
			double[] cAlgo = new double[n * p];
			for (int k = 0; k < r; k++) {
				W.axpyColumn(k, alpha[k] * beta[k], cAlgo);
			}

			// Naive matmul C[i,l] = Σ_j A[i,j] · B[j,l].
			double sumSq = 0;
			for (int i = 0; i < n; i++) {
				for (int l = 0; l < p; l++) {
					double c = 0;
					for (int j = 0; j < m; j++) c += aFlat[i * m + j] * bFlat[j * p + l];
					double diff = cAlgo[i * p + l] - c;
					sumSq += diff * diff;
				}
			}
			if (Math.sqrt(sumSq) > epsilon) return false;
		}
		return true;
	}

	/** Convenience: 5 samples, ε scaled to algorithm size for safe default. */
	public static boolean passesRandomMatmulSpotCheck(NonCubicBilinearAlgorithm alg) {
		// Heuristic ε: ~1e-9 per (n·p) entry, scaled by sqrt(r) for accumulated rounding.
		double eps = 1e-9 * Math.sqrt(alg.n * alg.p) * Math.sqrt(alg.r);
		return passesRandomMatmulSpotCheck(alg, 5, eps);
	}

	/**
	 * F_p analogue of {@link #passesRandomMatmulSpotCheck}: fast randomised
	 * verification over the prime field GF({@code p}) — for F₂/F₃ schemes where the
	 * char-0 (random-real) spot-check would WRONGLY reject a valid scheme (e.g.
	 * AlphaTensor ⟨4,4,4⟩=47 computes matmul only mod 2). The factor matrices are
	 * reduced mod {@code p} ({@link #reduceFactorModP} — false if any denominator is
	 * divisible by {@code p}, i.e. the scheme is not representable mod {@code p}), then
	 * the algorithm is run on random {@code A,B ∈ GF(p)} and compared to the naive
	 * product mod {@code p}. Exact integer arithmetic — no tolerance. Same
	 * O(samples·r·(nm+mp+np)) cost as the char-0 path, so feasible at large shapes
	 * (the full-tensor {@link #residualNonCubicFp} is not).
	 */
	public static boolean passesRandomMatmulSpotCheckFp(NonCubicBilinearAlgorithm alg, int p, int samples) {
		int n = alg.n, m = alg.m, pp = alg.p, r = alg.r;
		int dimU = n * m, dimV = m * pp, dimW = n * pp;
		int[][] uM = new int[dimU][r], vM = new int[dimV][r], wM = new int[dimW][r];
		if (!reduceFactorModP(alg.denseU(), uM, p) || !reduceFactorModP(alg.denseV(), vM, p)
				|| !reduceFactorModP(alg.denseW(), wM, p)) {
			return false;  // not representable mod p
		}
		java.util.Random rng = new java.util.Random(0xC0DE);
		for (int trial = 0; trial < samples; trial++) {
			int[] aFlat = new int[dimU], bFlat = new int[dimV];
			for (int x = 0; x < dimU; x++) aFlat[x] = rng.nextInt(p);
			for (int x = 0; x < dimV; x++) bFlat[x] = rng.nextInt(p);
			long[] alpha = new long[r], beta = new long[r];
			for (int k = 0; k < r; k++) {
				long sa = 0, sb = 0;
				for (int a = 0; a < dimU; a++) sa += (long) uM[a][k] * aFlat[a];
				for (int b = 0; b < dimV; b++) sb += (long) vM[b][k] * bFlat[b];
				alpha[k] = ((sa % p) + p) % p;
				beta[k] = ((sb % p) + p) % p;
			}
			for (int i = 0; i < n; i++) {
				for (int l = 0; l < pp; l++) {
					int c = i * pp + l;
					long got = 0;
					for (int k = 0; k < r; k++) got += (long) wM[c][k] * ((alpha[k] * beta[k]) % p);
					long want = 0;
					for (int j = 0; j < m; j++) want += (long) aFlat[i * m + j] * bFlat[j * pp + l];
					if ((((got - want) % p) + p) % p != 0) return false;
				}
			}
		}
		return true;
	}

	/** Convenience: 5 samples over GF(p). */
	public static boolean passesRandomMatmulSpotCheckFp(NonCubicBilinearAlgorithm alg, int p) {
		return passesRandomMatmulSpotCheckFp(alg, p, 5);
	}

	/**
	 * Single bilinear-coefficient discrepancy: in {@code C[i,l]} the
	 * algorithm puts coefficient {@code actual} on the term
	 * {@code A[αi,αj]·B[βj,βl]}, but the target matmul polynomial expects
	 * {@code expected} (= 1 iff αi=i ∧ αj=βj ∧ βl=l, else 0).
	 */
	public record SymbolicDiff(int outI, int outL, int aI, int aJ, int bJ, int bL,
			double actual, double expected) {
		public double residual() { return actual - expected; }
		@Override public String toString() {
			return String.format("C[%d,%d] : A[%d,%d]·B[%d,%d] -- actual=%g expected=%g (Δ=%g)",
					outI, outL, aI, aJ, bJ, bL, actual, expected, residual());
		}
	}

	/**
	 * Symbolic verification of a non-cubic bilinear algorithm: treat each
	 * {@code A[i,j]} and {@code B[j,l]} as a formal indeterminate, expand
	 * the algorithm's polynomial for every output {@code C[i,l]}, and
	 * compare term-by-term to the target {@code Σ_j A[i,j]·B[j,l]}.
	 *
	 * <p>Returns the LIST of discrepancies (empty list iff the algorithm
	 * is exact). Each entry pinpoints the exact bilinear term whose
	 * coefficient is off, which is what's needed to debug a broken
	 * recombined scheme — unlike the Frobenius residual which only tells
	 * you "something's wrong, somewhere".</p>
	 *
	 * <p>Time complexity: O((nm)·(mp)·(np)·r) — same as
	 * {@link #residualNonCubic} but emits typed diagnostics. For very
	 * large schemes prefer the fast {@link #passesRandomMatmulSpotCheck}
	 * for a yes/no answer, and use this only on smaller cases or to
	 * inspect the FIRST few discrepancies of a known-broken scheme.</p>
	 *
	 * @param alg          the scheme
	 * @param maxReport    cap the result list at this many discrepancies
	 *                     (set Integer.MAX_VALUE for full enumeration)
	 * @param tolerance    coefficient difference treated as zero
	 */
	public static java.util.List<SymbolicDiff> symbolicDiscrepancies(NonCubicBilinearAlgorithm alg,
			int maxReport, double tolerance) {
		java.util.List<SymbolicDiff> out = new java.util.ArrayList<>();
		forEachTerm(alg, (a, b, c, expected, actual) -> {
			if (Math.abs(actual - expected) > tolerance) {
				// Unflatten row-major indices back to cell coordinates.
				int aI = a / alg.m, aJ = a % alg.m;
				int bJ = b / alg.p, bL = b % alg.p;
				int i  = c / alg.p, l  = c % alg.p;
				out.add(new SymbolicDiff(i, l, aI, aJ, bJ, bL, actual, expected));
			}
			return out.size() < maxReport;
		});
		return out;
	}

	/** Convenience: first 20 discrepancies at default tolerance. */
	public static java.util.List<SymbolicDiff> symbolicDiscrepancies(NonCubicBilinearAlgorithm alg) {
		return symbolicDiscrepancies(alg, 20, 1e-10);
	}

	/**
	 * Verify a non-bilinear / quadratic algorithm for matmul over a
	 * COMMUTATIVE ring (per DIS09 §1.2). Each product
	 * {@code γ_k = α_k · β_k} where both {@code α_k} and {@code β_k}
	 * may carry coefficients on both A and B; correctness requires that
	 * scalar entries of A and B commute (so {@code A[i,j]·B[j',l'] =
	 * B[j',l']·A[i,j]} and similarly within A·A and B·B).
	 *
	 * <p>Three classes of polynomial-coefficient identities are checked:</p>
	 * <ol>
	 *   <li>A·B cross-terms (symmetrised over commutativity) must equal
	 *       the matmul-tensor coefficient {@code δ(α=i)·δ(β=γ)·δ(δ=l)}.</li>
	 *   <li>A·A cross-terms (symmetrised) must be zero.</li>
	 *   <li>B·B cross-terms (symmetrised) must be zero.</li>
	 * </ol>
	 *
	 * <p>Returns the Frobenius residual over all three constraint sets;
	 * residual = 0 iff the algorithm is exact (over the commutative ring).</p>
	 */
	public static double residualNonBilinear(NonBilinearAlgorithm alg) {
		int n = alg.n, m = alg.m, p = alg.p, r = alg.r;
		double sumSq = 0.0;

		// (1) A·B cross-terms (commutative-symmetrised).
		// For each (i, l, α, β, γ, δ): coefficient of A[α,β]·B[γ,δ] in C[i,l]
		// (from algorithm) = Σ_k W[il][k] · (Ua[αβ][k]·Vb[γδ][k] + Va[αβ][k]·Ub[γδ][k]).
		// Target = 1 iff (α=i, β=γ, δ=l), else 0.
		for (int i = 0; i < n; i++) {
			for (int l = 0; l < p; l++) {
				int ilIdx = i * p + l;
				for (int alpha = 0; alpha < n; alpha++) {
					for (int beta = 0; beta < m; beta++) {
						int abIdx = alpha * m + beta;
						for (int gamma = 0; gamma < m; gamma++) {
							for (int delta = 0; delta < p; delta++) {
								int gdIdx = gamma * p + delta;
								double approx = 0.0;
								for (int k = 0; k < r; k++) {
									approx += alg.W[ilIdx][k] *
											(alg.Ua[abIdx][k] * alg.Vb[gdIdx][k]
													+ alg.Va[abIdx][k] * alg.Ub[gdIdx][k]);
								}
								double target = (alpha == i && beta == gamma && delta == l) ? 1.0 : 0.0;
								double diff = target - approx;
								sumSq += diff * diff;
							}
						}
					}
				}
			}
		}

		// (2) A·A cross-terms must be zero (commutative-symmetrised over (α,β) ↔ (α',β')).
		// Σ_k W[il][k] · (Ua[αβ][k]·Va[α'β'][k] + Ua[α'β'][k]·Va[αβ][k]) = 0
		// We only check unordered pairs to avoid double-counting; the diagonal case
		// (α,β)=(α',β') is the unsymmetrised condition Σ_k W·Ua·Va = 0.
		for (int i = 0; i < n; i++) {
			for (int l = 0; l < p; l++) {
				int ilIdx = i * p + l;
				for (int ab = 0; ab < n * m; ab++) {
					for (int ab2 = ab; ab2 < n * m; ab2++) {
						double approx = 0.0;
						for (int k = 0; k < r; k++) {
							if (ab == ab2) {
								approx += alg.W[ilIdx][k] * alg.Ua[ab][k] * alg.Va[ab][k];
							} else {
								approx += alg.W[ilIdx][k] *
										(alg.Ua[ab][k] * alg.Va[ab2][k] + alg.Ua[ab2][k] * alg.Va[ab][k]);
							}
						}
						sumSq += approx * approx;
					}
				}
			}
		}

		// (3) B·B cross-terms must be zero.
		for (int i = 0; i < n; i++) {
			for (int l = 0; l < p; l++) {
				int ilIdx = i * p + l;
				for (int gd = 0; gd < m * p; gd++) {
					for (int gd2 = gd; gd2 < m * p; gd2++) {
						double approx = 0.0;
						for (int k = 0; k < r; k++) {
							if (gd == gd2) {
								approx += alg.W[ilIdx][k] * alg.Ub[gd][k] * alg.Vb[gd][k];
							} else {
								approx += alg.W[ilIdx][k] *
										(alg.Ub[gd][k] * alg.Vb[gd2][k] + alg.Ub[gd2][k] * alg.Vb[gd][k]);
							}
						}
						sumSq += approx * approx;
					}
				}
			}
		}

		return Math.sqrt(sumSq);
	}

	public static boolean isExactNonBilinear(NonBilinearAlgorithm alg) {
		return residualNonBilinear(alg) < 1e-10;
	}

	/**
	 * Fast random-input spot-check for {@link NonBilinearAlgorithm}:
	 * evaluates the algorithm on {@code samples} random {@code A,B}
	 * matrices (entries uniform in {@code [-magnitude, magnitude]}) and
	 * compares against the naive {@code A·B}. Returns {@code true} iff
	 * the Frobenius error is below {@code eps} for every sample.
	 *
	 * <p>O({@code samples · r · (n·m + m·p + n·p)}) — millions of times
	 * faster than {@link #residualNonBilinear} for large {@code n}, but
	 * a randomised correctness witness rather than a full algebraic proof.
	 * For schemes with integer or low-denominator rational coefficients,
	 * a single random sample is usually decisive.</p>
	 */
	public static boolean passesRandomMatmulSpotCheckNB(NonBilinearAlgorithm alg,
			int samples, double magnitude, double eps, long seed) {
		int n = alg.n, m = alg.m, p = alg.p, r = alg.r;
		java.util.Random rng = new java.util.Random(seed);
		for (int s = 0; s < samples; s++) {
			double[] aFlat = new double[n * m];
			double[] bFlat = new double[m * p];
			for (int i = 0; i < aFlat.length; i++) aFlat[i] = (rng.nextDouble() * 2 - 1) * magnitude;
			for (int i = 0; i < bFlat.length; i++) bFlat[i] = (rng.nextDouble() * 2 - 1) * magnitude;

			// Naive A·B
			double[] cNaive = new double[n * p];
			for (int i = 0; i < n; i++)
				for (int l = 0; l < p; l++) {
					double sum = 0;
					for (int j = 0; j < m; j++) sum += aFlat[i * m + j] * bFlat[j * p + l];
					cNaive[i * p + l] = sum;
				}

			// Algorithm output
			double[] alpha = new double[r];
			double[] beta = new double[r];
			for (int k = 0; k < r; k++) {
				double a = 0, b = 0;
				for (int ab = 0; ab < n * m; ab++) {
					a += alg.Ua[ab][k] * aFlat[ab];
					b += alg.Va[ab][k] * aFlat[ab];
				}
				for (int gd = 0; gd < m * p; gd++) {
					a += alg.Ub[gd][k] * bFlat[gd];
					b += alg.Vb[gd][k] * bFlat[gd];
				}
				alpha[k] = a;
				beta[k] = b;
			}
			double[] cAlgo = new double[n * p];
			for (int il = 0; il < n * p; il++) {
				double sum = 0;
				for (int k = 0; k < r; k++) sum += alg.W[il][k] * alpha[k] * beta[k];
				cAlgo[il] = sum;
			}

			double err = 0;
			for (int i = 0; i < n * p; i++) {
				double d = cAlgo[i] - cNaive[i];
				err += d * d;
			}
			err = Math.sqrt(err);
			double scale = 0;
			for (double v : cNaive) scale += v * v;
			scale = Math.sqrt(scale) + 1e-12;
			if (err / scale > eps) return false;
		}
		return true;
	}

	/** Default 5-sample, magnitude=1, eps=1e-9, seed=0 spot-check. */
	public static boolean passesRandomMatmulSpotCheckNB(NonBilinearAlgorithm alg) {
		return passesRandomMatmulSpotCheckNB(alg, 5, 1.0, 1e-9, 0L);
	}

	/**
	 * Complex-arithmetic Frobenius residual against the (real-valued) matmul
	 * tensor. The decomposition is exact when both the real-residual and the
	 * imaginary part of the sum vanish: {@code Σ_k U_k · V_k · W_k} must equal
	 * the matmul tensor entries (which are integers, hence have zero imaginary
	 * part).
	 *
	 * Used for AlphaEvolve-style algorithms with Gaussian-rational coefficients.
	 */
	public static double residualComplex(ComplexNonCubicBilinearAlgorithm alg) {
		int[][][] target = intMatmulTensor(alg.n, alg.m, alg.p);
		int dimU = alg.n * alg.m, dimV = alg.m * alg.p, dimW = alg.n * alg.p;
		double sumSq = 0.0;
		for (int a = 0; a < dimU; a++) {
			for (int b = 0; b < dimV; b++) {
				for (int c = 0; c < dimW; c++) {
					double approxRe = 0.0, approxIm = 0.0;
					for (int k = 0; k < alg.r; k++) {
						double ur = alg.uRe[a][k], ui = alg.uIm[a][k];
						double vr = alg.vRe[b][k], vi = alg.vIm[b][k];
						double wr = alg.wRe[c][k], wi = alg.wIm[c][k];
						// (u_re + i u_im) * (v_re + i v_im)
						double uvRe = ur * vr - ui * vi;
						double uvIm = ur * vi + ui * vr;
						// * (w_re + i w_im)
						approxRe += uvRe * wr - uvIm * wi;
						approxIm += uvRe * wi + uvIm * wr;
					}
					double diffRe = target[a][b][c] - approxRe;
					sumSq += diffRe * diffRe + approxIm * approxIm;
				}
			}
		}
		return Math.sqrt(sumSq);
	}

	public static boolean isExactComplex(ComplexNonCubicBilinearAlgorithm alg) {
		return residualComplex(alg) < 1e-9;
	}

	/**
	 * F₂-aware residual: the trilinear sum is computed in integer arithmetic
	 * and then reduced mod 2 before comparing against the matmul tensor. Many
	 * F₂-specific schemes (e.g. AlphaTensor) use {0, 1} coefficients where
	 * real-arithmetic sums can hit 2 or 3 — the F₂ verifier identifies those
	 * as 0 or 1 respectively.
	 */
	public static int residualNonCubicF2(NonCubicBilinearAlgorithm alg) {
		double[][] srcU = alg.denseU(), srcV = alg.denseV(), srcW = alg.denseW();
		int[][][] target = intMatmulTensor(alg.n, alg.m, alg.p);
		int dimU = alg.n * alg.m, dimV = alg.m * alg.p, dimW = alg.n * alg.p;
		int wrong = 0;
		for (int a = 0; a < dimU; a++) {
			for (int b = 0; b < dimV; b++) {
				for (int c = 0; c < dimW; c++) {
					int sum = 0;
					for (int k = 0; k < alg.r; k++) {
						int u = (int) Math.round(srcU[a][k]);
						int v = (int) Math.round(srcV[b][k]);
						int w = (int) Math.round(srcW[c][k]);
						// GF(2) bilinear product, masked to bit 0: the parity of a product
						// is the AND of the parities (incl. two's-complement negatives),
						// but for coefficients outside {0,1} (e.g. three -1s: -1&-1&-1 =
						// -1) the upper bits of u&v&w are garbage and would pollute the
						// XOR accumulator — falsely rejecting mod-2-exact Z schemes.
						sum ^= (u & v & w) & 1;
					}
					if (sum != target[a][b][c]) wrong++;
				}
			}
		}
		return wrong;
	}

	public static boolean isExactNonCubicF2(NonCubicBilinearAlgorithm alg) {
		return residualNonCubicF2(alg) == 0;
	}

	/**
	 * F₃-aware residual: trilinear sum computed in integer arithmetic and
	 * then reduced mod 3 before comparing against the matmul tensor.
	 *
	 * <p>Caller assumes integer (or 1/2-style) coefficients; this method
	 * rounds each entry to a long and reduces mod 3, so non-integer
	 * coefficients with denominators coprime-to-3 (e.g. 1/2 = 2 in F₃)
	 * silently round to the nearest integer — pre-filter via
	 * {@code firstNonIntegerEntry} if that distinction matters.</p>
	 */
	public static int residualNonCubicF3(NonCubicBilinearAlgorithm alg) {
		double[][] srcU = alg.denseU(), srcV = alg.denseV(), srcW = alg.denseW();
		int[][][] target = intMatmulTensor(alg.n, alg.m, alg.p);
		int dimU = alg.n * alg.m, dimV = alg.m * alg.p, dimW = alg.n * alg.p;
		int wrong = 0;
		for (int a = 0; a < dimU; a++) {
			for (int b = 0; b < dimV; b++) {
				for (int c = 0; c < dimW; c++) {
					int sum = 0;
					for (int k = 0; k < alg.r; k++) {
						int u = mod3((long) Math.round(srcU[a][k]));
						int v = mod3((long) Math.round(srcV[b][k]));
						int w = mod3((long) Math.round(srcW[c][k]));
						sum = (sum + u * v * w) % 3;
					}
					if (sum != mod3(target[a][b][c])) wrong++;
				}
			}
		}
		return wrong;
	}

	public static boolean isExactNonCubicF3(NonCubicBilinearAlgorithm alg) {
		return residualNonCubicF3(alg) == 0;
	}

	private static int mod3(long x) {
		int r = (int) (x % 3);
		return r < 0 ? r + 3 : r;
	}

	// ── rational reduction mod a prime ─────────────────────────────────────────

	/**
	 * Result of {@link #residualNonCubicFp}: {@code representable} is false when
	 * some coefficient's denominator is divisible by {@code p} (the scheme cannot
	 * be reduced mod p AS WRITTEN — an equivalent F_p scheme may still exist, but
	 * that would be a different scheme); {@code wrong} is the residual count when
	 * representable.
	 */
	public record FpReduction(boolean representable, int wrong) {
		public boolean exact() { return representable && wrong == 0; }
	}

	/**
	 * Reduce a RATIONAL-coefficient scheme mod prime {@code p} and verify it.
	 * Each coefficient {@code num/den} (recovered exactly from the double via
	 * continued fractions) maps to {@code num · den⁻¹ mod p}, which exists iff
	 * {@code gcd(den, p) = 1}. So {@code 1/2} is representable mod 3
	 * ({@code 2⁻¹ ≡ 2}) but NOT mod 2, and {@code 1/3} is representable mod 2
	 * ({@code 3 ≡ 1}) but NOT mod 3 — a Q-exact scheme is F_p-valid exactly when
	 * every denominator is coprime to p. Integer coefficients (den = 1) always
	 * reduce, so this subsumes the integer case.
	 */
	public static FpReduction residualNonCubicFp(NonCubicBilinearAlgorithm alg, int p) {
		double[][] srcU = alg.denseU(), srcV = alg.denseV(), srcW = alg.denseW();
		int dimU = alg.n * alg.m, dimV = alg.m * alg.p, dimW = alg.n * alg.p;
		int[][] uM = new int[dimU][alg.r], vM = new int[dimV][alg.r], wM = new int[dimW][alg.r];
		if (!reduceFactorModP(srcU, uM, p) || !reduceFactorModP(srcV, vM, p)
				|| !reduceFactorModP(srcW, wM, p)) {
			return new FpReduction(false, -1);
		}
		int wrong = 0;
		for (int a = 0; a < dimU; a++) {
			for (int b = 0; b < dimV; b++) {
				for (int c = 0; c < dimW; c++) {
					long sum = 0;
					for (int k = 0; k < alg.r; k++) {
						sum += (long) uM[a][k] * vM[b][k] * wM[c][k];
					}
					int got = (int) (((sum % p) + p) % p);
					int want = matmulTensorEntry(a, b, c, alg.n, alg.m, alg.p) % p;
					if (got != want) wrong++;
				}
			}
		}
		return new FpReduction(true, wrong);
	}

	/** True iff the scheme reduces mod p (all denominators coprime to p) AND verifies. */
	public static boolean isExactNonCubicFp(NonCubicBilinearAlgorithm alg, int p) {
		return residualNonCubicFp(alg, p).exact();
	}

	/** @return false iff some entry's denominator is divisible by {@code p}. */
	private static boolean reduceFactorModP(double[][] src, int[][] dst, int p) {
		for (int i = 0; i < src.length; i++) {
			for (int k = 0; k < src[i].length; k++) {
				int m = coefModP(src[i][k], p);
				if (m < 0) return false;
				dst[i][k] = m;
			}
		}
		return true;
	}

	/**
	 * Map a coefficient to its residue mod p, or -1 if not representable.
	 * Integers take the fast path; non-integers are rationalised by
	 * continued-fraction convergents (same approach as SchemeIO's writer).
	 */
	private static int coefModP(double v, int p) {
		long rounded = Math.round(v);
		if (Math.abs(v - rounded) < 1e-9) {
			return (int) (((rounded % p) + p) % p);
		}
		long[] frac = rationalizeToFraction(v);
		if (frac == null) return -1; // irrational / huge denominator: not reducible
		long num = frac[0], den = frac[1];
		if (den % p == 0) return -1; // denominator kills the reduction (e.g. 1/2 mod 2)
		int denInv = modInverse((int) (((den % p) + p) % p), p);
		long numMod = ((num % p) + p) % p;
		return (int) ((numMod * denInv) % p);
	}

	/** Continued-fraction recovery of the exact {@code [num, den]} a double encodes. */
	private static long[] rationalizeToFraction(double v) {
		final long DEN_CAP = 1_000_000L;
		final double TOL = 1e-12;
		boolean neg = v < 0;
		double x = Math.abs(v);
		long hPrev = 0, h = 1, kPrev = 1, k = 0;
		double frac = x;
		for (int iter = 0; iter < 64; iter++) {
			long a = (long) Math.floor(frac);
			long hNext = a * h + hPrev;
			long kNext = a * k + kPrev;
			if (kNext > DEN_CAP || kNext <= 0) return null;
			hPrev = h; h = hNext;
			kPrev = k; k = kNext;
			if (k > 0 && Math.abs((double) h / (double) k - x) <= TOL * Math.max(1.0, x)) {
				return new long[] { neg ? -h : h, k };
			}
			double rem = frac - a;
			if (rem < 1e-15) return null; // converged on an integer — handled upstream
			frac = 1.0 / rem;
		}
		return null;
	}

	/** Modular inverse for prime modulus via Fermat ({@code a^(p-2) mod p}). */
	private static int modInverse(int a, int p) {
		long result = 1, base = a, exp = p - 2L;
		while (exp > 0) {
			if ((exp & 1) == 1) result = result * base % p;
			base = base * base % p;
			exp >>= 1;
		}
		return (int) result;
	}

	/**
	 * Sampled residual: picks {@code samples} random {@code (a, b, c)} positions
	 * from the matmul tensor, computes the trilinear sum at each, and returns
	 * the count of mismatches. For sizes where the full `O(N⁶·r)` verifier is
	 * intractable (max-dim ≥ 12), this gives statistically-strong confidence at
	 * `O(samples · r)` cost.
	 *
	 * <p>Caller convention: a return value of {@code 0} after {@code 10_000+}
	 * samples on a Kronecker-constructed algorithm is effectively
	 * conclusive — the construction is provably correct, so any per-position
	 * mismatch would be detected with high probability.</p>
	 */
	public static int residualSampled(NonCubicBilinearAlgorithm alg, int samples, long seed) {
		// Row-major (CSR) snapshots: ~nnz memory instead of densifying the full
		// rows·r grid (413 MB → single-digit MB at ⟨32,32,32⟩). The trilinear
		// sum at one (a,b,c) is a three-row merge over shared product-columns.
		RowMajorFactor U = RowMajorFactor.of(alg.u());
		RowMajorFactor V = RowMajorFactor.of(alg.v());
		RowMajorFactor W = RowMajorFactor.of(alg.w());
		java.util.Random rnd = new java.util.Random(seed);
		int dimU = alg.n * alg.m, dimV = alg.m * alg.p, dimW = alg.n * alg.p;
		int wrong = 0;
		for (int s = 0; s < samples; s++) {
			int a = rnd.nextInt(dimU);
			int b = rnd.nextInt(dimV);
			int c = rnd.nextInt(dimW);
			double approx = RowMajorFactor.triProduct(U, a, V, b, W, c);
			if (Math.abs(matmulTensorEntry(a, b, c, alg.n, alg.m, alg.p) - approx) > 1e-9) wrong++;
		}
		return wrong;
	}

	/** Sampled complex residual — same idea, complex arithmetic. */
	public static int residualSampledComplex(ComplexNonCubicBilinearAlgorithm alg,
			int samples, long seed) {
		java.util.Random rnd = new java.util.Random(seed);
		int dimU = alg.n * alg.m, dimV = alg.m * alg.p, dimW = alg.n * alg.p;
		int wrong = 0;
		for (int s = 0; s < samples; s++) {
			int a = rnd.nextInt(dimU);
			int b = rnd.nextInt(dimV);
			int c = rnd.nextInt(dimW);
			double approxRe = 0.0, approxIm = 0.0;
			for (int k = 0; k < alg.r; k++) {
				double ur = alg.uRe[a][k], ui = alg.uIm[a][k];
				double vr = alg.vRe[b][k], vi = alg.vIm[b][k];
				double wr = alg.wRe[c][k], wi = alg.wIm[c][k];
				double uvRe = ur * vr - ui * vi;
				double uvIm = ur * vi + ui * vr;
				approxRe += uvRe * wr - uvIm * wi;
				approxIm += uvRe * wi + uvIm * wr;
			}
			double diffRe = matmulTensorEntry(a, b, c, alg.n, alg.m, alg.p) - approxRe;
			if (diffRe * diffRe + approxIm * approxIm > 1e-9) wrong++;
		}
		return wrong;
	}

	/** Sampled F₂ residual — same idea, XOR semantics. */
	public static int residualSampledF2(NonCubicBilinearAlgorithm alg, int samples, long seed) {
		// Row-major (CSR) snapshots — see residualSampled. ~nnz memory, no densify.
		RowMajorFactor U = RowMajorFactor.of(alg.u());
		RowMajorFactor V = RowMajorFactor.of(alg.v());
		RowMajorFactor W = RowMajorFactor.of(alg.w());
		java.util.Random rnd = new java.util.Random(seed);
		int dimU = alg.n * alg.m, dimV = alg.m * alg.p, dimW = alg.n * alg.p;
		int wrong = 0;
		for (int s = 0; s < samples; s++) {
			int a = rnd.nextInt(dimU);
			int b = rnd.nextInt(dimV);
			int c = rnd.nextInt(dimW);
			int sum = RowMajorFactor.triAndXor(U, a, V, b, W, c);
			if (sum != matmulTensorEntry(a, b, c, alg.n, alg.m, alg.p)) wrong++;
		}
		return wrong;
	}

	/**
	 * Count linear-combination additions in a real-valued bilinear algorithm.
	 *
	 * <p>Per the standard convention (Strassen 1969, Smirnov 2013):</p>
	 * <ul>
	 *   <li>Each rank-k multiplication consumes one linear combination of A
	 *       entries (the k-th column of U) and one of B (k-th column of V).
	 *       A column with {@code c} non-zero coefficients takes {@code c - 1}
	 *       additions to assemble — so input-side adds total
	 *       {@code nonzeros(U) - r + nonzeros(V) - r}.</li>
	 *   <li>Each output {@code C[i, j]} is then a linear combination of the
	 *       rank multiplications via the {@code (i·p + j)}-th row of W. A row
	 *       with {@code c} nonzeros costs {@code c - 1} additions — so
	 *       output-side adds total {@code nonzeros(W) - dimW}.</li>
	 * </ul>
	 *
	 * <p>Result: {@code nz(U) + nz(V) + nz(W) - 2·r - n·p}. Coefficients of
	 * magnitude other than 1 (e.g. AlphaEvolve's `0.5` or `2`) are NOT counted
	 * as extra ops here — that's a separate scaling-cost question.</p>
	 *
	 * <p>Strassen ⟨2,2,2⟩ r=7 example: nz(U) = nz(V) = nz(W) = 12, dimW = 4 →
	 * {@code 36 - 14 - 4 = 18 additions}, matching Strassen's published count.</p>
	 */
	public static int additionCount(NonCubicBilinearAlgorithm alg) {
		int dimW = alg.n * alg.p;
		return alg.u().nonZeros() + alg.v().nonZeros() + alg.w().nonZeros()
				- 2 * alg.r - dimW;
	}

	public static int additionCount(BilinearAlgorithm alg) {
		return additionCount(NonCubicBilinearAlgorithm.fromCubic(alg));
	}

	/**
	 * Count linear-combination additions in a non-bilinear algorithm
	 * (Waksman, Rosowski Algorithm 1, …). Each rank-1 product is
	 * {@code (Ua·vec(A) + Ub·vec(B)) · (Va·vec(A) + Vb·vec(B))} — so per
	 * product, the first factor's combo costs {@code nnz(Ua col) + nnz(Ub col) - 1}
	 * additions (or 0 if empty), and likewise the second factor. Per output
	 * position, {@code nnz(W row) - 1} additions to assemble (or 0 if empty).
	 *
	 * <p>Half-integer coefficients (Waksman's W carries ±0.5) are counted as
	 * non-zero contributions only — the scaling op is NOT counted as an
	 * addition. That's a separate cost question.</p>
	 */
	public static int additionCount(NonBilinearAlgorithm alg) {
		int dimA = alg.n * alg.m;
		int dimB = alg.m * alg.p;
		int dimC = alg.n * alg.p;
		int adds = 0;
		for (int k = 0; k < alg.r; k++) {
			int firstFactor = nonzerosInColumn(alg.Ua, k, dimA)
					+ nonzerosInColumn(alg.Ub, k, dimB);
			int secondFactor = nonzerosInColumn(alg.Va, k, dimA)
					+ nonzerosInColumn(alg.Vb, k, dimB);
			if (firstFactor > 1) adds += firstFactor - 1;
			if (secondFactor > 1) adds += secondFactor - 1;
		}
		for (int i = 0; i < dimC; i++) {
			int t = 0;
			for (int k = 0; k < alg.r; k++) if (alg.W[i][k] != 0) t++;
			if (t > 1) adds += t - 1;
		}
		return adds;
	}

	private static int nonzerosInColumn(double[][] matrix, int col, int rows) {
		int c = 0;
		for (int i = 0; i < rows; i++) if (matrix[i][col] != 0) c++;
		return c;
	}

	private static int countNonZero(double[][] matrix) {
		int c = 0;
		for (double[] row : matrix) {
			for (double x : row) {
				if (x != 0.0) c++;
			}
		}
		return c;
	}

	/**
	 * Complex-version of {@link #additionCount(NonCubicBilinearAlgorithm)}.
	 * A complex coefficient is "non-zero" iff either its real or imaginary
	 * part is non-zero.
	 */
	public static int additionCount(ComplexNonCubicBilinearAlgorithm alg) {
		int dimW = alg.n * alg.p;
		int nzU = countNonZeroComplex(alg.uRe, alg.uIm);
		int nzV = countNonZeroComplex(alg.vRe, alg.vIm);
		int nzW = countNonZeroComplex(alg.wRe, alg.wIm);
		return nzU + nzV + nzW - 2 * alg.r - dimW;
	}

	private static int countNonZeroComplex(double[][] re, double[][] im) {
		int c = 0;
		for (int i = 0; i < re.length; i++) {
			for (int k = 0; k < re[i].length; k++) {
				if (re[i][k] != 0.0 || im[i][k] != 0.0) c++;
			}
		}
		return c;
	}

	/**
	 * Returns a new {@link BilinearAlgorithm} whose W factor has its row indices
	 * transposed: row {@code i·n + j} ↔ row {@code j·n + i}. Converts a rank-r
	 * decomposition of {@link #matmulTensor} to a rank-r decomposition of
	 * {@link #trilinTensor}, and vice versa (the operation is self-inverse).
	 *
	 * Why: {@code trace(A·B·C) = Σ_{i,j} (A·B)[i,j] · C[j,i]}. A coefficient
	 * vector w that contracts as {@code C_flat[i·n + j] · w[c]} in the
	 * matmul-tensor convention contracts as {@code C_flat[j·n + i] · w[c]} in
	 * the trilinear-form convention — i.e. with row indices swapped.
	 */
	public static BilinearAlgorithm transposeW(BilinearAlgorithm alg) {
		int n = alg.n;
		int n2 = n * n;
		double[][] newW = new double[n2][alg.r];
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				int dst = i * n + j;
				int src = j * n + i;
				for (int k = 0; k < alg.r; k++) {
					newW[dst][k] = alg.W[src][k];
				}
			}
		}
		return new BilinearAlgorithm(n, alg.U, alg.V, newW);
	}
}
