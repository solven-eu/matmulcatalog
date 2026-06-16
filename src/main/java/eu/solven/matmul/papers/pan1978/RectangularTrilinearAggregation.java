package eu.solven.matmul.papers.pan1978;

import java.util.ArrayList;
import java.util.List;

import eu.solven.matmul.NonCubicBilinearAlgorithm;

/**
 * Rectangular trilinear aggregation (bucket (b), see
 * {@code references/RECTANGULAR_TA.md}). FMM fuses the two off-diagonal
 * cross-blocks of a symmetric peel — which are CYCLIC ROTATIONS of one shape —
 * via Pan TA (e.g. {@code TA(⟨26,3,26⟩,⟨26,26,3⟩)=2860} in ⟨26,29,29⟩=11693).
 *
 * <p>This class holds the SAFE, math-free part: detecting the disjoint
 * cyclic-rotation leaf-pairs a recombination produces. The fused RANK and the
 * factor-matrix CONSTRUCTION are an open derivation (rectangularising Islam's
 * SQUARE {@link PanTrilinearAggregationBuilder}); {@link #fusedRank} /
 * {@link #build} deliberately throw until that math is worked out AND
 * exact-verified — shipping an unverified TA construction is forbidden.</p>
 */
public final class RectangularTrilinearAggregation {

	private RectangularTrilinearAggregation() {}

	/** {@code true} iff {@code b} is a NON-trivial cyclic rotation of {@code a}
	 *  — i.e. {@code b ∈ {rot¹(a), rot²(a)}} and {@code a} is not already
	 *  rotation-symmetric (all three dims equal would make every rotation equal,
	 *  a trivial / non-fusable case). Order matters: ⟨a,b,c⟩→⟨b,c,a⟩→⟨c,a,b⟩. */
	public static boolean isCyclicRotation(int[] a, int[] b) {
		if (a.length != 3 || b.length != 3) {
			return false;
		}
		if (a[0] == a[1] && a[1] == a[2]) {
			return false;  // ⟨n,n,n⟩: rotations coincide — not a distinct pair
		}
		boolean rot1 = b[0] == a[1] && b[1] == a[2] && b[2] == a[0];
		boolean rot2 = b[0] == a[2] && b[1] == a[0] && b[2] == a[1];
		boolean identity = b[0] == a[0] && b[1] == a[1] && b[2] == a[2];
		return (rot1 || rot2) && !identity;
	}

	/**
	 * Index pairs {@code [i,j]} (i&lt;j) among the per-product leaf shapes that are
	 * cyclic rotations of each other — the candidates a recombination could
	 * TA-fuse. Each leaf comes from a distinct base product, so a returned pair is
	 * automatically a DISJOINT pair (the precondition Pan TA needs). Pure
	 * combinatorics — no rank/cost claim (that needs the open derivation).
	 */
	public static List<int[]> cyclicRotationPairs(int[][] leafShapes) {
		List<int[]> out = new ArrayList<>();
		for (int i = 0; i < leafShapes.length; i++) {
			for (int j = i + 1; j < leafShapes.length; j++) {
				if (isCyclicRotation(leafShapes[i], leafShapes[j])) {
					out.add(new int[] { i, j });
				}
			}
		}
		return out;
	}

	/** A signed rank-one trilinear term {@code sign · u⊗v⊗w} (the sign is folded
	 *  into {@code u}). Contribution to the tensor entry {@code (i,j,k)} is
	 *  {@code u[i]·v[j]·w[k]}. */
	public record Term(int[] u, int[] v, int[] w) {}

	/**
	 * Pan 1978 aggregation identity (3) — the ATOM of trilinear aggregation.
	 * Rewrites the sum of two rank-one terms {@code a⊗b⊗c + a'⊗b'⊗c'} as ONE
	 * aggregated product plus THREE correction products (4 terms total):
	 * <pre>
	 *   a⊗b⊗c + a'⊗b'⊗c'
	 *     =  (a+a')⊗(b+b')⊗(c+c')      // aggregated
	 *      − a'    ⊗(b+b')⊗ c          // correction 1
	 *      − a     ⊗ b'   ⊗(c+c')      // correction 2
	 *      −(a+a') ⊗ b    ⊗ c'         // correction 3
	 * </pre>
	 * On its own this is 4 terms for 2 products (a loss); the win comes only when
	 * the corrections UNITE across a kin family (Pan identity (4)) — that step is
	 * the open part. This atom is exact (verified by {@code TestRectangularTrilinearAggregation}).
	 *
	 * @param a,b,c   factor vectors of the first term ({@code u/v/w} spaces)
	 * @param a2,b2,c2 factor vectors of the second (kin) term, same spaces
	 * @return the 4 terms whose tensor sum equals {@code a⊗b⊗c + a2⊗b2⊗c2}
	 */
	public static List<Term> aggregatePair(int[] a, int[] b, int[] c,
			int[] a2, int[] b2, int[] c2) {
		List<Term> out = new ArrayList<>(4);
		out.add(new Term(add(a, a2), add(b, b2), add(c, c2)));   // (a+a')(b+b')(c+c')
		out.add(new Term(neg(a2), add(b, b2), c.clone()));       // − a'(b+b')c
		out.add(new Term(neg(a), b2.clone(), add(c, c2)));       // − a b'(c+c')
		out.add(new Term(neg(add(a, a2)), b.clone(), c2.clone()));// −(a+a')b c'
		return out;
	}

	private static int[] add(int[] x, int[] y) {
		int[] r = x.clone();
		for (int i = 0; i < r.length; i++) r[i] += y[i];
		return r;
	}

	private static int[] neg(int[] x) {
		int[] r = new int[x.length];
		for (int i = 0; i < r.length; i++) r[i] = -x[i];
		return r;
	}

	/**
	 * Disjoint-sum aggregation: given two rank-one decompositions {@code p1}
	 * (over its own u/v/w spaces of sizes {@code du1,dv1,dw1}) and {@code p2}
	 * ({@code du2,dv2,dw2}) — the two cross-blocks of a symmetric peel, which use
	 * DISJOINT a/b/c entries — embed both in the COMBINED space (p2 offset past
	 * p1's blocks), pair term-by-term and aggregate each pair via {@link
	 * #aggregatePair}. The result computes the disjoint sum {@code p1 ⊕ p2}
	 * exactly (by identity (3)), in {@code 4·r} terms.
	 *
	 * <p>This is the correct ASSEMBLY; it is a LOSS (4r vs 2r) until the
	 * corrections are united (Pan identity (4)) — which is the merge of terms
	 * sharing two factor-directions, the open step. Requires {@code |p1|==|p2|}
	 * (cyclic rotations of a base have equal rank).</p>
	 */
	public static List<Term> aggregateDisjoint(
			List<Term> p1, int du1, int dv1, int dw1,
			List<Term> p2, int du2, int dv2, int dw2) {
		if (p1.size() != p2.size()) {
			throw new IllegalArgumentException(
					"term-pairing needs equal rank, got " + p1.size() + " vs " + p2.size());
		}
		int du = du1 + du2, dv = dv1 + dv2, dw = dw1 + dw2;
		List<Term> out = new ArrayList<>(4 * p1.size());
		for (int i = 0; i < p1.size(); i++) {
			Term t1 = p1.get(i), t2 = p2.get(i);
			out.addAll(aggregatePair(
					embed(t1.u(), du, 0), embed(t1.v(), dv, 0), embed(t1.w(), dw, 0),
					embed(t2.u(), du, du1), embed(t2.v(), dv, dv1), embed(t2.w(), dw, dw1)));
		}
		return out;
	}

	private static int[] embed(int[] v, int dim, int offset) {
		int[] r = new int[dim];
		System.arraycopy(v, 0, r, offset, v.length);
		return r;
	}

	/**
	 * Fused rank of Pan-TA on the disjoint cyclic-rotation pair
	 * {@code ⟨n,r,p⟩ ⊕ ⟨p,n,r⟩}: {@code nrp + np + nr + rp}. Derived from the
	 * rotation kin-pairing {@code T1(i,j,k) ↔ T2(k,i,j)} — {@code nrp} aggregated
	 * products plus the three correction families that unite (over j / k / i) into
	 * {@code np / nr / rp} terms. For {@code ⟨26,3,26⟩} this is
	 * {@code 2028 + 676 + 78 + 78 = 2860} = FMM's TA term in ⟨26,29,29⟩=11693.
	 * (Beats the best separate schemes 1504+1504=3008.)
	 */
	public static long fusedRank(int n, int r, int p) {
		return (long) n * r * p + (long) n * p + (long) n * r + (long) r * p;
	}

	/** A fused TA construction over the combined (disjoint-sum) variable space. */
	public record TaScheme(int dimU, int dimV, int dimW, List<Term> terms) {}

	/**
	 * EXACT construction of Pan-TA on {@code ⟨n,r,p⟩ ⊕ ⟨p,n,r⟩} (the two cyclic
	 * cross-blocks). DETERMINISTIC in {@code (n,r,p)} — so replay is just
	 * re-running this; the lineage need only record {@code (n,r,p)} and the two
	 * cross-block leaf refs (bit-exact by construction). The combined space is
	 * {@code A1(n×r)⊕A2(p×n)}, {@code B1(r×p)⊕B2(n×r)}, {@code C1(n×p)⊕C2(p×r)};
	 * the term list computes the block-diagonal product of both — verified by
	 * {@code TestRectangularTrilinearAggregation}. Integer coefficients (±1).
	 */
	public static TaScheme build(int n, int r, int p) {
		int dimU = n * r + p * n;   // A1(n×r) then A2(p×n)
		int dimV = r * p + n * r;   // B1(r×p) then B2(n×r)
		int dimW = n * p + p * r;   // C1(n×p) then C2(p×r)
		List<Term> terms = new ArrayList<>();
		// aggregated: +(A1(i,j)+A2(k,i)) ⊗ (B1(j,k)+B2(i,j)) ⊗ (C1(i,k)+C2(k,j))
		for (int i = 0; i < n; i++) for (int j = 0; j < r; j++) for (int k = 0; k < p; k++) {
			int[] u = new int[dimU]; u[a1(i, j, r)] += 1; u[a2(k, i, n, r, n)] += 1;
			int[] v = new int[dimV]; v[b1(j, k, p)] += 1; v[b2(i, j, r, p, n)] += 1;
			int[] w = new int[dimW]; w[c1(i, k, p)] += 1; w[c2(k, j, r, n, p)] += 1;
			terms.add(new Term(u, v, w));
		}
		// corr1 (unite over j): −A2(k,i) ⊗ Σ_j(B1(j,k)+B2(i,j)) ⊗ C1(i,k), per (i,k)
		for (int i = 0; i < n; i++) for (int k = 0; k < p; k++) {
			int[] u = new int[dimU]; u[a2(k, i, n, r, n)] -= 1;
			int[] v = new int[dimV]; for (int j = 0; j < r; j++) { v[b1(j, k, p)] += 1; v[b2(i, j, r, p, n)] += 1; }
			int[] w = new int[dimW]; w[c1(i, k, p)] += 1;
			terms.add(new Term(u, v, w));
		}
		// corr2 (unite over k): −A1(i,j) ⊗ B2(i,j) ⊗ Σ_k(C1(i,k)+C2(k,j)), per (i,j)
		for (int i = 0; i < n; i++) for (int j = 0; j < r; j++) {
			int[] u = new int[dimU]; u[a1(i, j, r)] -= 1;
			int[] v = new int[dimV]; v[b2(i, j, r, p, n)] += 1;
			int[] w = new int[dimW]; for (int k = 0; k < p; k++) { w[c1(i, k, p)] += 1; w[c2(k, j, r, n, p)] += 1; }
			terms.add(new Term(u, v, w));
		}
		// corr3 (unite over i): −Σ_i(A1(i,j)+A2(k,i)) ⊗ B1(j,k) ⊗ C2(k,j), per (j,k)
		for (int j = 0; j < r; j++) for (int k = 0; k < p; k++) {
			int[] u = new int[dimU]; for (int i = 0; i < n; i++) { u[a1(i, j, r)] -= 1; u[a2(k, i, n, r, n)] -= 1; }
			int[] v = new int[dimV]; v[b1(j, k, p)] += 1;
			int[] w = new int[dimW]; w[c2(k, j, r, n, p)] += 1;
			terms.add(new Term(u, v, w));
		}
		return new TaScheme(dimU, dimV, dimW, terms);
	}

	/**
	 * Assemble the full peeled matmul {@code ⟨N, N+s, N+s⟩} via the TA-fused
	 * cross-pair: peel {@code M=P=N+s} into {@code N+s}, giving four blocks
	 * {@code ⟨N,N,N⟩ (diag→C1)}, {@code ⟨N,s,N⟩ & ⟨N,N,s⟩ (cross, TA-fused)},
	 * {@code ⟨N,s,s⟩ (corner→C2)}. Rank = {@code cube.r + fusedRank(N,s,N) +
	 * corner.r}. The TA terms' combined space is mapped onto the parent A/B/C
	 * sub-blocks (A1=cols[0,N), A2=cols[N,N+s); B11/B12/B21/B22; C1=cols[0,N),
	 * C2=cols[N,N+s)). DETERMINISTIC in {@code (N,s,cube,corner)} — replayable.
	 *
	 * @param cube   a ⟨N,N,N⟩ scheme; @param corner a ⟨N,s,s⟩ scheme
	 * @return a ⟨N, N+s, N+s⟩ scheme (verify with {@code Verifier})
	 */
	public static NonCubicBilinearAlgorithm buildPeeledViaTa(
			int N, int s, NonCubicBilinearAlgorithm cube, NonCubicBilinearAlgorithm corner) {
		if (!(cube.n == N && cube.m == N && cube.p == N)) {
			throw new IllegalArgumentException("cube must be ⟨" + N + "," + N + "," + N + "⟩");
		}
		if (!(corner.n == N && corner.m == s && corner.p == s)) {
			throw new IllegalArgumentException("corner must be ⟨" + N + "," + s + "," + s + "⟩");
		}
		int M = N + s, P = N + s;
		TaScheme ta = build(N, s, N);
		int rank = cube.r + ta.terms().size() + corner.r;
		double[][] U = new double[N * M][rank], V = new double[M * P][rank], W = new double[N * P][rank];
		int col = 0;
		// diag ⟨N,N,N⟩ → A1·B11→C1 (no offsets)
		col = appendBlock(cube, U, V, W, col, M, P, 0, 0, 0, 0, 0, 0);
		// corner ⟨N,s,s⟩ → A2·B22→C2 (A cols +N; B rows+N,cols+N; C cols+N)
		col = appendBlock(corner, U, V, W, col, M, P, 0, N, N, N, 0, N);
		// TA-fused cross-pair, mapped from build's combined space
		for (Term t : ta.terms()) {
			for (int idx = 0; idx < t.u().length; idx++) if (t.u()[idx] != 0) U[taU(idx, N, s, M)][col] += t.u()[idx];
			for (int idx = 0; idx < t.v().length; idx++) if (t.v()[idx] != 0) V[taV(idx, N, s, P)][col] += t.v()[idx];
			for (int idx = 0; idx < t.w().length; idx++) if (t.w()[idx] != 0) W[taW(idx, N, s, P)][col] += t.w()[idx];
			col++;
		}
		return new NonCubicBilinearAlgorithm(N, M, P, U, V, W);
	}

	/** Map a source ⟨a,b,c⟩ scheme's columns into the global space at the given
	 *  A/B/C (row,col) offsets; returns the next free column. */
	private static int appendBlock(NonCubicBilinearAlgorithm src,
			double[][] U, double[][] V, double[][] W, int col, int M, int P,
			int aRowOff, int aColOff, int bRowOff, int bColOff, int cRowOff, int cColOff) {
		double[][] su = src.denseU(), sv = src.denseV(), sw = src.denseW();
		int a = src.n, b = src.m, c = src.p;
		for (int t = 0; t < src.r; t++) {
			for (int i = 0; i < a; i++) for (int j = 0; j < b; j++) {
				double v = su[i * b + j][t];
				if (v != 0) U[(aRowOff + i) * M + (aColOff + j)][col] += v;
			}
			for (int j = 0; j < b; j++) for (int k = 0; k < c; k++) {
				double v = sv[j * c + k][t];
				if (v != 0) V[(bRowOff + j) * P + (bColOff + k)][col] += v;
			}
			for (int i = 0; i < a; i++) for (int k = 0; k < c; k++) {
				double v = sw[i * c + k][t];
				if (v != 0) W[(cRowOff + i) * P + (cColOff + k)][col] += v;
			}
			col++;
		}
		return col;
	}

	// map build(N,s,N) combined-space indices → global ⟨N,N+s,N+s⟩ flattened indices
	private static int taU(int idx, int N, int s, int M) {
		if (idx < N * s) { int i = idx / s, j = idx % s; return i * M + (N + j); }      // A1→A2 cols[N,N+s)
		int d = idx - N * s; int i = d / N, j = d % N; return i * M + j;                // A2→A1 cols[0,N)
	}
	private static int taV(int idx, int N, int s, int P) {
		if (idx < s * N) { int j = idx / N, k = idx % N; return (N + j) * P + k; }      // B1→B21
		int d = idx - s * N; int j = d / s, k = d % s; return j * P + (N + k);          // B2→B12
	}
	private static int taW(int idx, int N, int s, int P) {
		if (idx < N * N) { int i = idx / N, k = idx % N; return i * P + k; }            // C1→C1
		int d = idx - N * N; int i = d / s, k = d % s; return i * P + (N + k);          // C2→C2 cols[N,N+s)
	}

	// flattened indices into the combined space. A2/B2/C2 are offset past block 1.
	private static int a1(int i, int j, int r) { return i * r + j; }                       // A1(n×r)
	private static int a2(int ip, int jp, int n, int r, int nCols) { return n * r + ip * nCols + jp; } // A2(p×n): cols=n
	private static int b1(int j, int k, int p) { return j * p + k; }                       // B1(r×p)
	private static int b2(int jp, int kp, int r, int p, int n) { return r * p + jp * r + kp; }         // B2(n×r): cols=r
	private static int c1(int i, int k, int p) { return i * p + k; }                       // C1(n×p)
	private static int c2(int ip, int kp, int r, int n, int p) { return n * p + ip * r + kp; }         // C2(p×r): cols=r
}
