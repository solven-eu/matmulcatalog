package eu.solven.matmul.search.als;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

/**
 * Enumerate ternary {@code ⟨2,2,2⟩} rank-7 schemes by sweeping the
 * <b>ternary-isotropy orbit</b> of one or more seed schemes, and bucket them by
 * their <b>canonical recombination multiset</b>.
 *
 * <h2>Why an orbit sweep rather than blind brute force</h2>
 * A blind covering-DFS over {@code ({-1,0,1}⁴)³} branches ~39k ways per term with
 * no effective pruning at the top levels ({@code 39k⁷ ≈ 10³²}) — measured
 * infeasible (see {@link Ternary2x2x2MultisetEnumerator}). De Groote (1978)
 * proves <b>all</b> rank-7 {@code ⟨2,2,2⟩} decompositions form a single
 * {@code GL₂(ℚ)³} orbit, so every ternary scheme is a change-of-basis image of
 * Strassen. We therefore enumerate the orbit directly under <b>ternary</b>
 * changes of basis {@code (X,Y,Z)} (each an invertible {@code 2×2} matrix with
 * entries in {@code {-1,0,1}}), keeping only the images that stay ternary.
 *
 * <p>Isotropy action (sandwiching {@code Â=XAY, B̂=Y⁻¹BZ ⟹ Ĉ=XCZ}), per product
 * {@code k}, with {@code U_k,V_k,W_k} read as {@code 2×2} matrices:
 * <pre>
 *   U'_k = Xᵀ · U_k · Yᵀ
 *   V'_k = Y⁻ᵀ · V_k · Zᵀ
 *   W'_k = X⁻¹ · W_k · Z⁻¹
 * </pre>
 *
 * <h2>Optimality tier</h2>
 * The result is <b>exhaustive within the ternary-isotropy scope</b>
 * (optimal-within-scope, per CLAUDE.md): every ternary scheme reachable from a
 * seed by a ternary {@code (X,Y,Z)}. It is <i>not</i> certified-exhaustive over
 * all ternary schemes (some might need a rational {@code (X,Y,Z)}). Robustness is
 * evidenced by seeding from every catalog scheme and by optionally widening the
 * change-of-basis alphabet; if the canonical-multiset set is stable, that is
 * strong evidence of completeness.
 */
public final class Ternary2x2x2Orbit {

	private Ternary2x2x2Orbit() {}

	/** A {@code 2×2} integer change-of-basis matrix and its (rational) inverse. */
	private static final class Mat2 {
		final int[] a;       // row-major [a00,a01,a10,a11]
		final double[] inv;  // inverse, row-major
		Mat2(int[] a) {
			this.a = a;
			int det = a[0] * a[3] - a[1] * a[2];
			this.inv = new double[] { a[3] / (double) det, -a[1] / (double) det,
					-a[2] / (double) det, a[0] / (double) det };
		}
	}

	/** All invertible ternary {@code 2×2} matrices. */
	private static List<Mat2> invertibleTernary(int lo, int hi) {
		List<Mat2> out = new ArrayList<>();
		for (int a = lo; a <= hi; a++)
			for (int b = lo; b <= hi; b++)
				for (int c = lo; c <= hi; c++)
					for (int d = lo; d <= hi; d++)
						if (a * d - b * c != 0) out.add(new Mat2(new int[] { a, b, c, d }));
		return out;
	}

	// ---- 2×2 helpers (row-major double[4]) ----------------------------------

	private static double[] mmId(int[] x) { // int → double
		return new double[] { x[0], x[1], x[2], x[3] };
	}

	private static double[] transpose(double[] m) {
		return new double[] { m[0], m[2], m[1], m[3] };
	}

	private static double[] mul(double[] a, double[] b) {
		return new double[] {
				a[0] * b[0] + a[1] * b[2], a[0] * b[1] + a[1] * b[3],
				a[2] * b[0] + a[3] * b[2], a[2] * b[1] + a[3] * b[3] };
	}

	/** Round to nearest int and verify exactly ternary; return null otherwise. */
	private static int[] toTernary(double[] m) {
		int[] out = new int[4];
		for (int i = 0; i < 4; i++) {
			long r = Math.round(m[i]);
			if (Math.abs(m[i] - r) > 1e-7 || r < -1 || r > 1) return null;
			out[i] = (int) r;
		}
		return out;
	}

	/**
	 * Apply isotropy {@code (X,Y,Z)} to a seed scheme {@code seed[7][3][4]}
	 * (per product: {@code [u,v,w]}, each a row-major {@code 2×2}). Returns the
	 * transformed scheme if it is ternary everywhere, else {@code null}.
	 */
	static int[][][] transform(int[][][] seed, Mat2 X, Mat2 Y, Mat2 Z) {
		double[] Xt = transpose(mmId(X.a)), Yt = transpose(mmId(Y.a)), Zt = transpose(mmId(Z.a));
		double[] Xinv = X.inv, Zinv = Z.inv, YinvT = transpose(Y.inv);
		int[][][] out = new int[7][3][];
		for (int k = 0; k < 7; k++) {
			int[] uTer = toTernary(mul(mul(Xt, mmId(seed[k][0])), Yt));
			if (uTer == null) return null;
			int[] vTer = toTernary(mul(mul(YinvT, mmId(seed[k][1])), Zt));
			if (vTer == null) return null;
			int[] wTer = toTernary(mul(mul(Xinv, mmId(seed[k][2])), Zinv));
			if (wTer == null) return null;
			out[k][0] = uTer; out[k][1] = vTer; out[k][2] = wTer;
		}
		return out;
	}

	// ---- canonical multiset --------------------------------------------------

	/** The 6 permutations of 3 coordinate positions (the scheme's S₃ symmetry). */
	private static final int[][] PERMS3 = {
			{ 0, 1, 2 }, { 0, 2, 1 }, { 1, 0, 2 }, { 1, 2, 0 }, { 2, 0, 1 }, { 2, 1, 0 } };

	private static final int BIG = 9, SMALL = 8;

	private static int maxBlk(int bits) { return (bits & 1) != 0 ? BIG : SMALL; }

	/** Per-product sub-shapes at the (9,8)³ reference split. */
	static int[][] shapes(int[][][] scheme) {
		int[][] shapes = new int[7][3];
		for (int k = 0; k < 7; k++) {
			int[] u = scheme[k][0], v = scheme[k][1], w = scheme[k][2];
			int uRowN = 0, uColM = 0, vRowM = 0, vColP = 0, wRowN = 0, wColP = 0;
			for (int a = 0; a < 4; a++) if (u[a] != 0) { uRowN |= 1 << (a >> 1); uColM |= 1 << (a & 1); }
			for (int b = 0; b < 4; b++) if (v[b] != 0) { vRowM |= 1 << (b >> 1); vColP |= 1 << (b & 1); }
			for (int c = 0; c < 4; c++) if (w[c] != 0) { wRowN |= 1 << (c >> 1); wColP |= 1 << (c & 1); }
			shapes[k][0] = Math.min(maxBlk(uRowN), maxBlk(wRowN));
			shapes[k][1] = Math.min(maxBlk(uColM), maxBlk(vRowM));
			shapes[k][2] = Math.min(maxBlk(vColP), maxBlk(wColP));
		}
		return shapes;
	}

	/**
	 * Canonical multiset key: minimise the sorted bag of sub-shapes over the 6
	 * global axis-permutations (the scheme's S₃ symmetry). Two schemes share this
	 * key iff their recombination multisets are equal up to rotating/transposing
	 * the whole scheme.
	 */
	static String canonicalMultisetKey(int[][][] scheme) {
		int[][] base = shapes(scheme);
		String best = null;
		for (int[] perm : PERMS3) {
			String[] s = new String[7];
			for (int k = 0; k < 7; k++) {
				int[] sh = base[k];
				s[k] = sh[perm[0]] + "," + sh[perm[1]] + "," + sh[perm[2]];
			}
			Arrays.sort(s);
			String key = String.join("|", s);
			if (best == null || key.compareTo(best) < 0) best = key;
		}
		return best;
	}

	/** maxBlk under a per-axis block-flip: when flipped, the "big" block is block1. */
	private static int maxBlkFlip(int bits, boolean flip) {
		int bigBit = flip ? 2 : 1; // bit0=block0, bit1=block1
		return (bits & bigBit) != 0 ? BIG : SMALL;
	}

	/** Per-product sub-shapes under axis block-flip mask (bit0=A/n, bit1=B/m, bit2=C/p). */
	private static int[][] shapesMasked(int[][][] scheme, int flipMask) {
		boolean fA = (flipMask & 1) != 0, fB = (flipMask & 2) != 0, fC = (flipMask & 4) != 0;
		int[][] out = new int[7][3];
		for (int k = 0; k < 7; k++) {
			int[] u = scheme[k][0], v = scheme[k][1], w = scheme[k][2];
			int uRowN = 0, uColM = 0, vRowM = 0, vColP = 0, wRowN = 0, wColP = 0;
			for (int a = 0; a < 4; a++) if (u[a] != 0) { uRowN |= 1 << (a >> 1); uColM |= 1 << (a & 1); }
			for (int b = 0; b < 4; b++) if (v[b] != 0) { vRowM |= 1 << (b >> 1); vColP |= 1 << (b & 1); }
			for (int c = 0; c < 4; c++) if (w[c] != 0) { wRowN |= 1 << (c >> 1); wColP |= 1 << (c & 1); }
			out[k][0] = Math.min(maxBlkFlip(uRowN, fA), maxBlkFlip(wRowN, fA));
			out[k][1] = Math.min(maxBlkFlip(uColM, fB), maxBlkFlip(vRowM, fB));
			out[k][2] = Math.min(maxBlkFlip(vColP, fC), maxBlkFlip(wColP, fC));
		}
		return out;
	}

	/**
	 * Fully canonical multiset key: minimise over the scheme's S₃ axis-role
	 * permutations <b>and</b> the per-axis block-flip masks (Z₂³ — relabelling
	 * which of {@code n₁,n₂} is the "first" block). This is the strongest
	 * equivalence: two schemes share this key iff their recombination multisets
	 * agree up to any rotation/transpose of the scheme and any reordering of the
	 * 2-part split on each axis.
	 */
	static String flipCanonicalMultisetKey(int[][][] scheme) {
		String best = null;
		for (int mask = 0; mask < 8; mask++) {
			int[][] base = shapesMasked(scheme, mask);
			for (int[] perm : PERMS3) {
				String[] s = new String[7];
				for (int k = 0; k < 7; k++) {
					int[] sh = base[k];
					s[k] = sh[perm[0]] + "," + sh[perm[1]] + "," + sh[perm[2]];
				}
				Arrays.sort(s);
				String key = String.join("|", s);
				if (best == null || key.compareTo(best) < 0) best = key;
			}
		}
		return best;
	}

	/** Pretty counted form of a canonical key, e.g. {@code 1·⟨8,8,8⟩ + 3·⟨8,8,9⟩ …}. */
	static String pretty(String canonicalKey) {
		String[] parts = canonicalKey.split("\\|");
		Map<String, Integer> counts = new TreeMap<>();
		for (String p : parts) counts.merge(p, 1, Integer::sum);
		StringBuilder sb = new StringBuilder();
		boolean first = true;
		for (var e : counts.entrySet()) {
			if (!first) sb.append(" + ");
			sb.append(e.getValue()).append("·⟨").append(e.getKey()).append('⟩');
			first = false;
		}
		return sb.toString();
	}

	// ---- driver --------------------------------------------------------------

	public static final class Result {
		/** canonical multiset key -> representative scheme. */
		public final Map<String, int[][][]> representatives = new LinkedHashMap<>();
		/** canonical multiset key -> how many distinct ternary schemes realise it. */
		public final Map<String, Long> schemeCounts = new LinkedHashMap<>();
		public long transformsTried;
		public long ternaryHits;
	}

	/**
	 * Sweep the ternary-isotropy orbit of every seed and collect distinct
	 * canonical multisets.
	 *
	 * @param seeds          seed schemes ({@code [7][3][4]})
	 * @param cobAlphabetLo  change-of-basis entry lower bound (e.g. -1)
	 * @param cobAlphabetHi  upper bound (e.g. +1)
	 */
	public static Result sweep(List<int[][][]> seeds, int cobAlphabetLo, int cobAlphabetHi) {
		return sweep(seeds, cobAlphabetLo, cobAlphabetHi, true);
	}

	/**
	 * @param countDistinctSchemes when {@code true}, dedups schemes (memory ∝ #
	 *        distinct ternary schemes) and fills {@code schemeCounts}; when
	 *        {@code false}, only collects distinct canonical multisets — O(1)
	 *        memory, suitable for wide change-of-basis alphabets.
	 */
	public static Result sweep(List<int[][][]> seeds, int cobAlphabetLo, int cobAlphabetHi,
			boolean countDistinctSchemes) {
		List<Mat2> cob = invertibleTernary(cobAlphabetLo, cobAlphabetHi);
		Result res = new Result();
		java.util.Set<String> seenScheme = countDistinctSchemes ? new java.util.HashSet<>() : null;
		long progressEvery = 50_000_000L;
		long start = System.nanoTime();
		for (int[][][] seed : seeds) {
			for (Mat2 X : cob)
				for (Mat2 Y : cob)
					for (Mat2 Z : cob) {
						res.transformsTried++;
						if (res.transformsTried % progressEvery == 0) {
							long ms = (System.nanoTime() - start) / 1_000_000;
							System.out.printf("[progress] transforms=%,d ternary-hits=%,d distinct-multisets=%d %,d ms%n",
									res.transformsTried, res.ternaryHits, res.representatives.size(), ms);
							System.out.flush();
						}
						int[][][] img = transform(seed, X, Y, Z);
						if (img == null) continue;
						res.ternaryHits++;
						if (seenScheme != null && !seenScheme.add(signature(img))) continue;
						String key = canonicalMultisetKey(img);
						res.schemeCounts.merge(key, 1L, Long::sum);
						res.representatives.putIfAbsent(key, img);
					}
		}
		return res;
	}

	/** A ternary scheme realising a target multiset, plus the exact change-of-basis from the seed. */
	public static final class Realization {
		public final int[][][] scheme;
		public final int[] X, Y, Z; // row-major 2×2 change-of-basis matrices
		Realization(int[][][] scheme, int[] X, int[] Y, int[] Z) {
			this.scheme = scheme; this.X = X; this.Y = Y; this.Z = Z;
		}
	}

	/**
	 * Find a ternary scheme in the isotropy orbit of {@code seed} whose canonical
	 * multiset equals {@code targetCanonicalKey}, returning it with the exact
	 * {@code (X,Y,Z)} change-of-basis that produced it (for a reproducible
	 * lineage). Returns {@code null} if none within the alphabet.
	 */
	public static Realization findRealization(int[][][] seed, String targetCanonicalKey, int cobAlphabet) {
		List<Mat2> cob = invertibleTernary(-cobAlphabet, cobAlphabet);
		for (Mat2 X : cob)
			for (Mat2 Y : cob)
				for (Mat2 Z : cob) {
					int[][][] img = transform(seed, X, Y, Z);
					if (img == null) continue;
					if (canonicalMultisetKey(img).equals(targetCanonicalKey))
						return new Realization(img, X.a.clone(), Y.a.clone(), Z.a.clone());
				}
		return null;
	}

	/** Stable signature of a scheme (for distinct-scheme counting). */
	private static String signature(int[][][] s) {
		StringBuilder sb = new StringBuilder();
		for (int[][] term : s)
			for (int[] vec : term)
				for (int x : vec) sb.append(x + 1); // 0/1/2, single char
		return sb.toString();
	}

	// ========================================================================
	//  Symbolic / exact enumeration over the FULL GL₂(ℚ)³ orbit
	// ========================================================================
	//
	//  By de Groote (1978) every rank-7 ⟨2,2,2⟩ scheme is (X,Y,Z)·seed for some
	//  X,Y,Z ∈ GL₂(ℚ). The recombination sub-dimension on each axis depends on
	//  ONLY ONE of X,Y,Z, and only through its two COLUMN DIRECTIONS:
	//
	//    subAₖ = BIG  ⟺  d₀ᵀUₖ ≠ 0  ∧  perp(d₁)ᵀWₖ ≠ 0      (X = [d₀ | d₁])
	//    subBₖ = BIG  ⟺  Uₖ·e₀ ≠ 0  ∧  perp(e₁)ᵀVₖ ≠ 0      (Y rows e₀,e₁)
	//    subCₖ = BIG  ⟺  Vₖ·f₀ ≠ 0  ∧  Wₖ·perp(f₁) ≠ 0       (Z rows f₀,f₁)
	//
	//  The pattern is piecewise-constant in the directions, jumping only when a
	//  direction equals a null-space of a fixed integer Strassen factor. So a
	//  finite direction set covering all integer null-spaces PLUS one generic
	//  direction enumerates EVERY realizable pattern — exactly, with no sampling.
	//  Axes are independent, so the full multiset set is the (zipped-by-product)
	//  combination of the three per-axis pattern sets.

	private static boolean leftNZ(int[] d, int[] M) { // dᵀ·M ≠ 0  (M row-major 2×2)
		return (d[0] * M[0] + d[1] * M[2]) != 0 || (d[0] * M[1] + d[1] * M[3]) != 0;
	}

	private static boolean rightNZ(int[] M, int[] e) { // M·e ≠ 0
		return (M[0] * e[0] + M[1] * e[1]) != 0 || (M[2] * e[0] + M[3] * e[1]) != 0;
	}

	private static int[] perp(int[] d) { return new int[] { -d[1], d[0] }; }

	/** Projective-distinct integer directions with entries in [-bound, bound]. */
	private static List<int[]> directions(int bound) {
		java.util.Set<String> seen = new java.util.HashSet<>();
		List<int[]> out = new ArrayList<>();
		for (int x = -bound; x <= bound; x++)
			for (int y = -bound; y <= bound; y++) {
				if (x == 0 && y == 0) continue;
				int g = gcd(Math.abs(x), Math.abs(y));
				int cx = x / g, cy = y / g;
				if (cx < 0 || (cx == 0 && cy < 0)) { cx = -cx; cy = -cy; } // sign-canonical
				if (seen.add(cx + "," + cy)) out.add(new int[] { cx, cy });
			}
		return out;
	}

	private static int gcd(int a, int b) { return b == 0 ? Math.max(a, 1) : gcd(b, a % b); }

	private static boolean sameDir(int[] a, int[] b) { return a[0] * b[1] - a[1] * b[0] == 0; }

	/** The seed's 7 factor matrices as row-major {@code int[4]} (U=index0, V=1, W=2). */
	private static int[][][] factors(int[][][] seed) {
		int[][][] f = new int[3][7][4];
		for (int k = 0; k < 7; k++)
			for (int t = 0; t < 3; t++)
				f[t][k] = seed[k][t].clone();
		return f;
	}

	/** Add every realizable 7-bit axis pattern to {@code out}, sweeping two directions. */
	private interface BitFn { boolean big(int[] dir0, int[] dir1, int k); }

	private static void collectPatterns(List<int[]> D, BitFn fn, java.util.Set<Integer> out) {
		for (int[] d0 : D)
			for (int[] d1 : D) {
				if (sameDir(d0, d1)) continue; // columns/rows must be independent
				int pat = 0;
				for (int k = 0; k < 7; k++) if (fn.big(d0, d1, k)) pat |= 1 << k;
				out.add(pat);
			}
	}

	/**
	 * Exact, certified-complete set of canonical recombination multisets over the
	 * entire {@code GL₂(ℚ)³} orbit of {@code seed} (= all rank-7 ⟨2,2,2⟩ schemes).
	 * @param dirBound direction-coordinate bound (≥2 to include a generic direction).
	 */
	public static java.util.Set<String> certifiedQMultisets(int[][][] seed, int dirBound) {
		List<int[]> D = directions(dirBound);
		int[][][] f = factors(seed);
		int[][] U = f[0], V = f[1], W = f[2];

		java.util.Set<Integer> nPat = new java.util.LinkedHashSet<>();
		java.util.Set<Integer> mPat = new java.util.LinkedHashSet<>();
		java.util.Set<Integer> pPat = new java.util.LinkedHashSet<>();
		collectPatterns(D, (d0, d1, k) -> leftNZ(d0, U[k]) && leftNZ(perp(d1), W[k]), nPat);
		collectPatterns(D, (e0, e1, k) -> rightNZ(U[k], e0) && leftNZ(perp(e1), V[k]), mPat);
		collectPatterns(D, (f0, f1, k) -> rightNZ(V[k], f0) && rightNZ(W[k], perp(f1)), pPat);

		java.util.Set<String> multisets = new java.util.LinkedHashSet<>();
		for (int np : nPat)
			for (int mp : mPat)
				for (int pp : pPat) {
					int[][] shapes = new int[7][3];
					for (int k = 0; k < 7; k++) {
						shapes[k][0] = ((np >> k) & 1) != 0 ? BIG : SMALL;
						shapes[k][1] = ((mp >> k) & 1) != 0 ? BIG : SMALL;
						shapes[k][2] = ((pp >> k) & 1) != 0 ? BIG : SMALL;
					}
					multisets.add(canonicalMultisetKeyOfShapes(shapes));
				}
		return multisets;
	}

	/** S₃-canonical key directly from a shape array (same convention as {@link #canonicalMultisetKey}). */
	static String canonicalMultisetKeyOfShapes(int[][] base) {
		String best = null;
		for (int[] perm : PERMS3) {
			String[] s = new String[7];
			for (int k = 0; k < 7; k++) {
				int[] sh = base[k];
				s[k] = sh[perm[0]] + "," + sh[perm[1]] + "," + sh[perm[2]];
			}
			Arrays.sort(s);
			String key = String.join("|", s);
			if (best == null || key.compareTo(best) < 0) best = key;
		}
		return best;
	}

	/** Per-axis pattern-set sizes (diagnostics). */
	public static int[] perAxisPatternCounts(int[][][] seed, int dirBound) {
		List<int[]> D = directions(dirBound);
		int[][][] f = factors(seed);
		int[][] U = f[0], V = f[1], W = f[2];
		java.util.Set<Integer> nPat = new java.util.LinkedHashSet<>();
		java.util.Set<Integer> mPat = new java.util.LinkedHashSet<>();
		java.util.Set<Integer> pPat = new java.util.LinkedHashSet<>();
		collectPatterns(D, (d0, d1, k) -> leftNZ(d0, U[k]) && leftNZ(perp(d1), W[k]), nPat);
		collectPatterns(D, (e0, e1, k) -> rightNZ(U[k], e0) && leftNZ(perp(e1), V[k]), mPat);
		collectPatterns(D, (f0, f1, k) -> rightNZ(V[k], f0) && rightNZ(W[k], perp(f1)), pPat);
		return new int[] { nPat.size(), mPat.size(), pPat.size() };
	}
}
