package eu.solven.matmul;

import eu.solven.matmul.recombination.BlockSplitSearch;

import java.util.ArrayList;
import java.util.List;

/**
 * Symmetry transforms on bilinear algorithms — the action of the symmetric
 * group {@code S₃} on tensor slots. Generators:
 *
 * <ul>
 *   <li>{@link NonCubicBilinearAlgorithm#cyclicShift} — order-3 cyclic
 *       relabelling {@code ⟨n,m,p⟩ → ⟨m,p,n⟩}</li>
 *   <li>{@link NonCubicBilinearAlgorithm#transpose} — order-2 transpose
 *       involution {@code ⟨n,m,p⟩ → ⟨p,m,n⟩}</li>
 * </ul>
 *
 * <p>DIS09 §3 (Layer 3) uses these "{@code S_{X,Y,Z}(U,V,W)}" transforms to
 * generate equivalent algorithms with different U/V/W column-support
 * distributions per axis. When such a variant is used as an outer base in
 * a recursive composition with min-reduction over sub-product shapes
 * (e.g. {@link eu.solven.matmul.recombination.BlockSplitSearch#findBestMultiBaseSplit}),
 * each variant can yield a different total rank — the framework picks the
 * best.</p>
 */
public final class SymmetryTransforms {

	private SymmetryTransforms() {}

	/**
	 * Returns the {@code S₃} orbit of {@code alg} under tensor-slot
	 * relabelling. Up to 6 distinct algorithms (fewer if {@code alg} has
	 * a non-trivial stabiliser in {@code S₃}, e.g. cubic algorithms that
	 * are cyclic-symmetric).
	 *
	 * <p>Dedup is by-content over (n, m, p, U, V, W); two orbit elements
	 * with identical factor matrices are returned only once.</p>
	 */
	/** One S₃ orbit element together with the EXACT axis-relabel ("ABC-&gt;XYZ")
	 *  that produced it from the canonical scheme. Unlike a size-matched
	 *  {@code Transpose} perm, this records the precise permutation even when the
	 *  base has two equal-sized axes (⟨4,4,3⟩ from ⟨3,4,4⟩) — which is what makes
	 *  the recorded lineage reproduce the SAME orientation the search scored
	 *  (otherwise replay re-orients ambiguously → a worse multiset → the
	 *  predict/build divergence guard fires). */
	public record S3Variant(NonCubicBilinearAlgorithm alg, String perm) {}

	/**
	 * Like {@link #s3Orbit} but each variant carries the exact axis-relabel
	 * "ABC-&gt;XYZ" that produced it (deterministic composition of
	 * {@link NonCubicBilinearAlgorithm#transpose()} /
	 * {@link NonCubicBilinearAlgorithm#cyclicShift()}). Same traversal ORDER and
	 * by-content dedup as {@link #s3Orbit}, so the per-shape representative is
	 * identical — only now the orientation is recordable as
	 * {@link eu.solven.matmul.catalog.Lineage.Transpose}, replayed bit-exactly by
	 * {@link #permuteAxes}.
	 */
	public static List<S3Variant> s3OrbitWithPerms(NonCubicBilinearAlgorithm alg) {
		List<S3Variant> out = new ArrayList<>(6);
		List<String> seen = new ArrayList<>(6);
		NonCubicBilinearAlgorithm cur = alg;
		String t = "ABC"; // current axis triple relative to canonical ABC
		for (int c = 0; c < 3; c++) {
			addVariantIfNew(cur, "ABC->" + t, out, seen);
			addVariantIfNew(cur.transpose(), "ABC->" + reverse3(t), out, seen);
			cur = cur.cyclicShift();
			t = rotateLeft3(t);
		}
		return out;
	}

	private static void addVariantIfNew(NonCubicBilinearAlgorithm a, String perm,
			List<S3Variant> out, List<String> seen) {
		String sig = signature(a);
		if (!seen.contains(sig)) {
			seen.add(sig);
			out.add(new S3Variant(a, perm));
		}
	}

	private static String reverse3(String s) {
		return "" + s.charAt(2) + s.charAt(1) + s.charAt(0);
	}

	private static String rotateLeft3(String s) {
		return "" + s.charAt(1) + s.charAt(2) + s.charAt(0);
	}

	/**
	 * Apply the axis-relabel {@code perm} ("ABC-&gt;XYZ" or "NMP-&gt;XYZ") to
	 * {@code alg}, producing the EXACT variant by composing
	 * {@link NonCubicBilinearAlgorithm#transpose()} /
	 * {@link NonCubicBilinearAlgorithm#cyclicShift()} — bit-exact and unambiguous
	 * even for equal-sized axes (unlike a shape-based {@code orientAs}). Returns
	 * {@code null} when {@code perm} is not a parseable bijection of {A,B,C}
	 * (e.g. a legacy non-bijective string), so callers can fall back.
	 */
	public static NonCubicBilinearAlgorithm permuteAxes(NonCubicBilinearAlgorithm alg, String perm) {
		int arrow = perm.indexOf("->");
		if (arrow < 0 || perm.length() - arrow - 2 < 3) return null;
		String tgt = perm.substring(arrow + 2, arrow + 5)
				.replace('N', 'A').replace('M', 'B').replace('P', 'C');
		return switch (tgt) {
			case "ABC" -> alg;
			case "CBA" -> alg.transpose();
			case "BCA" -> alg.cyclicShift();
			case "CAB" -> alg.cyclicShift().cyclicShift();
			case "ACB" -> alg.cyclicShift().transpose();
			case "BAC" -> alg.transpose().cyclicShift();
			default -> null; // non-bijective / unrecognised
		};
	}

	public static List<NonCubicBilinearAlgorithm> s3Orbit(NonCubicBilinearAlgorithm alg) {
		List<NonCubicBilinearAlgorithm> out = new ArrayList<>(6);
		List<String> seen = new ArrayList<>(6);

		NonCubicBilinearAlgorithm cur = alg;
		for (int c = 0; c < 3; c++) {
			addIfNew(cur, out, seen);
			addIfNew(cur.transpose(), out, seen);
			cur = cur.cyclicShift();
		}
		return out;
	}

	/**
	 * Axis-flip orbit: up to 8 variants from independently reversing
	 * the index order on each of the three matrix axes. DIS09 §3
	 * S_{X,Y,Z} transform restricted to the row-reverse-J subgroup —
	 * the cheap, shape-uniform subset of the full per-axis-permutation
	 * group (which scales as n!·m!·p! and is intractable for big n).
	 * Each variant computes the same matmul but assigns the r rank-1
	 * products to different output positions; used as outer templates
	 * in recombination, they give different sub-shape allocations on
	 * the same target.
	 */
	public static List<NonCubicBilinearAlgorithm> axisFlipOrbit(NonCubicBilinearAlgorithm alg) {
		List<NonCubicBilinearAlgorithm> out = new ArrayList<>(8);
		List<String> seen = new ArrayList<>(8);
		for (int mask = 0; mask < 8; mask++) {
			boolean swapA = (mask & 1) != 0;
			boolean swapB = (mask & 2) != 0;
			boolean swapC = (mask & 4) != 0;
			addIfNew(applyAxisFlip(alg, swapA, swapB, swapC), out, seen);
		}
		return out;
	}

	/**
	 * Same as {@link #axisFlipOrbit} but excludes the canonical (mask=0)
	 * scheme — returns ONLY the flipped variants. Used as an A/B probe:
	 * running the search with this orbit measures whether axis-flipped
	 * outer bases find candidates the canonical scheme misses, without
	 * letting the canonical scheme "win by default".
	 */
	public static List<NonCubicBilinearAlgorithm> axisFlipOrbitExcludingCanonical(NonCubicBilinearAlgorithm alg) {
		List<NonCubicBilinearAlgorithm> out = new ArrayList<>(7);
		List<String> seen = new ArrayList<>(7);
		// Seed seen with the canonical signature so it's not re-added if some
		// flipped variant collapses back to it (e.g. fully-symmetric Strassen).
		seen.add(signature(alg));
		for (int mask = 1; mask < 8; mask++) {
			boolean swapA = (mask & 1) != 0;
			boolean swapB = (mask & 2) != 0;
			boolean swapC = (mask & 4) != 0;
			addIfNew(applyAxisFlip(alg, swapA, swapB, swapC), out, seen);
		}
		return out;
	}

	private static NonCubicBilinearAlgorithm applyAxisFlip(NonCubicBilinearAlgorithm alg,
			boolean swapA, boolean swapB, boolean swapC) {
		int a = alg.n, b = alg.m, c = alg.p, r = alg.r;
		double[][] srcU = alg.denseU();
		double[][] srcV = alg.denseV();
		double[][] srcW = alg.denseW();
		double[][] U2 = new double[a * b][r];
		double[][] V2 = new double[b * c][r];
		double[][] W2 = new double[a * c][r];
		for (int i = 0; i < a; i++) for (int j = 0; j < b; j++) {
			int iP = swapA ? (a - 1 - i) : i;
			int jP = swapB ? (b - 1 - j) : j;
			for (int k = 0; k < r; k++) U2[iP * b + jP][k] = srcU[i * b + j][k];
		}
		for (int j = 0; j < b; j++) for (int l = 0; l < c; l++) {
			int jP = swapB ? (b - 1 - j) : j;
			int lP = swapC ? (c - 1 - l) : l;
			for (int k = 0; k < r; k++) V2[jP * c + lP][k] = srcV[j * c + l][k];
		}
		for (int i = 0; i < a; i++) for (int l = 0; l < c; l++) {
			int iP = swapA ? (a - 1 - i) : i;
			int lP = swapC ? (c - 1 - l) : l;
			for (int k = 0; k < r; k++) W2[iP * c + lP][k] = srcW[i * c + l][k];
		}
		return new NonCubicBilinearAlgorithm(a, b, c, U2, V2, W2);
	}

	/**
	 * Full cheap orbit: every distinct (U, V, W) variant produced by
	 * composing {@link #s3Orbit} (shape relabeling) with
	 * {@link #axisFlipOrbit} (axis-row-reverse). Upper bound 6·8 = 48
	 * variants; in practice 4-24 after content-dedup. Canonical
	 * "all the cheap variants of a scheme" enumeration for pool
	 * construction in {@code BlockSplitSearch}.
	 */
	public static List<NonCubicBilinearAlgorithm> fullCheapOrbit(NonCubicBilinearAlgorithm alg) {
		List<NonCubicBilinearAlgorithm> out = new ArrayList<>();
		List<String> seen = new ArrayList<>();
		for (NonCubicBilinearAlgorithm shape : s3Orbit(alg)) {
			for (NonCubicBilinearAlgorithm flipped : axisFlipOrbit(shape)) {
				addIfNew(flipped, out, seen);
			}
		}
		return out;
	}

	/**
	 * Three "internal" orbit modes, in increasing coverage / cost,
	 * matching the DIS09 §3 {@code S_{X,Y,Z}(U, V, W)} hierarchy
	 * restricted to the search-feasible subgroups.
	 *
	 * <ul>
	 *   <li>{@link #CANONICAL} — the scheme as-stored, no orbit
	 *       expansion. 1 variant.</li>
	 *   <li>{@link #AXIS_FLIP} — DIS09 with {@code X, Y, Z ∈ {I, J}}
	 *       (J = anti-diagonal). At most {@code 2³ = 8} variants
	 *       <em>regardless of shape</em>. Always tractable. Coincides
	 *       with the full permutation orbit for ⟨2,2,2⟩ schemes.</li>
	 *   <li>{@link #PERMUTATION_BOUNDED} — DIS09 with full per-axis
	 *       permutation matrices, {@code (n! · m! · p!)} variants per
	 *       scheme. Tractable for ⟨3,3,3⟩ (216 variants) and below;
	 *       borderline at ⟨4,4,4⟩ (13 824); intractable at ⟨5,5,5⟩
	 *       (≈ 1.7M). Auto-falls-back to {@link #AXIS_FLIP} when the
	 *       per-scheme variant count exceeds a configurable cap.</li>
	 * </ul>
	 *
	 * <p>The DIS09 framework's most general form takes
	 * {@code (X, Y, Z) ∈ GL_n × GL_m × GL_p} — a continuous space, not
	 * enumerable. Probert-Fischer 1980 and similar constructions
	 * USE specific invertibles from the algebra of the problem, but
	 * even DIS09's own search restricted to permutation matrices
	 * — so the {@link #PERMUTATION_BOUNDED} mode is the practical
	 * upper limit, not the theoretical one.</p>
	 *
	 * <p><strong>Caveat for practice</strong>: even
	 * {@code PERMUTATION_BOUNDED} on ⟨4,4,4⟩ adds 13 824 variants per
	 * scheme — multiplied across a pool of {@code k} schemes that's
	 * {@code 13.8k · k} pool entries. Recombination search cost is
	 * roughly linear in pool size, so this multiplies wall-clock by
	 * the same factor. Useful for one-off audits ("can we close
	 * gap X with full Laderman orbit?") but not for routine sweeps.
	 * We may find experimentally that the extra variants almost never
	 * surface a better result, in which case the cap should stay low.</p>
	 */
	public enum InternalOrbitMode {
		CANONICAL,
		AXIS_FLIP,
		AXIS_FLIP_ONLY,
		PERMUTATION_BOUNDED;
	}

	/**
	 * Dispatch to one of the {@link InternalOrbitMode} generators, with
	 * an optional cap on the variant count above which
	 * {@code PERMUTATION_BOUNDED} falls back to {@code AXIS_FLIP}.
	 *
	 * @param maxVariants ignored unless {@code mode == PERMUTATION_BOUNDED};
	 *                    must be ≥ 8 to make the fallback meaningful.
	 *                    {@code Long.MAX_VALUE} means "no cap".
	 */
	public static List<NonCubicBilinearAlgorithm> internalOrbit(NonCubicBilinearAlgorithm alg,
			InternalOrbitMode mode, long maxVariants) {
		return switch (mode) {
			case CANONICAL -> List.of(alg);
			case AXIS_FLIP -> axisFlipOrbit(alg);
			case AXIS_FLIP_ONLY -> axisFlipOrbitExcludingCanonical(alg);
			case PERMUTATION_BOUNDED -> permutationOrbit(alg, maxVariants);
		};
	}

	/**
	 * Per-axis permutation orbit: enumerate every
	 * {@code (P_x, P_y, P_z) ∈ S_n × S_m × S_p} and apply DIS09's
	 * {@code S_{P_x, P_y, P_z}(U, V, W)}. Bounded variant count
	 * {@code n! · m! · p!}; auto-falls-back to
	 * {@link #axisFlipOrbit} if the count exceeds {@code maxVariants}.
	 * Dedup by full {@code (U, V, W)} content signature — many
	 * permutations collapse to the same factor matrices for symmetric
	 * schemes (e.g. Strassen).
	 */
	public static List<NonCubicBilinearAlgorithm> permutationOrbit(NonCubicBilinearAlgorithm alg,
			long maxVariants) {
		long count = (long) factorial(alg.n) * factorial(alg.m) * factorial(alg.p);
		if (count > maxVariants) {
			// Fallback: too many to enumerate. Caller asked for the bounded
			// version, so honour the bound by stepping down rather than
			// throwing.
			return axisFlipOrbit(alg);
		}
		List<int[]> permsX = allPermutations(alg.n);
		List<int[]> permsY = allPermutations(alg.m);
		List<int[]> permsZ = allPermutations(alg.p);
		List<NonCubicBilinearAlgorithm> out = new ArrayList<>();
		List<String> seen = new ArrayList<>();
		for (int[] px : permsX)
			for (int[] py : permsY)
				for (int[] pz : permsZ)
					addIfNew(applyAxisPerm(alg, px, py, pz), out, seen);
		return out;
	}

	/** An orbit variant together with the EXACT per-axis permutation that produced it
	 *  from the canonical scheme. Axis-flips are reversal permutations, so this captures
	 *  flips AND permutations uniformly — letting a lineage record
	 *  {@code AxisPermute(canonical, permX, permY, permZ)} that replays to the variant
	 *  precisely (task #91: the base's exact orientation is no longer lost). */
	public record PermutedVariant(NonCubicBilinearAlgorithm alg, int[] permX, int[] permY, int[] permZ) {}

	/** {@link #internalOrbit} but each variant carries the per-axis permutation that
	 *  reconstructs it from {@code alg} (variant 0 is always the canonical with identity
	 *  perms). Used by the pool builder to record a precise, replayable base orientation. */
	public static List<PermutedVariant> internalOrbitWithPerms(NonCubicBilinearAlgorithm alg,
			InternalOrbitMode mode, long maxVariants) {
		return switch (mode) {
			case CANONICAL -> List.of(new PermutedVariant(alg, identityPerm(alg.n), identityPerm(alg.m), identityPerm(alg.p)));
			case AXIS_FLIP -> flipOrbitWithPerms(alg, false);
			case AXIS_FLIP_ONLY -> flipOrbitWithPerms(alg, true);
			case PERMUTATION_BOUNDED -> permutationOrbitWithPerms(alg, maxVariants);
		};
	}

	private static List<PermutedVariant> flipOrbitWithPerms(NonCubicBilinearAlgorithm alg,
			boolean excludeCanonical) {
		List<PermutedVariant> out = new ArrayList<>(8);
		List<String> seen = new ArrayList<>(8);
		for (int mask = 0; mask < 8; mask++) {
			if (excludeCanonical && mask == 0) continue;
			boolean swapA = (mask & 1) != 0, swapB = (mask & 2) != 0, swapC = (mask & 4) != 0;
			NonCubicBilinearAlgorithm v = applyAxisFlip(alg, swapA, swapB, swapC);
			String key = signature(v);
			if (seen.contains(key)) continue;
			seen.add(key);
			out.add(new PermutedVariant(v, reversalOrIdentity(alg.n, swapA),
					reversalOrIdentity(alg.m, swapB), reversalOrIdentity(alg.p, swapC)));
		}
		return out;
	}

	private static List<PermutedVariant> permutationOrbitWithPerms(NonCubicBilinearAlgorithm alg,
			long maxVariants) {
		long count = (long) factorial(alg.n) * factorial(alg.m) * factorial(alg.p);
		if (count > maxVariants) {
			return flipOrbitWithPerms(alg, false);  // bounded fallback (flips are reversal perms)
		}
		List<int[]> permsX = allPermutations(alg.n);
		List<int[]> permsY = allPermutations(alg.m);
		List<int[]> permsZ = allPermutations(alg.p);
		List<PermutedVariant> out = new ArrayList<>();
		List<String> seen = new ArrayList<>();
		for (int[] px : permsX)
			for (int[] py : permsY)
				for (int[] pz : permsZ) {
					NonCubicBilinearAlgorithm v = applyAxisPerm(alg, px, py, pz);
					String key = signature(v);
					if (seen.contains(key)) continue;
					seen.add(key);
					out.add(new PermutedVariant(v, px.clone(), py.clone(), pz.clone()));
				}
		return out;
	}

	private static int[] identityPerm(int n) {
		int[] a = new int[n];
		for (int i = 0; i < n; i++) a[i] = i;
		return a;
	}

	private static int[] reversalOrIdentity(int n, boolean reverse) {
		int[] a = new int[n];
		for (int i = 0; i < n; i++) a[i] = reverse ? (n - 1 - i) : i;
		return a;
	}

	private static int factorial(int n) {
		int r = 1;
		for (int i = 2; i <= n; i++) r *= i;
		return r;
	}

	private static List<int[]> allPermutations(int n) {
		List<int[]> out = new ArrayList<>(factorial(n));
		int[] cur = new int[n];
		for (int i = 0; i < n; i++) cur[i] = i;
		permRec(cur, 0, out);
		return out;
	}

	private static void permRec(int[] cur, int k, List<int[]> out) {
		if (k == cur.length - 1) {
			out.add(cur.clone());
			return;
		}
		for (int i = k; i < cur.length; i++) {
			int t = cur[k]; cur[k] = cur[i]; cur[i] = t;
			permRec(cur, k + 1, out);
			t = cur[k]; cur[k] = cur[i]; cur[i] = t;
		}
	}

	private static NonCubicBilinearAlgorithm applyAxisPerm(NonCubicBilinearAlgorithm alg,
			int[] px, int[] py, int[] pz) {
		int a = alg.n, b = alg.m, c = alg.p, r = alg.r;
		double[][] srcU = alg.denseU();
		double[][] srcV = alg.denseV();
		double[][] srcW = alg.denseW();
		double[][] U2 = new double[a * b][r];
		double[][] V2 = new double[b * c][r];
		double[][] W2 = new double[a * c][r];
		for (int i = 0; i < a; i++) for (int j = 0; j < b; j++)
			for (int k = 0; k < r; k++)
				U2[px[i] * b + py[j]][k] = srcU[i * b + j][k];
		for (int j = 0; j < b; j++) for (int l = 0; l < c; l++)
			for (int k = 0; k < r; k++)
				V2[py[j] * c + pz[l]][k] = srcV[j * c + l][k];
		for (int i = 0; i < a; i++) for (int l = 0; l < c; l++)
			for (int k = 0; k < r; k++)
				W2[px[i] * c + pz[l]][k] = srcW[i * c + l][k];
		return new NonCubicBilinearAlgorithm(a, b, c, U2, V2, W2);
	}

	private static void addIfNew(NonCubicBilinearAlgorithm a,
			List<NonCubicBilinearAlgorithm> out, List<String> seen) {
		String key = signature(a);
		if (!seen.contains(key)) {
			seen.add(key);
			out.add(a);
		}
	}

	/**
	 * Apply axis-flip mask {@code (mask&1, mask&2, mask&4)} = (swapA, swapB, swapC)
	 * to {@code alg}. Mask=0 returns the canonical (no flips); mask=7 flips all
	 * three axes. Use this when a caller has chosen a mask analytically (e.g.
	 * via {@link eu.solven.matmul.recombination.AnalyticalMaskSearch}) and now needs
	 * the actual scheme to materialise.
	 */
	public static NonCubicBilinearAlgorithm applyAxisFlipMask(NonCubicBilinearAlgorithm alg, int mask) {
		boolean sA = (mask & 1) != 0;
		boolean sB = (mask & 2) != 0;
		boolean sC = (mask & 4) != 0;
		return applyAxisFlip(alg, sA, sB, sC);
	}

	/**
	 * Canonical axis-flip representative of {@code alg}: the orbit element
	 * with the lex-minimum {@link #signature}. Two schemes that are
	 * axis-flip variants of each other have the same canonical form.
	 *
	 * <p>Use case: dedup pool entries that are axis-flip-equivalent. With
	 * {@link eu.solven.matmul.recombination.AnalyticalMaskSearch} at search time,
	 * keeping only the canonical avoids redundancy (the search re-derives
	 * the 8 axis-flip variants analytically per allocation).
	 */
	public static NonCubicBilinearAlgorithm axisFlipCanonical(NonCubicBilinearAlgorithm alg) {
		NonCubicBilinearAlgorithm best = alg;
		String bestSig = signature(alg);
		for (int mask = 1; mask < 8; mask++) {
			boolean sA = (mask & 1) != 0;
			boolean sB = (mask & 2) != 0;
			boolean sC = (mask & 4) != 0;
			NonCubicBilinearAlgorithm flipped = applyAxisFlip(alg, sA, sB, sC);
			String sig = signature(flipped);
			if (sig.compareTo(bestSig) < 0) {
				best = flipped;
				bestSig = sig;
			}
		}
		return best;
	}

	/** Lex-min signature over the axis-flip orbit — useful as a hash for dedup. */
	public static String axisFlipCanonicalSignature(NonCubicBilinearAlgorithm alg) {
		return signature(axisFlipCanonical(alg));
	}

	/**
	 * True iff {@code a} and {@code b} are in the same axis-flip orbit. Two
	 * schemes that pass this test produce the same set of shape multisets
	 * across all (allocation, mask) pairs — keeping both in the pool is
	 * cost-redundant if mask enumeration is done at search time.
	 */
	public static boolean axisFlipEquivalent(NonCubicBilinearAlgorithm a, NonCubicBilinearAlgorithm b) {
		if (a.n != b.n || a.m != b.m || a.p != b.p || a.r != b.r) return false;
		return axisFlipCanonicalSignature(a).equals(axisFlipCanonicalSignature(b));
	}

	private static String signature(NonCubicBilinearAlgorithm a) {
		double[][] srcU = a.denseU();
		double[][] srcV = a.denseV();
		double[][] srcW = a.denseW();
		StringBuilder sb = new StringBuilder();
		sb.append(a.n).append(',').append(a.m).append(',').append(a.p).append('|');
		appendMatrix(sb, srcU);
		sb.append('|');
		appendMatrix(sb, srcV);
		sb.append('|');
		appendMatrix(sb, srcW);
		return sb.toString();
	}

	private static void appendMatrix(StringBuilder sb, double[][] M) {
		for (double[] row : M) {
			for (double v : row) {
				sb.append(v).append(',');
			}
			sb.append(';');
		}
	}
}
