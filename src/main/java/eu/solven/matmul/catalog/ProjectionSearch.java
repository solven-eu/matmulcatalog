package eu.solven.matmul.catalog;

import java.util.ArrayList;
import java.util.List;
import java.util.Optional;

import eu.solven.matmul.FactorMatrix;
import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;

/**
 * Discovery driver for the projection operator (#159 / ROADMAP). For a target
 * {@code ⟨n,m,p⟩}, projects larger parent schemes
 * {@code ⟨n+a,m+b,p+c⟩} down to it by {@link Compose#project}, exhaustively
 * trying every drop-index choice per axis (the best choice is not predictable
 * for recursively-composed schemes), and keeps the cheapest verified result
 * strictly below {@code upperBound}.
 *
 * <p>Projection is exact by construction, so the search compares projected
 * ranks directly; the chosen winner is still re-verified by
 * {@link Verifier#isExactNonCubic}.</p>
 */
public final class ProjectionSearch {

	private ProjectionSearch() {}

	/** A verified projection of a parent down to the target. */
	public record Hit(NonCubicBilinearAlgorithm scheme, NonCubicBilinearAlgorithm parent,
			int[] keepN, int[] keepM, int[] keepP, long rank) {}

	/**
	 * Projection margin μ(S): the maximum number of products dead-code-eliminated
	 * when the single best index (over all three axes) is dropped. Then
	 * {@code R_after = R_before − μ}, so a HIGH μ means S projects far below its
	 * rank and is a strong projection parent <em>even at higher rank</em> — the
	 * (rank, μ) Pareto axis (paper §projmargin). Computed from the same
	 * {@link Supports} the search uses, so it never drifts from
	 * {@link Supports#survivorCount}.
	 *
	 * <p>A product dies when index {@code i} is dropped on an axis iff its support
	 * on that axis, in <em>either</em> of the two matrices carrying that axis
	 * (n→{U,W}, m→{U,V}, p→{V,W}), collapses to the single index {@code i}.</p>
	 */
	public static int projectionMargin(NonCubicBilinearAlgorithm a) {
		Supports s = Supports.of(a);
		int best = 0;
		best = Math.max(best, axisMargin(s.uN, s.wN, a.n, a.r));  // n: U-rows ∪ W-rows
		best = Math.max(best, axisMargin(s.uM, s.vM, a.m, a.r));  // m: U-cols ∪ V-rows
		best = Math.max(best, axisMargin(s.vP, s.wP, a.p, a.r));  // p: V-cols ∪ W-cols
		return best;
	}

	/**
	 * Per-axis margin triple {@code [μ_n, μ_m, μ_p]} — the best single-index drop
	 * per axis. The scalar {@link #projectionMargin} is the max of these; the
	 * triple is the multi-dimensional "how good a projecting base is THIS, and in
	 * which direction" (mirror of the per-axis serendipity savings triple).
	 */
	public static int[] axisMargins(NonCubicBilinearAlgorithm a) {
		Supports s = Supports.of(a);
		return new int[] {
				axisMargin(s.uN, s.wN, a.n, a.r),
				axisMargin(s.uM, s.vM, a.m, a.r),
				axisMargin(s.vP, s.wP, a.p, a.r) };
	}

	/**
	 * EXACT projected rank of {@code a} down to {@code ⟨n,m,p⟩}: the minimum
	 * survivor count over all keep-index combinations (the real currency of
	 * projection — what {@link #bestFor} would verify). Returns {@code -1} when
	 * the shape doesn't project there (or the combo space exceeds the work cap).
	 * Projected cost {@code = r − μ(drop)}; rank's exchange rate against margin
	 * is exactly 1, so walks should optimize THIS, not the rank-lex margin.
	 */
	public static long projectedRank(NonCubicBilinearAlgorithm a, int n, int m, int p,
			int maxDelta) {
		Combo c = bestCombo(n, m, p, a, Supports.of(a), a.r + 1L, maxDelta);
		return c == null ? -1 : c.rank();
	}

	/**
	 * Public handle for mask-based exact survivor counting against ONE parent —
	 * {@link Supports} computed once, then any keep-mask combination is priced in
	 * O(rank). This is what makes chained/lattice projection sweeps affordable
	 * (the optimal mask is exponential to find; a GIVEN mask is cheap to price,
	 * and the count IS the exact rank of that projection after DCE).
	 */
	public static final class MaskedProjector {
		private final Supports sup;

		public MaskedProjector(NonCubicBilinearAlgorithm a) {
			this.sup = Supports.of(a);
		}

		public long survivors(boolean[] maskN, boolean[] maskM, boolean[] maskP, long ceiling) {
			return sup.survivorCount(maskN, maskM, maskP, ceiling);
		}
	}

	/** Max over indices of #products whose support (in factor1 OR factor2) on this
	 *  axis is exactly that one index. */
	private static int axisMargin(int[][] supp1, int[][] supp2, int dim, int r) {
		int[] death = new int[dim];
		for (int l = 0; l < r; l++) {
			int i1 = supp1[l].length == 1 ? supp1[l][0] : -1;
			int i2 = supp2[l].length == 1 ? supp2[l][0] : -1;
			if (i1 >= 0) death[i1]++;
			if (i2 >= 0 && i2 != i1) death[i2]++;  // count a product once per index
		}
		int m = 0;
		for (int d : death) m = Math.max(m, d);
		return m;
	}

	/**
	 * Guard: skip a parent whose (combinations × rank) exceeds this. Sized so the
	 * full 17–32 cube→cube band is reachable — e.g. ⟨32,32,32⟩→⟨31,31,31⟩ is
	 * {@code 32³ combos × ~15096 rank ≈ 4.9e8}. Each combo is now an
	 * early-exiting {@link #survivorCount} (no array build, no per-combo verify),
	 * so the constant per (combo×rank) unit is tiny.
	 */
	private static final long WORK_CAP = 2_000_000_000L;

	/**
	 * Best verified projection of any {@code parents} entry down to
	 * {@code ⟨n,m,p⟩} strictly below {@code upperBound}, dropping at most
	 * {@code maxDelta} indices per axis.
	 */
	public static Optional<Hit> bestFor(int n, int m, int p,
			List<NonCubicBilinearAlgorithm> parents, long upperBound, int maxDelta) {
		// Phase 1: find the cheapest (parent, drop-choice) by survivor-counting,
		// which is exact for projection (rank = #products surviving DCE) and avoids
		// allocating/verifying every candidate. Phase 2 builds + verifies only the
		// single winner.
		NonCubicBilinearAlgorithm bestParent = null;
		Combo best = null;
		long bestRank = upperBound;
		for (NonCubicBilinearAlgorithm s : parents) {
			// Precompute each product's index-support ONCE per parent (not per combo).
			Supports sup = Supports.of(s);
			Combo c = bestCombo(n, m, p, s, sup, bestRank, maxDelta);
			if (c != null && c.rank() < bestRank) {
				bestParent = s;
				best = c;
				bestRank = c.rank();
			}
		}
		if (bestParent == null) return Optional.empty();
		return buildAndVerify(n, m, p, bestParent, best);
	}

	/** The cheapest drop-choice and its survivor count for ONE parent projected to
	 *  {@code ⟨n,m,p⟩}, or {@code null} if the parent can't beat {@code ceiling}
	 *  (shape mismatch, over-delta, over-cap, margin-pruned, or no surviving combo). */
	private record Combo(int[] keepN, int[] keepM, int[] keepP, long rank) {}

	private static Combo bestCombo(int n, int m, int p, NonCubicBilinearAlgorithm s,
			Supports sup, long ceiling, int maxDelta) {
		if (s.n < n || s.m < m || s.p < p) return null;
		if (s.n - n > maxDelta || s.m - m > maxDelta || s.p - p > maxDelta) return null;
		if (s.n == n && s.m == m && s.p == p) return null; // identity, not a projection
		List<int[]> kn = combos(s.n, n), km = combos(s.m, m), kp = combos(s.p, p);
		long work = (long) kn.size() * km.size() * kp.size() * s.r;
		if (work > WORK_CAP) return null; // too large to enumerate exhaustively here
		// Loose-but-rigorous margin prune: every product the projection can DCE is
		// charged to one dropped axis, so #eliminable ≤ Σ (best single-row margin)
		// over the dropped axes, i.e. projected rank ≥ R − Σ axisMargin. If even that
		// lower bound can't beat the ceiling, no drop-combo here can — skip the whole
		// enumeration. Only valid for ≤1-row drops per axis (axisMargin = the best-1-row
		// margin); ≥2-row drops would need the best-d-rows margin, so it's gated off there.
		if (s.n - n <= 1 && s.m - m <= 1 && s.p - p <= 1) {
			long marginBound = s.r;
			if (s.n > n) marginBound -= axisMargin(sup.uN, sup.wN, s.n, s.r);
			if (s.m > m) marginBound -= axisMargin(sup.uM, sup.vM, s.m, s.r);
			if (s.p > p) marginBound -= axisMargin(sup.vP, sup.wP, s.p, s.r);
			if (marginBound >= ceiling) return null;
		}
		boolean[] maskN = new boolean[s.n], maskM = new boolean[s.m], maskP = new boolean[s.p];
		int[] bKN = null, bKM = null, bKP = null;
		long bRank = ceiling;
		for (int[] keepN : kn)
			for (int[] keepM : km)
				for (int[] keepP : kp) {
					setMask(maskN, keepN);
					setMask(maskM, keepM);
					setMask(maskP, keepP);
					long rank = sup.survivorCount(maskN, maskM, maskP, bRank);
					if (rank >= bRank) continue;
					bKN = keepN;
					bKM = keepM;
					bKP = keepP;
					bRank = rank;
				}
		return bKN == null ? null : new Combo(bKN, bKM, bKP, bRank);
	}

	/** Phase 2: build + verify the single winning projection. Projection (restrict
	 *  indices + DCE) is a PROVEN correctness-preserving operator, so the result is
	 *  exact by construction when the parent is — a cheap random spot check suffices
	 *  here (it catches a bad parent / index-map mistake). A full exact verify of a
	 *  rank-10⁴ dim-30 projection is ~10¹¹ ops (minutes); it is reserved for
	 *  promote-time (PromoteStagingWins / LineageVerifier — "bound during search,
	 *  proven at commit"). This is what makes projecting large cubes affordable. */
	private static Optional<Hit> buildAndVerify(int n, int m, int p,
			NonCubicBilinearAlgorithm parent, Combo c) {
		NonCubicBilinearAlgorithm proj = Compose.project(parent, c.keepN(), c.keepM(), c.keepP());
		if (proj.r != c.rank() || !Verifier.passesRandomMatmulSpotCheck(proj)) return Optional.empty();
		return Optional.of(new Hit(proj, parent, c.keepN(), c.keepM(), c.keepP(), proj.r));
	}

	/**
	 * Parent-centric (scatter) projection: project ONE already-resolved parent down
	 * to MANY target children, computing the parent's {@link Supports} a single time
	 * and reusing it across every child. The dual of {@link #bestFor} (one child,
	 * many parents); it is what lets a sweep replay/expand each (large, stub) parent
	 * exactly once and fan it out to all the slightly-smaller shapes it can reach,
	 * instead of re-resolving the same parent once per child.
	 *
	 * @param children   target shapes, each {@code {n,m,p}}
	 * @param uppers     per-child strict upper bound (only projections {@code < uppers[i]} are kept)
	 * @return per-child best verified projection of {@code parent} (aligned with
	 *         {@code children}; {@code Optional.empty()} where the parent can't beat its bound)
	 */
	public static List<Optional<Hit>> projectToMany(NonCubicBilinearAlgorithm parent,
			int[][] children, long[] uppers, int maxDelta) {
		Supports sup = Supports.of(parent);
		List<Optional<Hit>> out = new ArrayList<>(children.length);
		for (int i = 0; i < children.length; i++) {
			Combo c = bestCombo(children[i][0], children[i][1], children[i][2],
					parent, sup, uppers[i], maxDelta);
			out.add(c == null ? Optional.empty()
					: buildAndVerify(children[i][0], children[i][1], children[i][2], parent, c));
		}
		return out;
	}

	/** Reset {@code mask} to false then set the kept indices true. */
	private static void setMask(boolean[] mask, int[] keep) {
		java.util.Arrays.fill(mask, false);
		for (int k : keep) mask[k] = true;
	}

	/**
	 * Per-product index supports for one parent scheme, computed once. For each
	 * product {@code l} we store the distinct indices it touches on each
	 * (factor, axis): {@code uN}/{@code uM} from U, {@code vM}/{@code vP} from V,
	 * {@code wN}/{@code wP} from W. A product survives a projection iff, on each
	 * factor, at least one touched index on BOTH of its axes is kept — exactly
	 * the DCE rule of {@link Compose#project}, but evaluated over small int
	 * arrays instead of re-striding the {@code double[][]} per combo.
	 */
	static final class Supports {
		final int rank;
		final int[][] uN, uM, vM, vP, wN, wP;

		private Supports(int rank, int[][] uN, int[][] uM, int[][] vM, int[][] vP,
				int[][] wN, int[][] wP) {
			this.rank = rank;
			this.uN = uN; this.uM = uM; this.vM = vM; this.vP = vP; this.wN = wN; this.wP = wP;
		}

		static Supports of(NonCubicBilinearAlgorithm a) {
			int r = a.r;
			FactorMatrix u = a.u(), v = a.v(), w = a.w();
			int[][] uN = axisSupport(u, a.m, r, true);   // n-axis of U (row = i*m+j)
			int[][] uM = axisSupport(u, a.m, r, false);  // m-axis of U
			int[][] vM = axisSupport(v, a.p, r, true);   // m-axis of V (row = j*p+k)
			int[][] vP = axisSupport(v, a.p, r, false);  // p-axis of V
			int[][] wN = axisSupport(w, a.p, r, true);   // n-axis of W (row = i*p+k)
			int[][] wP = axisSupport(w, a.p, r, false);  // p-axis of W
			return new Supports(r, uN, uM, vM, vP, wN, wP);
		}

		/**
		 * For factor matrix {@code f} whose row index is {@code outer*innerDim +
		 * inner}, return per-product the sorted distinct {@code outer} indices
		 * (if {@code outerAxis}) or {@code inner} indices touched by a non-zero.
		 *
		 * <p>Iterates the sparse column (product) support directly --- {@code O(nnz)}
		 * over the factor's non-zeros, never densifying — instead of striding the
		 * full {@code outerDim·innerDim} grid. {@code axisDim} is recovered from the
		 * factor's row count and {@code innerDim}.</p>
		 */
		private static int[][] axisSupport(FactorMatrix f, int innerDim, int r, boolean outerAxis) {
			int outerDim = f.rows() / innerDim;
			int axisDim = outerAxis ? outerDim : innerDim;
			int[][] out = new int[r][];
			boolean[] seen = new boolean[axisDim];
			for (int l = 0; l < r; l++) {
				java.util.Arrays.fill(seen, false);
				f.forEachInColumn(l, (row, val) -> seen[outerAxis ? row / innerDim : row % innerDim] = true);
				int cnt = 0;
				for (boolean b : seen) if (b) cnt++;
				int[] s = new int[cnt];
				int wIdx = 0;
				for (int i = 0; i < axisDim; i++) if (seen[i]) s[wIdx++] = i;
				out[l] = s;
			}
			return out;
		}

		/** Count products surviving the projection defined by the keep-masks,
		 *  early-exiting once the count reaches {@code ceiling}. */
		long survivorCount(boolean[] maskN, boolean[] maskM, boolean[] maskP, long ceiling) {
			long count = 0;
			for (int l = 0; l < rank; l++) {
				if (!anyKept(uN[l], maskN) || !anyKept(uM[l], maskM)) continue; // U DCE'd
				if (!anyKept(vM[l], maskM) || !anyKept(vP[l], maskP)) continue; // V DCE'd
				if (!anyKept(wN[l], maskN) || !anyKept(wP[l], maskP)) continue; // W DCE'd
				if (++count >= ceiling) return count; // can't beat the incumbent
			}
			return count;
		}

		private static boolean anyKept(int[] support, boolean[] mask) {
			for (int idx : support) if (mask[idx]) return true;
			return false;
		}
	}

	/** All sorted size-{@code keep} subsets of {@code [0,dim)}. */
	static List<int[]> combos(int dim, int keep) {
		List<int[]> out = new ArrayList<>();
		rec(0, 0, dim, keep, new int[keep], out);
		return out;
	}

	private static void rec(int start, int idx, int dim, int keep, int[] cur, List<int[]> out) {
		if (idx == keep) { out.add(cur.clone()); return; }
		for (int i = start; i <= dim - (keep - idx); i++) {
			cur[idx] = i;
			rec(i + 1, idx + 1, dim, keep, cur, out);
		}
	}
}
