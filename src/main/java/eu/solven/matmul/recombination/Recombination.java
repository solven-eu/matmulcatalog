package eu.solven.matmul.recombination;

import eu.solven.matmul.catalog.KnownAlgorithmCatalog;
import eu.solven.matmul.catalog.FieldAwareLookup;

import eu.solven.matmul.catalog.KnownAlgorithm;

import eu.solven.matmul.search.CitedBound;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Optional;

import eu.solven.matmul.NonCubicBilinearAlgorithm;

/**
 * Java port of DeepMind AlphaTensor's recombination algorithm
 * ({@code recombination.py}), generalizing Sedoglavic 2017 and Drevet–Islam–Schost 2009.
 *
 * <p>Given a target matmul problem {@code ⟨a, b, c⟩} and a base algorithm
 * for {@code ⟨base_a, base_b, base_c⟩}, this searches over all ways to
 * allocate target dimensions onto the base algorithm's structure, finds the
 * best resulting decomposition, and reports the rank cost (using a SOTA
 * resolver to look up best-known ranks of the sub-problems).</p>
 *
 * <p>Unlike pure Kronecker product (see {@link Compose}), recombination can
 * find compositions with non-cubic intermediate sub-problems and exploit
 * "addition redundancy" via {@link #processAdditions} — the size of a sum
 * `M_a + M_b` only needs to be as big as the larger summand, not their
 * concatenation. This sometimes yields lower total rank than `r_outer ·
 * r_inner` would.</p>
 *
 * <p><b>Provenance</b>: faithful port from
 * <a href="https://github.com/google-deepmind/alphatensor/blob/main/recombination/recombination.py">alphatensor/recombination/recombination.py</a>
 * (Apache 2.0).</p>
 */
public final class Recombination {

	private Recombination() {}

	/** Lookup: best known rank for {@code ⟨a, b, c⟩} over some field. */
	@FunctionalInterface
	public interface SotaResolver {
		/**
		 * Sentinel for "no scheme in catalog", historically returned by lookups. Chosen as
		 * {@code Integer.MAX_VALUE / 100} (not {@code Integer.MAX_VALUE}) so that summing a
		 * handful of them in a bound does NOT overflow {@code int}.
		 *
		 * <p><b>Contract:</b> a {@link SotaResolver} (best <i>achievable</i> rank) must NEVER
		 * return this — it must fall back to the naive {@code a·b·c} (always constructible). A
		 * raw sota that returns this poisons sum-based lower bounds
		 * ({@code AllocationOptimizer}'s root LB sums getRank over relaxation blocks): one
		 * sentinel makes the LB ≥ any incumbent and a good base is dropped (the ⟨12,13,13⟩=1274
		 * loss). {@code FieldAwareLookup.findRank} now returns naive instead of this; the
		 * constant remains the canonical "unknown" marker for the few callers that test it.</p>
		 */
		int UNKNOWN_RANK = Integer.MAX_VALUE / 100;

		int getRank(int a, int b, int c);
	}

	/** Catalog-backed resolver; falls back to the naive {@code a·b·c} bound if no entry exists. */
	/**
	 * SOTA sub-rank oracle for recombination costing. For NON-commutative algebras this consults the
	 * FULL on-disk catalog ({@link FieldAwareLookup#findRank}), memoised — NOT the small curated
	 * {@link KnownAlgorithmCatalog}, which only carries a hand-picked starting set and so returned the
	 * CUBIC product {@code a·b·c} for any shape it lacked (e.g. ⟨9,9,9⟩→729 instead of the catalog's
	 * 486). That cubic fallback silently inflated every recombination cost and quietly broke
	 * frontier/allocation experiments. Commutative algebras stay on the curated set, since
	 * {@code FieldAwareLookup} is non-commutative-oriented and would mis-serve them.
	 */
	public static SotaResolver catalogResolver(eu.solven.matmul.algebra.Algebra algebra) {
		if (algebra.commutative()) {
			return (a, b, c) -> {
				if (a == 0 || b == 0 || c == 0) return 0;
				return KnownAlgorithmCatalog.bestKnown(a, b, c, algebra).map(k -> k.rank).orElse(a * b * c);
			};
		}
		FieldAwareLookup lookup = new FieldAwareLookup(algebra.field());
		java.util.Map<Long, Integer> memo = new java.util.concurrent.ConcurrentHashMap<>();
		return (a, b, c) -> {
			if (a == 0 || b == 0 || c == 0) return 0;
			if (a == 1) return b * c;
			if (b == 1) return a * c;
			if (c == 1) return a * b;
			return memo.computeIfAbsent((long) a * 1_000_000L + b * 1_000L + c, k -> {
				int v = lookup.findRank(a, b, c);
				return v >= SotaResolver.UNKNOWN_RANK ? a * b * c : v;
			});
		};
	}

	public static final class Result {
		public final long totalRank;
		/** {@code [3][]}: per-axis allocation across base slots (a, b, c). */
		public final int[][] allocations;
		/** {@code [base_rank][3]}: dimensions of the sub-problem at each base-rank position. */
		public final int[][] smallMatrixSizes;

		public Result(long totalRank, int[][] allocations, int[][] smallMatrixSizes) {
			this.totalRank = totalRank;
			this.allocations = allocations;
			this.smallMatrixSizes = smallMatrixSizes;
		}

		@Override
		public String toString() {
			return String.format("rank=%d alloc=[%s,%s,%s] subSizes=%s",
					totalRank,
					Arrays.toString(allocations[0]),
					Arrays.toString(allocations[1]),
					Arrays.toString(allocations[2]),
					Arrays.deepToString(smallMatrixSizes));
		}
	}

	/**
	 * Decompose {@code ⟨targetA, targetB, targetC⟩} using {@code base} on the
	 * first level of recursion. Searches over all allocations of target
	 * dimensions onto the base's {@code (base_a, base_b, base_c)} slots and
	 * returns the configuration with minimal total rank under {@code sota}.
	 *
	 * <p>Note: the returned {@link Result} describes the best decomposition's
	 * shape (allocation + per-sub-problem sizes); it does <b>not</b> build the
	 * composed factor matrices. To realise an actual algorithm, recursively
	 * recombine each sub-problem then compose constructively (future work).</p>
	 */
	/**
	 * Score one specific allocation triple, without searching. Useful to
	 * force a non-degenerate split (e.g. {@code [4, 3]} on a base
	 * {@code ⟨2,2,2⟩}) instead of letting the search collapse to
	 * {@code [0, 7]} when the catalog has a strong direct entry for the
	 * full target. Allocations must sum to {@code (targetA, targetB, targetC)}.
	 */
	public static Result recombineWithAllocation(NonCubicBilinearAlgorithm base, SotaResolver sota,
			int[] allocA, int[] allocB, int[] allocC) {
		return recombineWithAllocation(base, sota, allocA, allocB, allocC, null, null, null);
	}

	/**
	 * Output-zero-peel-aware variant of {@link #recombineWithAllocation}.
	 *
	 * <p>When the outer split is on a PADDED target (e.g. {@code ⟨17,17,17⟩}
	 * padded to {@code ⟨18,18,18⟩} via Strassen ⟨2,2,2⟩ on {@code [9,9]³}),
	 * some blocks of A/B/C are larger than what the original problem
	 * needs — the extra rows/columns are padding zeros that will be
	 * peeled off the result. The {@code peelA/B/C} arrays give the
	 * per-block count of those wasted units. Sub-products that read
	 * from or write to a padded block can be computed at the smaller
	 * effective shape, since the wasted positions contribute nothing
	 * to the final un-padded result. This is the
	 * <strong>Islam 2009 MSc Ch. 4 γ5 reduction</strong> — see
	 * {@code docs/PEELING_ZEROS.md} for the plain-English explanation
	 * and {@code docs/diagnostics/17x17x17_pair_fuse.md} for the
	 * numerical accounting that lands at {@code R(⟨17,17,17⟩) ≤ 2934}
	 * via this mechanism.
	 *
	 * <p>If {@code peelA == null} the legacy semantics apply (no peel,
	 * sum of {@code allocA} must equal the target axis size).
	 *
	 * @param allocA over-allocation of axis A (block sizes including padding)
	 * @param peelA  per-block padding count to peel off (0 = full block,
	 *               {@code allocA[i]} = entirely-padded block); must
	 *               have the same length as {@code allocA} and satisfy
	 *               {@code 0 ≤ peelA[i] ≤ allocA[i]} elementwise. Null
	 *               means zero peel everywhere.
	 */
	/**
	 * Fast variant of {@link #recombineWithAllocation} that reuses
	 * pre-extracted {@link eu.solven.matmul.recombination.AnalyticalMaskSearch.SchemeSupports}
	 * instead of re-scanning the U/V/W factor matrices on every call.
	 *
	 * <p>Cost per call drops from {@code O(r · 3 · dim²)} (full factor scan)
	 * to {@code O(r · dim)} (one per-product support lookup + max over
	 * effective-allocation values). For an outer base with r=93 products
	 * on a 5×5×5 scheme this is roughly a 15× speedup in the inner loop
	 * of {@link eu.solven.matmul.recombination.BlockSplitSearch#findBestMultiBaseSplit}.
	 *
	 * <p>Caller must extract the supports ONCE per base (outside the alloc
	 * loops) and pass the same reference on every call. The supports are
	 * read-only and thread-safe.
	 *
	 * @param supports pre-extracted block-support sets for {@code base}.
	 *                 Must satisfy {@code supports.n == base.n}, etc.
	 */
	public static Result recombineWithAllocationFast(
			NonCubicBilinearAlgorithm base,
			eu.solven.matmul.recombination.AnalyticalMaskSearch.SchemeSupports supports,
			SotaResolver sota,
			int[] allocA, int[] allocB, int[] allocC,
			int[] peelA, int[] peelB, int[] peelC) {
		int baseRank = base.r;
		int[] effA = applyPeel(allocA, peelA);
		int[] effB = applyPeel(allocB, peelB);
		int[] effC = applyPeel(allocC, peelC);
		int[][] shapes = eu.solven.matmul.recombination.AnalyticalMaskSearch.shapesAt(
				supports, effA, effB, effC);
		long totalRank = 0;
		for (int r = 0; r < baseRank; r++) {
			totalRank += sota.getRank(shapes[r][0], shapes[r][1], shapes[r][2]);
		}
		return new Result(totalRank,
				new int[][] { allocA.clone(), allocB.clone(), allocC.clone() },
				shapes);
	}

	public static Result recombineWithAllocation(NonCubicBilinearAlgorithm base, SotaResolver sota,
			int[] allocA, int[] allocB, int[] allocC,
			int[] peelA, int[] peelB, int[] peelC) {
		int baseA = base.n, baseB = base.m, baseC = base.p, baseRank = base.r;
		if (allocA.length != baseA || allocB.length != baseB || allocC.length != baseC) {
			throw new IllegalArgumentException("allocation length must equal base dim");
		}
		// Effective (peeled) per-block sizes. Used wherever the sub-product
		// touches a padded block — both INPUT side (input padding reads as
		// zero, so the linear-form across blocks effectively narrows) and
		// OUTPUT side (output padding entries get peeled, so they don't
		// need to be computed).
		int[] effA = applyPeel(allocA, peelA);
		int[] effB = applyPeel(allocB, peelB);
		int[] effC = applyPeel(allocC, peelC);
		long totalRank = 0;
		int[][] subSizes = new int[baseRank][3];
		for (int r = 0; r < baseRank; r++) {
			int[] uu = processAdditions(base.denseU(), r, baseA, baseB, effA, effB);
			int[] vv = processAdditions(base.denseV(), r, baseB, baseC, effB, effC);
			int[] ww = processAdditions(base.denseW(), r, baseA, baseC, effA, effC);
			int subA = Math.min(uu[0], ww[0]);
			int subB = Math.min(uu[1], vv[0]);
			int subC = Math.min(vv[1], ww[1]);
			subSizes[r] = new int[] { subA, subB, subC };
			totalRank += sota.getRank(subA, subB, subC);
		}
		return new Result(totalRank, new int[][] { allocA.clone(), allocB.clone(), allocC.clone() }, subSizes);
	}

	private static int[] applyPeel(int[] alloc, int[] peel) {
		if (peel == null) return alloc;
		if (peel.length != alloc.length) {
			throw new IllegalArgumentException(
					"peel length " + peel.length + " must equal alloc length " + alloc.length);
		}
		int[] eff = new int[alloc.length];
		for (int i = 0; i < alloc.length; i++) {
			if (peel[i] < 0 || peel[i] > alloc[i]) {
				throw new IllegalArgumentException(
						"peel[" + i + "]=" + peel[i] + " out of range [0, alloc=" + alloc[i] + "]");
			}
			eff[i] = alloc[i] - peel[i];
		}
		return eff;
	}

	public static Result recombine(int targetA, int targetB, int targetC,
			NonCubicBilinearAlgorithm base, SotaResolver sota) {
		int baseA = base.n, baseB = base.m, baseC = base.p, baseRank = base.r;

		// Reshape U/V/W from [rows·cols][r] to [rows][cols][r] (lazy: we keep them flat).
		// U: [baseA·baseB][r];  V: [baseB·baseC][r];  W: [baseA·baseC][r]
		// (No transpose needed: our W is already row-major C-flatten, matching the
		// post-transpose form in the original Python.)

		List<int[]> fillingsA = blockFillings(baseA, targetA);
		List<int[]> fillingsB = blockFillings(baseB, targetB);
		List<int[]> fillingsC = blockFillings(baseC, targetC);

		long bestRank = Long.MAX_VALUE;
		int[][] bestAlloc = null;
		int[][] bestSubSizes = null;

		for (int[] allocA : fillingsA) {
			for (int[] allocB : fillingsB) {
				for (int[] allocC : fillingsC) {
					long totalRank = 0;
					int[][] subSizes = new int[baseRank][3];

					for (int r = 0; r < baseRank; r++) {
						int[] uu = processAdditions(base.denseU(), r, baseA, baseB, allocA, allocB);
						int[] vv = processAdditions(base.denseV(), r, baseB, baseC, allocB, allocC);
						int[] ww = processAdditions(base.denseW(), r, baseA, baseC, allocA, allocC);

						// The sub-problem only needs the part of the product that
						// reaches the non-zero output region — hence the mins.
						int subA = Math.min(uu[0], ww[0]);
						int subB = Math.min(uu[1], vv[0]);
						int subC = Math.min(vv[1], ww[1]);

						subSizes[r] = new int[] { subA, subB, subC };
						totalRank += sota.getRank(subA, subB, subC);
						if (totalRank >= bestRank) break; // early prune
					}

					if (totalRank < bestRank) {
						bestRank = totalRank;
						bestAlloc = new int[][] {
								allocA.clone(), allocB.clone(), allocC.clone()
						};
						bestSubSizes = subSizes;
					}
				}
			}
		}

		return new Result(bestRank, bestAlloc, bestSubSizes);
	}

	/**
	 * Enumerate all ways to distribute {@code budget} units across {@code numBlocks}
	 * non-negative integer blocks. Count is {@code C(budget + numBlocks - 1, numBlocks - 1)}.
	 */
	public static List<int[]> blockFillings(int numBlocks, int budget) {
		List<int[]> out = new ArrayList<>();
		int[] current = new int[numBlocks];
		fillRec(numBlocks, budget, current, 0, out);
		return out;
	}

	private static void fillRec(int numBlocks, int budget, int[] current, int idx, List<int[]> out) {
		if (idx == numBlocks - 1) {
			current[idx] = budget;
			out.add(current.clone());
			return;
		}
		for (int i = 0; i <= budget; i++) {
			current[idx] = i;
			fillRec(numBlocks, budget - i, current, idx + 1, out);
		}
	}

	/**
	 * For the {@code rank}-th column of {@code factor} (slot U/V/W in
	 * {@code [rows·cols][r]} flatten), find the maximum row-allocation count over
	 * its non-zero entries' row-indices, and similarly for columns. Captures the
	 * effective sub-problem size after the additions implied by the factor vector.
	 */
	public static int[] processAdditions(double[][] factor, int rank, int rows, int cols,
			int[] rowAlloc, int[] colAlloc) {
		int maxRows = 0, maxCols = 0;
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				if (factor[i * cols + j][rank] != 0.0) {
					if (rowAlloc[i] > maxRows) maxRows = rowAlloc[i];
					if (colAlloc[j] > maxCols) maxCols = colAlloc[j];
				}
			}
		}
		return new int[] { maxRows, maxCols };
	}

	// ───────────────────────────────────────────────────────────────────────────
	// Constructive recombination — actually builds the composed algorithm
	// ───────────────────────────────────────────────────────────────────────────

	/** Look up the best known sub-algorithm for a given format. */
	@FunctionalInterface
	public interface AlgorithmLookup {
		Optional<NonCubicBilinearAlgorithm> find(int n, int m, int p);

		/**
		 * Best known rank for {@code ⟨n,m,p⟩}, INCLUDING lineage-only stubs that
		 * {@link #find} skips (it returns no matrices for maxDim&gt;16 stubs). The
		 * default derives it from {@code find} (so a bare-lambda lookup stays
		 * stub-blind, as before); {@link eu.solven.matmul.catalog.FieldAwareLookup}
		 * overrides it to read stub ranks straight from its filename index — which is
		 * what lets sota value stub-based decompositions correctly. Sentinel
		 * {@code Recombination.SotaResolver.UNKNOWN_RANK} when unknown (matches FieldAwareLookup).
		 */
		default int findRank(int n, int m, int p) {
			return find(n, m, p).map(a -> a.r).orElse(Recombination.SotaResolver.UNKNOWN_RANK);
		}
	}

	/**
	 * Constructive sibling of {@link #recombine}: returns an actual
	 * {@link NonCubicBilinearAlgorithm} for the target format, not just a rank
	 * bound. The search picks the best allocation under {@code lookup}; each
	 * base multiplication is then replaced by the corresponding sub-algorithm
	 * (which must exist in {@code lookup}). For uniform allocations the result
	 * is identical to {@link Compose#kroneckerGeneral} with the same factors;
	 * for asymmetric allocations the result captures the AlphaTensor-style
	 * recombination construction.
	 *
	 * <p>Constraint: when {@code recombine}'s {@code min()} reductions produce
	 * a sub-size strictly smaller than the available block, the lookup must
	 * still return an algorithm for exactly that sub-size. (No padding /
	 * truncation is applied — the construction faithfully embeds the smaller
	 * sub-algorithm into the larger block.)</p>
	 *
	 * @throws IllegalStateException if any sub-algorithm needed by the
	 *                               selected allocation is missing from
	 *                               {@code lookup}.
	 */
	public static NonCubicBilinearAlgorithm construct(int targetA, int targetB, int targetC,
			NonCubicBilinearAlgorithm base, AlgorithmLookup lookup) {

		SotaResolver sota = (a, b, c) -> {
			if (a == 0 || b == 0 || c == 0) return 0;
			return lookup.find(a, b, c).map(alg -> alg.r).orElse(Recombination.SotaResolver.UNKNOWN_RANK);
		};
		Result rec = recombine(targetA, targetB, targetC, base, sota);
		return constructFromResult(targetA, targetB, targetC, base, lookup, rec);
	}

	/**
	 * Constructive sibling of {@link #recombineWithAllocation}: builds the
	 * algorithm for the given fixed allocation, no search. Useful to force a
	 * specific Sedoglavic-style split (e.g. {@code [u, v]} on Strassen base
	 * for {@code ⟨u+v,u+v,u+v⟩}) when the default search would collapse to a
	 * degenerate allocation by looking up a strong direct catalog entry.
	 */
	public static NonCubicBilinearAlgorithm constructWithAllocation(
			NonCubicBilinearAlgorithm base, AlgorithmLookup lookup,
			int[] allocA, int[] allocB, int[] allocC) {
		SotaResolver sota = (a, b, c) -> {
			if (a == 0 || b == 0 || c == 0) return 0;
			return lookup.find(a, b, c).map(alg -> alg.r).orElse(Recombination.SotaResolver.UNKNOWN_RANK);
		};
		return constructWithAllocation(base, lookup, sota, allocA, allocB, allocC);
	}

	/**
	 * Variant of {@link #constructWithAllocation(NonCubicBilinearAlgorithm,
	 * AlgorithmLookup, int[], int[], int[])} that uses a CALLER-PROVIDED
	 * {@link SotaResolver} for the rank search, while still using
	 * {@code lookup} for the actual scheme materialisation.
	 *
	 * <p>Use this to inject a formula-aware resolver (e.g.
	 * {@link CitedBound}) — the search will discover allocations
	 * that take advantage of closed-form bounds (Pan TA, Hopcroft-Kerr)
	 * even when the catalog has no scheme for the leaf, but the
	 * materialisation will still fall back to whatever scheme {@code
	 * lookup} returns. If {@code sota} predicts a rank lower than the
	 * scheme {@code lookup} can materialise, the materialised scheme
	 * will be larger than the predicted rank — but it is the best we
	 * can construct with the current catalog.</p>
	 */
	public static NonCubicBilinearAlgorithm constructWithAllocation(
			NonCubicBilinearAlgorithm base, AlgorithmLookup lookup, SotaResolver sota,
			int[] allocA, int[] allocB, int[] allocC) {
		Result rec = recombineWithAllocation(base, sota, allocA, allocB, allocC);
		int targetA = sum(allocA), targetB = sum(allocB), targetC = sum(allocC);
		return constructFromResult(targetA, targetB, targetC, base, lookup, rec);
	}

	private static int sum(int[] a) {
		int s = 0;
		for (int x : a) s += x;
		return s;
	}

	private static NonCubicBilinearAlgorithm constructFromResult(int targetA, int targetB, int targetC,
			NonCubicBilinearAlgorithm base, AlgorithmLookup lookup, Result rec) {

		int[] allocA = rec.allocations[0];
		int[] allocB = rec.allocations[1];
		int[] allocC = rec.allocations[2];
		int[] cumA = cumulative(allocA);
		int[] cumB = cumulative(allocB);
		int[] cumC = cumulative(allocC);

		int baseR = base.r;
		int[] subStarts = new int[baseR + 1];
		NonCubicBilinearAlgorithm[] subs = new NonCubicBilinearAlgorithm[baseR];
		for (int k = 0; k < baseR; k++) {
			int[] sz = rec.smallMatrixSizes[k];
			if (sz[0] == 0 || sz[1] == 0 || sz[2] == 0) {
				subs[k] = null;
				subStarts[k + 1] = subStarts[k];
			} else {
				NonCubicBilinearAlgorithm sub = lookup.find(sz[0], sz[1], sz[2]).orElse(null);
				if (sub == null) {
					// A width-1 sub-block ⟨1,m,p⟩/⟨n,1,p⟩/⟨n,m,1⟩ is the trivial naive
					// scheme (rank = n·m·p, optimal — a single row/col has no bilinear
					// savings); it is never catalogued, so build it directly rather than
					// throwing. Mirrors LineageReplayer.resolveLeaf's width-1 handling.
					// Without this, a projection-derived allocation that yields a size-1
					// block produces an unreplayable stub (the ⟨1,7,7⟩ "base mult 0" bug).
					if (sz[0] == 1 || sz[1] == 1 || sz[2] == 1) {
						sub = NonCubicBilinearAlgorithm.naive(sz[0], sz[1], sz[2]);
					} else {
						throw new IllegalStateException(String.format(
								"construct: missing sub-algorithm for ⟨%d,%d,%d⟩ (base mult %d)",
								sz[0], sz[1], sz[2], k));
					}
				}
				subs[k] = sub;
				subStarts[k + 1] = subStarts[k] + sub.r;
			}
		}
		int totalRank = subStarts[baseR];
		if (totalRank == 0) {
			throw new IllegalStateException("construct: empty composed algorithm");
		}

		double[][] U = new double[targetA * targetB][totalRank];
		double[][] V = new double[targetB * targetC][totalRank];
		double[][] W = new double[targetA * targetC][totalRank];

		for (int kBase = 0; kBase < baseR; kBase++) {
			if (subs[kBase] == null) continue;
			NonCubicBilinearAlgorithm sub = subs[kBase];
			int kStart = subStarts[kBase];

			// U: base shape (base.n × base.m), sub shape (sub.n × sub.m), block placement uses (cumA, cumB).
			embedFactor(U, base.denseU(), sub.denseU(), kBase, kStart,
					base.n, base.m, sub.n, sub.m, cumA, cumB, targetB);
			// V: base shape (base.m × base.p), sub shape (sub.m × sub.p), placement (cumB, cumC).
			embedFactor(V, base.denseV(), sub.denseV(), kBase, kStart,
					base.m, base.p, sub.m, sub.p, cumB, cumC, targetC);
			// W: base shape (base.n × base.p), sub shape (sub.n × sub.p), placement (cumA, cumC).
			embedFactor(W, base.denseW(), sub.denseW(), kBase, kStart,
					base.n, base.p, sub.n, sub.p, cumA, cumC, targetC);
		}

		return new NonCubicBilinearAlgorithm(targetA, targetB, targetC, U, V, W);
	}

	private static int[] cumulative(int[] alloc) {
		int[] cum = new int[alloc.length + 1];
		for (int i = 0; i < alloc.length; i++) cum[i + 1] = cum[i] + alloc[i];
		return cum;
	}

	// ===================================================================
	// Generic Pan trilinear-aggregation (TA) over an arbitrary base.
	//
	// For a NAÏVE-GRID base (each base product (i,j,l) is a single-block
	// sub-matmul A_block(i,j)·B_block(j,l)→C_block(i,l), coeff +1 — e.g. the
	// ⟨1,2,2⟩ peel carrier or FMM's ⟨2,3,3⟩ block grid), two products whose
	// sub-shapes are cyclic rotations and whose A/B/C blocks are all distinct
	// (so the products are disjoint) can be FUSED via Pan TA at
	// fusedRank = abc+ab+bc+ca instead of the two leaves' summed rank. This is
	// the generic version of RectangularTrilinearAggregation.buildPeeledViaTa
	// (which only did the symmetric ⟨1,2,2⟩ peel): the combined-space TA block
	// build(n,r,p) is mapped onto the two products' global block positions.
	// ===================================================================

	/** A composed scheme together with the base-product index pairs that were TA-fused
	 *  (each {@code int[]{k1,k2}} oriented so product {@code k2}'s shape is the rot² of
	 *  {@code k1}'s — the orientation {@code build} expects). The pair list is what a
	 *  lineage records so replay re-fuses identically. */
	public record TaFusedConstruction(NonCubicBilinearAlgorithm alg, List<int[]> fusedPairs) {}

	/** {@code true} iff {@code base} is a NAÏVE GRID — every product is a single-block
	 *  sub-matmul (coeff +1), i.e. {@code base.r == n·m·p} and each {@link #singleBlockGrid}
	 *  is non-null. Only naïve grids can carry TA-fusable disjoint cyclic pairs (FMM's
	 *  ⟨1,2,2⟩ peel, ⟨2,3,3⟩ grids, …); a block-combining base like Strassen cannot. */
	public static boolean isNaiveGrid(NonCubicBilinearAlgorithm base) {
		if (base.r != base.n * base.m * base.p) return false;
		for (int k = 0; k < base.r; k++) {
			if (singleBlockGrid(base, k) == null) return false;
		}
		return true;
	}

	/** Grid coords {@code (i,j,l)} of a SINGLE-BLOCK base product {@code k}
	 *  (A-block(i,j)·B-block(j,l)→C-block(i,l)), or {@code null} when product
	 *  {@code k} is not single-block / not coeff +1 (so not TA-fusable). */
	private static int[] singleBlockGrid(NonCubicBilinearAlgorithm base, int k) {
		double[][] u = base.denseU(), v = base.denseV(), w = base.denseW();
		int[] uc = onlyOnePlusOne(u, k, base.m);   // (i,j)
		int[] vc = onlyOnePlusOne(v, k, base.p);   // (j,l)
		int[] wc = onlyOnePlusOne(w, k, base.p);   // (i,l)
		if (uc == null || vc == null || wc == null) return null;
		// Consistency of a naïve matmul product: A-col == B-row (j), A-row == C-row (i),
		// B-col == C-col (l).
		if (uc[1] != vc[0] || uc[0] != wc[0] || vc[1] != wc[1]) return null;
		return new int[] { uc[0], uc[1], vc[1] };  // (i, j, l)
	}

	/** If column {@code k} of factor {@code f} (flattened rows×cols, cols given) has
	 *  EXACTLY one non-zero and it equals +1, return its {@code (row,col)}; else null. */
	private static int[] onlyOnePlusOne(double[][] f, int k, int cols) {
		int found = -1;
		for (int a = 0; a < f.length; a++) {
			double c = f[a][k];
			if (c == 0.0) continue;
			if (c != 1.0 || found != -1) return null;  // ≠+1, or a second non-zero
			found = a;
		}
		if (found == -1) return null;
		return new int[] { found / cols, found % cols };
	}

	/** {@code true} iff {@code s2} is the rot² of {@code s1} ({@code ⟨p,n,r⟩} of
	 *  {@code ⟨n,r,p⟩}) — the orientation {@link
	 *  eu.solven.matmul.papers.pan1978.RectangularTrilinearAggregation#build} fuses. */
	private static boolean isRot2(int[] s1, int[] s2) {
		return s2[0] == s1[2] && s2[1] == s1[0] && s2[2] == s1[1];
	}

	/**
	 * Greedy max-savings matching of disjoint cyclic-rotation SINGLE-BLOCK product
	 * pairs (the {@link PairedSubProducts#applyPairing} matching, but returning the
	 * matched index pairs and requiring single-block + block-disjointness so the
	 * fusion is constructible). Each returned pair {@code {k1,k2}} is oriented so
	 * product {@code k2}'s sub-shape is the rot² of {@code k1}'s.
	 */
	private static List<int[]> findTaPairs(int[][] subShapes, int[][] grids, long[] subRanks) {
		int n = subShapes.length;
		List<long[]> cand = new ArrayList<>();   // {savings, k1, k2}
		for (int i = 0; i < n; i++) {
			if (grids[i] == null) continue;
			int[] si = subShapes[i];
			for (int j = i + 1; j < n; j++) {
				if (grids[j] == null) continue;
				int[] sj = subShapes[j];
				// orient so (k1=a, k2=b) has sb = rot²(sa)
				int k1, k2;
				if (isRot2(si, sj)) { k1 = i; k2 = j; }
				else if (isRot2(sj, si)) { k1 = j; k2 = i; }
				else continue;  // not a (non-trivial) cyclic pair build can fuse
				// blocks must be pairwise distinct (disjoint products)
				int[] g1 = grids[k1], g2 = grids[k2];
				boolean distinctA = !(g1[0] == g2[0] && g1[1] == g2[1]);
				boolean distinctB = !(g1[1] == g2[1] && g1[2] == g2[2]);
				boolean distinctC = !(g1[0] == g2[0] && g1[2] == g2[2]);
				if (!(distinctA && distinctB && distinctC)) continue;
				long fused = PairedSubProducts.pairCost(si[0], si[1], si[2]);
				long savings = subRanks[i] + subRanks[j] - fused;
				if (savings <= 0) continue;
				cand.add(new long[] { savings, k1, k2 });
			}
		}
		cand.sort((x, y) -> Long.compare(y[0], x[0]));
		boolean[] used = new boolean[n];
		List<int[]> pairs = new ArrayList<>();
		for (long[] c : cand) {
			int k1 = (int) c[1], k2 = (int) c[2];
			if (used[k1] || used[k2]) continue;
			used[k1] = used[k2] = true;
			pairs.add(new int[] { k1, k2 });
		}
		return pairs;
	}

	/**
	 * Constructive recombination WITH generic Pan-TA fusion of disjoint
	 * cyclic-rotation single-block product pairs. Returns the composed scheme (each
	 * fused pair contributes {@code fusedRank} columns instead of the two leaves'
	 * summed rank) plus the fused-pair list (for the lineage). Falls back to exactly
	 * {@link #constructFromResult} when no pair fuses.
	 */
	public static TaFusedConstruction constructWithTaFusion(
			NonCubicBilinearAlgorithm base, AlgorithmLookup lookup, SotaResolver sota,
			int[] allocA, int[] allocB, int[] allocC) {
		// find()-based sub resolver (skips stubs); callers with stub leaves (cube/corner
		// of a peel) must use the SubResolver overload with a replay-capable resolver.
		return constructWithTaFusion(base, (sz) ->
				lookup.find(sz[0], sz[1], sz[2]).orElse(null), sota, allocA, allocB, allocC);
	}

	/** Resolves a sub-shape {@code ⟨n,m,p⟩} (oriented) to an actual scheme — possibly via
	 *  REPLAY of a lineage-only stub. Returns {@code null} if unresolvable. */
	@FunctionalInterface
	public interface SubResolver {
		NonCubicBilinearAlgorithm resolve(int[] shape);
	}

	/**
	 * Generic TA-fused construction with a caller-supplied {@link SubResolver} for the
	 * UNPAIRED leaves — so a peel whose cube/corner are lineage-only STUBS (skipped by
	 * {@code find()}) can still be built (the materialiser passes a replay-capable
	 * resolver). The {@code sota} resolver (rank-only) drives the pairing decision, so
	 * the fused pair's leaves are never materialised at all (TA replaces them).
	 */
	public static TaFusedConstruction constructWithTaFusion(
			NonCubicBilinearAlgorithm base, SubResolver resolveSub, SotaResolver sota,
			int[] allocA, int[] allocB, int[] allocC) {
		Result rec = recombineWithAllocation(base, sota, allocA, allocB, allocC);
		int targetA = sum(allocA), targetB = sum(allocB), targetC = sum(allocC);
		int[] cumA = cumulative(allocA), cumB = cumulative(allocB), cumC = cumulative(allocC);
		int baseR = base.r;

		// Per-product: single-block grid coords + rank (rank-only via sota — no scheme yet).
		int[][] grids = new int[baseR][];
		long[] subRanks = new long[baseR];
		for (int k = 0; k < baseR; k++) {
			int[] sz = rec.smallMatrixSizes[k];
			if (sz[0] == 0 || sz[1] == 0 || sz[2] == 0) continue;
			subRanks[k] = sota.getRank(sz[0], sz[1], sz[2]);
			grids[k] = singleBlockGrid(base, k);
		}

		List<int[]> pairs = findTaPairs(rec.smallMatrixSizes, grids, subRanks);
		boolean[] inPair = new boolean[baseR];
		for (int[] pr : pairs) { inPair[pr[0]] = true; inPair[pr[1]] = true; }

		// Resolve the UNPAIRED leaves to actual schemes (fused leaves need none).
		NonCubicBilinearAlgorithm[] subs = new NonCubicBilinearAlgorithm[baseR];
		long total = 0;
		for (int k = 0; k < baseR; k++) {
			int[] sz = rec.smallMatrixSizes[k];
			if (sz[0] == 0 || sz[1] == 0 || sz[2] == 0 || inPair[k]) continue;
			NonCubicBilinearAlgorithm sub = (sz[0] == 1 || sz[1] == 1 || sz[2] == 1)
					? NonCubicBilinearAlgorithm.naive(sz[0], sz[1], sz[2])
					: resolveSub.resolve(sz);
			if (sub == null) {
				throw new IllegalStateException(String.format(
						"constructWithTaFusion: missing sub-algorithm for ⟨%d,%d,%d⟩ (base mult %d)",
						sz[0], sz[1], sz[2], k));
			}
			subs[k] = sub;
			total += sub.r;
		}
		for (int[] pr : pairs) {
			int[] s = rec.smallMatrixSizes[pr[0]];
			total += PairedSubProducts.pairCost(s[0], s[1], s[2]);
		}
		if (total == 0 || total > Integer.MAX_VALUE) {
			throw new IllegalStateException("constructWithTaFusion: bad rank " + total);
		}
		int totalRank = (int) total;

		double[][] U = new double[targetA * targetB][totalRank];
		double[][] V = new double[targetB * targetC][totalRank];
		double[][] W = new double[targetA * targetC][totalRank];
		int col = 0;

		// Fused pairs first.
		for (int[] pr : pairs) {
			col = embedTaPair(U, V, W, col, pr[0], pr[1], rec.smallMatrixSizes,
					grids, cumA, cumB, cumC, targetB, targetC);
		}
		// Then the unpaired products, via the normal embed.
		for (int kBase = 0; kBase < baseR; kBase++) {
			if (subs[kBase] == null) continue;
			NonCubicBilinearAlgorithm sub = subs[kBase];
			embedFactor(U, base.denseU(), sub.denseU(), kBase, col,
					base.n, base.m, sub.n, sub.m, cumA, cumB, targetB);
			embedFactor(V, base.denseV(), sub.denseV(), kBase, col,
					base.m, base.p, sub.m, sub.p, cumB, cumC, targetC);
			embedFactor(W, base.denseW(), sub.denseW(), kBase, col,
					base.n, base.p, sub.n, sub.p, cumA, cumC, targetC);
			col += sub.r;
		}
		return new TaFusedConstruction(
				new NonCubicBilinearAlgorithm(targetA, targetB, targetC, U, V, W), pairs);
	}

	/** Embed the fused TA block {@code build(n,r,p)} for the (rot²-oriented) product
	 *  pair {@code (k1,k2)} into the global factors starting at column {@code col};
	 *  returns the next free column. {@code k1}'s sub-shape is {@code ⟨n,r,p⟩},
	 *  {@code k2}'s is {@code ⟨p,n,r⟩}. */
	private static int embedTaPair(double[][] U, double[][] V, double[][] W, int col,
			int k1, int k2, int[][] subShapes, int[][] grids,
			int[] cumA, int[] cumB, int[] cumC, int targetB, int targetC) {
		int[] s1 = subShapes[k1];
		int n = s1[0], r = s1[1], p = s1[2];
		int[] g1 = grids[k1], g2 = grids[k2];   // (i,j,l) each
		int i1 = g1[0], j1 = g1[1], l1 = g1[2];
		int i2 = g2[0], j2 = g2[1], l2 = g2[2];
		var ta = eu.solven.matmul.papers.pan1978.RectangularTrilinearAggregation.build(n, r, p);
		int nr = n * r, rp = r * p, np = n * p;
		for (var t : ta.terms()) {
			int[] tu = t.u(), tv = t.v(), tw = t.w();
			for (int idx = 0; idx < tu.length; idx++) {
				if (tu[idx] == 0) continue;
				int gRow, gCol;
				if (idx < nr) {           // A1(n×r) → P1 A-block (i1,j1)
					gRow = cumA[i1] + idx / r; gCol = cumB[j1] + idx % r;
				} else {                  // A2(p×n) → P2 A-block (i2,j2)
					int d = idx - nr; gRow = cumA[i2] + d / n; gCol = cumB[j2] + d % n;
				}
				U[gRow * targetB + gCol][col] += tu[idx];
			}
			for (int idx = 0; idx < tv.length; idx++) {
				if (tv[idx] == 0) continue;
				int gRow, gCol;
				if (idx < rp) {           // B1(r×p) → P1 B-block (j1,l1)
					gRow = cumB[j1] + idx / p; gCol = cumC[l1] + idx % p;
				} else {                  // B2(n×r) → P2 B-block (j2,l2)
					int d = idx - rp; gRow = cumB[j2] + d / r; gCol = cumC[l2] + d % r;
				}
				V[gRow * targetC + gCol][col] += tv[idx];
			}
			for (int idx = 0; idx < tw.length; idx++) {
				if (tw[idx] == 0) continue;
				int gRow, gCol;
				if (idx < np) {           // C1(n×p) → P1 C-block (i1,l1)
					gRow = cumA[i1] + idx / p; gCol = cumC[l1] + idx % p;
				} else {                  // C2(p×r) → P2 C-block (i2,l2)
					int d = idx - np; gRow = cumA[i2] + d / r; gCol = cumC[l2] + d % r;
				}
				W[gRow * targetC + gCol][col] += tw[idx];
			}
			col++;
		}
		return col;
	}

	/**
	 * Place {@code base[aBase][kBase] · sub[aSub][kSub]} into
	 * {@code dst[globalRow * targetCols + globalCol][kStart + kSub]} for every
	 * non-zero base/sub pair, where {@code (globalRow, globalCol)} is the block
	 * offset given by the allocation {@code cumulatives}.
	 */
	private static void embedFactor(double[][] dst, double[][] base, double[][] sub,
			int kBase, int kStart,
			int baseRows, int baseCols, int subRows, int subCols,
			int[] cumRows, int[] cumCols, int targetCols) {
		int subRank = sub[0].length;
		int baseDim = baseRows * baseCols;
		int subDim = subRows * subCols;
		for (int aBase = 0; aBase < baseDim; aBase++) {
			double c = base[aBase][kBase];
			if (c == 0.0) continue;
			int iBase = aBase / baseCols;
			int jBase = aBase % baseCols;
			int rowOff = cumRows[iBase];
			int colOff = cumCols[jBase];
			// Block's actual size — may be smaller than sub.{rows,cols} for non-uniform
			// allocations (e.g. block (1,1) is 3×3 but the sub-algorithm shape is 4×4
			// when Strassen's M1 stitches A_uu + A_vv via padding). Out-of-range
			// (iSub, jSub) corresponds to a padded zero in the input — skip it.
			int blockRowSize = cumRows[iBase + 1] - cumRows[iBase];
			int blockColSize = cumCols[jBase + 1] - cumCols[jBase];
			for (int aSub = 0; aSub < subDim; aSub++) {
				int iSub = aSub / subCols;
				int jSub = aSub % subCols;
				if (iSub >= blockRowSize || jSub >= blockColSize) continue;
				int aGlob = (rowOff + iSub) * targetCols + (colOff + jSub);
				double[] srcRow = sub[aSub];
				double[] dstRow = dst[aGlob];
				for (int kSub = 0; kSub < subRank; kSub++) {
					double s = srcRow[kSub];
					if (s == 0.0) continue;
					dstRow[kStart + kSub] += c * s;
				}
			}
		}
	}
}
