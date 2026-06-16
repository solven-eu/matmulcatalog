package eu.solven.matmul.papers.hopcroftkerr1971;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import eu.solven.matmul.NonCubicBilinearAlgorithm;

/**
 * Constructive implementation of Hopcroft-Kerr 1971's algorithm for
 * {@code ⟨p, 2, n⟩} matrix multiplication of rank
 * {@code ⌈(3pn + max(p,n))/2⌉}.
 *
 * <p>The paper presents a <em>meta-algorithm</em> with two parameter
 * choices: (i) linear functionals (Lemma 1) used to extend {@code A}
 * for asymmetric cases, and (ii) a method sequence (Lemma 3) assigning
 * each diagonal {@code y_{ii}} to one of three Strassen-style
 * decompositions. This class fixes <em>specific</em> choices for both
 * and emits an explicit {@link NonCubicBilinearAlgorithm}.</p>
 *
 * <p>Reference: Hopcroft &amp; Kerr 1971, <em>SIAM J. Applied Math
 * 20(1):30–36</em>. PDF archived locally at
 * {@code references/papers/hopcroft_kerr_1971_2bc_2n2.pdf}. The 1969
 * Cornell Tech Report TR 69-44 contains the full 34-page proof; the
 * SIAM article is a 7-page condensation.</p>
 *
 * <p>Implementation roadmap (this class):</p>
 * <table>
 *   <caption>HK2bc implementation status</caption>
 *   <tr><th>Component</th><th>Paper section</th><th>Status</th></tr>
 *   <tr><td>7 product templates (A..G)</td><td>p3</td><td>✓ done</td></tr>
 *   <tr><td>3 diagonal methods</td><td>p9</td><td>✓ done</td></tr>
 *   <tr><td>Lemma 2: 3 joint-pair cases (1↔2, 1↔3, 2↔3)</td><td>p7</td><td>✓ done</td></tr>
 *   <tr><td>Linear-form bookkeeping + product sharing</td><td>—</td><td>✓ done</td></tr>
 *   <tr><td>Square odd-n, n=3</td><td>Case 1, p11</td><td>✓ done (test verifies)</td></tr>
 *   <tr><td>Square odd-n, n ≥ 5 (same-method pair fallback)</td><td>p10</td><td>TODO — explicit page-10 case 1+1+bridge-2 only; cases 2+2 and 3+3 need analogous derivation</td></tr>
 *   <tr><td>Square even-n: Lemma 3 sequence builder</td><td>p7-8</td><td>TODO</td></tr>
 *   <tr><td>Square even-n: Steps 1+2+3</td><td>p11-12</td><td>TODO</td></tr>
 *   <tr><td>Asymmetric p &lt; n ≤ 2p: Lemma 1 augmentation</td><td>p5-6, p9</td><td>TODO</td></tr>
 *   <tr><td>Decomposition for n &gt; 2p</td><td>p8 Theorem 1 proof</td><td>TODO (trivial once asymmetric works)</td></tr>
 * </table>
 *
 * <p>To close the catalog's 3 remaining FMM-better gaps
 * ({@code ⟨2,10,15⟩}, {@code ⟨2,10,16⟩}, {@code ⟨2,12,16⟩}) we need
 * the asymmetric row (Lemma 1 augmentation). The square-even-n case
 * doesn't add catalog value (FMM Maple files cover those) but is a
 * useful stepping stone.</p>
 */
public final class HopcroftKerr2bc {

	private HopcroftKerr2bc() {}

	/** Closed-form rank: {@code ⌈(3pn + max(p,n))/2⌉}. */
	public static int rank(int p, int n) {
		long num = 3L * p * n + Math.max(p, n);
		return (int) ((num + 1) / 2);
	}

	/**
	 * Build the {@code ⟨n, 2, n⟩} scheme for <strong>odd</strong>
	 * {@code n}. Follows paper §3 Case 1 (page 11): diagonal methods
	 * alternate between {@code A+B} (odd rows except last), {@code -C+D}
	 * (even rows), and {@code E+F} (last row if odd-indexed).
	 */
	public static NonCubicBilinearAlgorithm buildSquareOdd(int n) {
		if (n < 1 || (n & 1) == 0) {
			throw new IllegalArgumentException("buildSquareOdd: n must be odd, got " + n);
		}
		Builder b = new Builder(n, 2, n);

		// Diagonal method assignment per page 11 Case 1.
		// Indexing: rows i ∈ [1..n]. Use 1-based to match the paper.
		// method[i] ∈ {1, 2, 3} per:
		//   i odd, i ≠ n  → 1 (A + B)
		//   i even        → 2 (-C + D)
		//   i odd, i = n  → 3 (E + F)
		int[] method = new int[n + 1];
		for (int i = 1; i <= n; i++) {
			if ((i & 1) == 1 && i != n) method[i] = 1;
			else if ((i & 1) == 0)       method[i] = 2;
			else                         method[i] = 3;
		}

		// Step 1: compute diagonals y_{ii}.
		for (int i = 1; i <= n; i++) {
			emitDiagonal(b, i, method[i]);
		}

		// Step 2: pair processing in DEPENDENCY ORDER (by increasing
		// j - i distance):
		//   - Adjacent pairs (distance 1) emit E and F products that
		//     same-method pairs reuse.
		//   - A same-method pair (i, j) at distance d also needs the
		//     bridge pair (i+1, j) at distance d-1 to be already
		//     computed (page-10 subtraction step). Processing pairs
		//     by ascending distance guarantees this.
		for (int dist = 1; dist < n; dist++) {
			for (int i = 1; i + dist <= n; i++) {
				int j = i + dist;
				if (dist == 1 || method[i] != method[j]) {
					emitOffDiagonalPair(b, i, j, method[i], method[j]);
				} else {
					emitSameMethodPair(b, i, j, method[i], method[i + 1]);
				}
			}
		}

		return b.build();
	}

	/**
	 * Build the {@code ⟨n, 2, n⟩} scheme for <strong>even</strong>
	 * {@code n}. Follows paper §3 Case 2 (page 11) with the simplest
	 * Lemma-3-conformant sequence: strict alternation
	 * {@code (1, 2, 1, 2, …, 1, 2)}. For {@code n = 2k+2}, every
	 * same-method pair is either {@code (1,1)} with a method-2 bridge
	 * or {@code (2,2)} with a method-1 bridge — both implemented.
	 *
	 * <p>The paper's Step 2 boundary pairs (at cyclic distance
	 * {@code k+1 = n/2}) end up with differing methods under this
	 * coloring, so they're handled by the standard Lemma 2 emission
	 * without needing a separate boundary case.</p>
	 *
	 * <p>The Lemma 3 complications (variable {@code ℓ} / {@code m}
	 * boundary regions) only matter for {@code n > 2k+2}, i.e. when
	 * the asymmetric augmentation in task #48 inflates the target
	 * row-count past the natural square. For pure square
	 * {@code ⟨n,2,n⟩} with {@code n = 2k+2}, this implementation is
	 * complete.</p>
	 */
	public static NonCubicBilinearAlgorithm buildSquareEven(int n) {
		if (n < 2 || (n & 1) != 0) {
			throw new IllegalArgumentException("buildSquareEven: n must be even, got " + n);
		}
		Builder b = new Builder(n, 2, n);
		int[] method = new int[n + 1];
		for (int i = 1; i <= n; i++) {
			method[i] = ((i & 1) == 1) ? 1 : 2;
		}

		for (int i = 1; i <= n; i++) {
			emitDiagonal(b, i, method[i]);
		}
		// Same dependency-ordered pair processing as the odd case — the
		// distance-by-distance order ensures every same-method pair's
		// bridge has been computed before subtraction.
		for (int dist = 1; dist < n; dist++) {
			for (int i = 1; i + dist <= n; i++) {
				int j = i + dist;
				if (dist == 1 || method[i] != method[j]) {
					emitOffDiagonalPair(b, i, j, method[i], method[j]);
				} else {
					emitSameMethodPair(b, i, j, method[i], method[i + 1]);
				}
			}
		}
		return b.build();
	}

	/** Convenience dispatcher: {@link #buildSquareOdd} or {@link #buildSquareEven}. */
	public static NonCubicBilinearAlgorithm buildSquare(int n) {
		return (n & 1) == 1 ? buildSquareOdd(n) : buildSquareEven(n);
	}

	/**
	 * Band-restricted variant of {@link #buildSquareOdd} for the
	 * asymmetric HK Case 1 (paper §4 Case 1, p9-11). Computes only the
	 * cells {@code (i', j) ∈ S} where the cyclic distance
	 * {@code |i' - j|_cyc ≤ k = (p-1)/2}; out-of-band cells have
	 * {@code W = 0}. Used as the internal scheme inside
	 * {@code HopcroftKerr2bcAsymmetric.buildOdd}, which augments
	 * {@code A → Ā = M·A} via Lemma 1 and back-substitutes {@code AX}
	 * via per-column Vandermonde inverse.
	 *
	 * <p>Method assignment per Lemma 3 simplest sequence:
	 * {@code (1, 2, 1, 2, ..., 1, 2, 3)} for n odd, with cyclic adjacency
	 * preserved (row n→1 has methods 3→1, different).</p>
	 *
	 * <p>Internal cost: {@code n} diagonals (2 mults each) +
	 * {@code n·k} cyclic-band pairs (3 mults each) =
	 * {@code 2n + 3nk = (3pn + n)/2} for {@code p = 2k+1}.</p>
	 *
	 * @param n must be ≥ p and the working n-dimension
	 * @param k band half-width, ≥ 1; must satisfy {@code 2k+1 ≤ n} (band fits)
	 */
	public static NonCubicBilinearAlgorithm buildOddBanded(int n, int k) {
		if (n < 3) throw new IllegalArgumentException("n must be ≥ 3");
		if (k < 1) throw new IllegalArgumentException("k must be ≥ 1");
		if (2 * k + 1 > n) throw new IllegalArgumentException("band 2k+1 must fit in n");
		Builder b = new Builder(n, 2, n);

		int[] method = new int[n + 1];
		if ((n & 1) == 1) {
			// Odd n: (1, 2, 1, 2, …, 1, 2, 3) — last row method 3.
			for (int i = 1; i <= n; i++) {
				if ((i & 1) == 1 && i != n) method[i] = 1;
				else if ((i & 1) == 0)       method[i] = 2;
				else                         method[i] = 3;
			}
		} else {
			// Even n: strict alternation (1, 2, 1, 2, ..., 1, 2).
			for (int i = 1; i <= n; i++) {
				method[i] = ((i & 1) == 1) ? 1 : 2;
			}
		}

		for (int i = 1; i <= n; i++) {
			emitDiagonal(b, i, method[i]);
		}

		// Cyclic-band pair processing. For each cyclic distance d ∈ [1, k]
		// and each i ∈ [1, n], the pair partner j = ((i-1+d) % n) + 1.
		// Each unordered pair {i, j} at cyclic distance d ≤ k appears
		// exactly once because k < n/2 (no diametric overlap).
		// Bridge SELECTION (task #7) happens inside emitBandPair: same-method
		// pairs use an arc-interior position with the derivable OPPOSITE method
		// (1↔2), so only the proven (1,1,b2)/(2,2,b1) emitters run and the
		// impossible (2,2,bridge-3) configuration never arises. The historical
		// hardwired "cyclic next row" could land on the lone method-3 position.
		for (int dist = 1; dist <= k; dist++) {
			for (int i = 1; i <= n; i++) {
				int j = ((i - 1 + dist) % n) + 1;
				emitBandPair(b, method, n, i, j, dist);
			}
		}

		return b.build();
	}

	/**
	 * Result of {@link #buildEvenBanded}: the banded scheme plus, per output
	 * column {@code j0} (0-based), whether its EXTRA cell (the one completing the
	 * even window {@code p = 2k+2}) sits at row {@code j0+(k+1)} (up, window
	 * {@code [j0−k, j0+k+1]}) or row {@code j0−(k+1)} (down, window
	 * {@code [j0−k−1, j0+k]}). The asymmetric back-sub needs this to pick the
	 * matching contiguous Lemma-1 window.
	 */
	/**
	 * A repaired-Step-3 "Z-pair" (HK 1971 p12–13, our (3,3)-free variant): two
	 * leftover columns {@code i2} (extra cell UP at row {@code i1 = i2+k+1}) and
	 * {@code i4} (extra cell DOWN at row {@code i3 = i4−k−1}), completed jointly
	 * by a virtual Lemma-2 cross pair on the virtual rows
	 * {@code α = a_{i2}+a_{i3}} (method 3 — its E/F products are the band
	 * (1,2)-pair (i2,i3) cross products) and {@code β = a_{i1}−a_{i4}}
	 * (method {@code mBeta} ∈ {1,2} — its products are the band (·,3)-pair
	 * (i4,i1) virtuals, i1 being recolored to method 3). All indices 1-based.
	 */
	public record ZPair(int i1, int i2, int i3, int i4, int mBeta) {}

	public record BandedEven(NonCubicBilinearAlgorithm alg, boolean[] extraIsUp,
			java.util.List<ZPair> zPairs) {}

	/**
	 * HK §4 <b>Case 2</b> internal computation (task #7): the odd band of
	 * half-width {@code k} (Step 1) PLUS one extra cell per column at cyclic
	 * distance {@code k+1} (Step 2), giving each column {@code p = 2k+2}
	 * contiguous cells. The budget — {@code ⌈(3pn+n)/2⌉ − n(3k+2) ≈ 3n/2} —
	 * is exactly {@code n/2} full pairs at distance {@code k+1}: a perfect
	 * matching in the circulant graph {@code i ↔ i+(k+1) (mod n)}. The matching
	 * is built per orbit cycle of {@code +(k+1) mod n}; odd cycles leave one
	 * leftover column whose extra cell is computed naively (2 products — Step 3,
	 * costing {@code +(g−1)/2}-ish above formula when the circulant decomposes
	 * into {@code g>1} odd cycles; exact formula whenever every cycle is even,
	 * or {@code g=1} odd cycle and {@code n} odd, which the ceiling absorbs).
	 */
	public static BandedEven buildEvenBanded(int n, int k) {
		if (n < 4) throw new IllegalArgumentException("n must be ≥ 4");
		if (k < 1) throw new IllegalArgumentException("k must be ≥ 1");
		if (2 * k + 2 > n) throw new IllegalArgumentException("even band 2k+2 must fit in n");
		Builder b = new Builder(n, 2, n);

		int[] method = new int[n + 1];
		if ((n & 1) == 1) {
			for (int i = 1; i <= n; i++) {
				if ((i & 1) == 1 && i != n) method[i] = 1;
				else if ((i & 1) == 0)       method[i] = 2;
				else                         method[i] = 3;
			}
		} else {
			for (int i = 1; i <= n; i++) {
				method[i] = ((i & 1) == 1) ? 1 : 2;
			}
		}

		// ── Plan Step 2 (circulant matching) + Step 3 (Z-pairs) BEFORE any
		//    emission: Z-pairs recolor their i1 row to method 3, which must be
		//    known when the diagonals and band pairs are emitted. ──
		int s = k + 1;
		boolean[] extraIsUp = new boolean[n];
		java.util.List<int[]> matchEdges = new java.util.ArrayList<>(); // {i0, j0} 0-based
		java.util.List<java.util.List<Integer>> cycles = new java.util.ArrayList<>();
		{
			boolean[] visited = new boolean[n];
			for (int start = 0; start < n; start++) {
				if (visited[start]) continue;
				java.util.List<Integer> cycle = new java.util.ArrayList<>();
				int v = start;
				while (!visited[v]) {
					visited[v] = true;
					cycle.add(v);
					v = (v + s) % n;
				}
				cycles.add(cycle);
			}
		}
		// Even cycles: full alternate-edge matching, no leftovers. Odd cycles: the
		// leftover VERTEX is free — chosen by the Step-3 placement search below.
		java.util.List<java.util.List<Integer>> oddCycles = new java.util.ArrayList<>();
		for (java.util.List<Integer> cycle : cycles) {
			if ((cycle.size() & 1) == 0) {
				for (int t = 0; t + 1 < cycle.size(); t += 2) {
					matchEdges.add(new int[] { cycle.get(t), cycle.get(t + 1) });
				}
			} else {
				oddCycles.add(cycle);
			}
		}
		Step3Plan plan = planStep3(n, k, method, oddCycles);
		matchEdges.addAll(plan.extraEdges);
		// The plan owns the coloring: arc-phase around the recolored 3-rows (or
		// the unchanged base coloring when there are no odd cycles).
		method = plan.method;

		for (int i = 1; i <= n; i++) {
			emitDiagonal(b, i, method[i]);
		}

		// Step 1: the odd band, distances 1..k (identical to buildOddBanded).
		for (int dist = 1; dist <= k; dist++) {
			for (int i = 1; i <= n; i++) {
				int j = ((i - 1 + dist) % n) + 1;
				emitBandPair(b, method, n, i, j, dist);
			}
		}

		// Step 2: distance-(k+1) full pairs on the matching edges.
		for (int[] e : matchEdges) {
			emitBandPair(b, method, n, e[0] + 1, e[1] + 1, s);
			extraIsUp[e[0]] = true;   // column e0's extra = row e1 (up)
			extraIsUp[e[1]] = false;  // column e1's extra = row e0 (down)
		}

		// Step 3a: Z-pairs — emit the virtual Lemma-2 cross pair; the missing
		// cells' W rows receive the Z-combination here; the asymmetric wrapper
		// adds the rational corrections from complete columns.
		for (ZPair z : plan.zPairs) {
			emitVirtualCrossPair(b, n, z);
			extraIsUp[z.i2() - 1] = true;   // column i2's extra = row i1 (up)
			extraIsUp[z.i4() - 1] = false;  // column i4's extra = row i3 (down)
		}

		// Step 3b: an odd leftover count leaves ONE column (n odd — the ceiling
		// absorbs it): naive DOWN extra cell with 2 plain products.
		for (int j0 : plan.naiveLeftovers) {
			int r1 = Math.floorMod(j0 - s, n) + 1;
			int col = j0 + 1;
			LinearA aR1 = LinearA.aEntry(r1, 1, b.p, b.n);
			LinearA aR2 = LinearA.aEntry(r1, 2, b.p, b.n);
			LinearB x1 = LinearB.xEntry(1, col, b.p, b.n);
			LinearB x2 = LinearB.xEntry(2, col, b.p, b.n);
			int k1 = b.addProduct(LinearA.of(aR1), x1);
			int k2 = b.addProduct(LinearA.of(aR2), x2);
			int out = b.outIdx(r1, col);
			b.addToOutput(out, k1, 1.0);
			b.addToOutput(out, k2, 1.0);
			extraIsUp[j0] = false;
		}

		return new BandedEven(b.build(), extraIsUp, plan.zPairs);
	}

	private record Step3Plan(java.util.List<ZPair> zPairs, java.util.List<int[]> extraEdges,
			java.util.List<Integer> naiveLeftovers, int[] method) {}

	/**
	 * Placement search + ARC-PHASE COLORING for the repaired Step 3.
	 *
	 * <p>Structure search: pick a leftover vertex per odd cycle (free choice),
	 * pair the leftovers into Z-pairs (≤1 naive single), choose role orders,
	 * subject to: pair separation δ ∈ [1,k]; rows {i1,i2,i3,i4} distinct;
	 * i1/i3 complete (non-leftover) columns; all method-3 rows pairwise more
	 * than k apart (no (3,3) band pair — the provably-underivable case).</p>
	 *
	 * <p>Coloring: the recolored i1 rows (plus a movable SPARE 3-row when
	 * {@code #3rows ≢ n (mod 2)}) cut the cycle into arcs; each arc alternates
	 * methods 1,2 with phases chained so {@code m(t−1) ≠ m(t+1)} at every
	 * 3-row t — this removes the dist-2 same-method pair whose only interior
	 * position is the 3-row itself (the old hardwired lone-3-at-n could not).
	 * Consistency of the phase chain around the cycle holds iff
	 * {@code Σ arc lengths = n − #3rows} is even — guaranteed by the spare.</p>
	 *
	 * <p>Validation simulates the emission: Z-pair method constraints
	 * ({@code m(i2) ≠ m(i3)} in {1,2}, {@code m(i4)} in {1,2}) and an
	 * interior opposite-method bridge for every same-method pair the band +
	 * matching will emit. Deterministic seeded retries.</p>
	 */
	private static Step3Plan planStep3(int n, int k, int[] baseMethod,
			java.util.List<java.util.List<Integer>> oddCycles) {
		if (oddCycles.isEmpty()) {
			return new Step3Plan(java.util.List.of(), java.util.List.of(), java.util.List.of(),
					baseMethod);
		}
		long seed = 0x9E3779B97F4A7C15L * ((long) n << 16 ^ (k + 1)) + 1;
		int gAll = oddCycles.size();
		for (int attempt = 0; attempt < 120_000; attempt++) {
			seed = seed * 6364136223846793005L + 1442695040888963407L;
			java.util.Random rnd = new java.util.Random(seed);
			// Graceful degradation: spend the first attempts on the full Z-pair
			// count ⌊g/2⌋ (formula-exact); progressively allow fewer pairs (each
			// drop costs +1 over formula via 2 naive cells). The triangle family
			// n = 3(k+1) with g ≥ 6 provably cannot place ≥3 pairwise->k 3-rows
			// (the three arcs would need to sum past n), so it degrades.
			int maxPairs = Math.max(0, gAll / 2 - attempt / 20_000);
			// 1) Leftover vertex per odd cycle, random pairing order.
			java.util.List<Integer> leftovers = new java.util.ArrayList<>();
			for (java.util.List<Integer> cycle : oddCycles) {
				leftovers.add(cycle.get(rnd.nextInt(cycle.size())));
			}
			java.util.Collections.shuffle(leftovers, rnd);
			java.util.Set<Integer> leftoverSet = new java.util.HashSet<>(leftovers);

			// 2) Structural Z-pairs + 3-row placement.
			java.util.List<int[]> zRows = new java.util.ArrayList<>(); // {i1,i2,i3,i4} 0-based
			java.util.Set<Integer> threeRows = new java.util.HashSet<>();
			java.util.List<Integer> naive = new java.util.ArrayList<>();
			boolean ok = true;
			for (int t = 0; t + 1 < leftovers.size(); t += 2) {
				int c = leftovers.get(t);
				int cp = leftovers.get(t + 1);
				if (zRows.size() >= maxPairs) {
					naive.add(c);
					naive.add(cp);
					continue;
				}
				int[] z = tryZRows(n, k, c, cp, leftoverSet, threeRows);
				if (z == null) z = tryZRows(n, k, cp, c, leftoverSet, threeRows);
				if (z == null) { ok = false; break; }
				zRows.add(z);
				threeRows.add(z[0]);
			}
			if (!ok) continue;
			if ((leftovers.size() & 1) == 1) {
				naive.add(leftovers.get(leftovers.size() - 1));
			}

			// 3) Coloring candidates: arc coloring around the 3-rows (no parity
			//    spare — boundary same-method pairs are bridge-3-servable), plus
			//    the base-parity fallback.
			java.util.List<int[]> colorings = new java.util.ArrayList<>();
			if (!threeRows.isEmpty()) {
				colorings.add(colorArcs(n, threeRows, rnd.nextBoolean()));
			}
			base:
			{
				int[] base = new int[n + 1];
				java.util.Set<Integer> all3 = new java.util.HashSet<>(threeRows);
				if ((n & 1) == 1) {
					// odd n: the lone-3 at position n must respect 3-row spacing.
					for (int t : threeRows) {
						int d = Math.floorMod((n - 1) - t, n);
						if (Math.min(d, n - d) <= k) break base;
					}
					all3.add(n - 1);
				}
				for (int v = 1; v <= n; v++) {
					base[v] = all3.contains(v - 1) ? 3 : (((v & 1) == 1) ? 1 : 2);
				}
				// adjacency check (dist-1 pairs need distinct methods)
				for (int v = 1; v <= n; v++) {
					int w = (v % n) + 1;
					if (base[v] == base[w]) break base;
				}
				colorings.add(base);
			}
			if (colorings.isEmpty()) continue;

			// 5) Matching edges per odd cycle (leftover last, alternate) — these
			//    depend only on the leftover choice, not the coloring.
			java.util.List<int[]> extraEdges = new java.util.ArrayList<>();
			for (java.util.List<Integer> cycle : oddCycles) {
				int lv = -1;
				for (int v : cycle) {
					if (leftoverSet.contains(v)) { lv = v; break; }
				}
				int pos = cycle.indexOf(lv);
				int L = cycle.size();
				for (int t = 0; t + 1 < L; t += 2) {
					extraEdges.add(new int[] {
							cycle.get((pos + 1 + t) % L), cycle.get((pos + 1 + t + 1) % L) });
				}
			}

			// 6) Pick the first coloring that satisfies BOTH the Z-pair method
			//    constraints AND bridge availability for every emitted pair.
			int[] method = null;
			java.util.List<ZPair> zPairs = null;
			for (int[] cand : colorings) {
				boolean zOk = true;
				java.util.List<ZPair> zp = new java.util.ArrayList<>();
				for (int[] z : zRows) {
					int m2 = cand[z[1] + 1], m3 = cand[z[2] + 1], m4 = cand[z[3] + 1];
					if (m2 == 3 || m3 == 3 || m4 == 3 || m2 == m3) { zOk = false; break; }
					zp.add(new ZPair(z[0] + 1, z[1] + 1, z[2] + 1, z[3] + 1, (m4 == 1) ? 2 : 1));
				}
				if (!zOk || !bridgesAvailable(n, k, cand, extraEdges)) continue;
				method = cand;
				zPairs = zp;
				break;
			}
			if (method == null) continue;

			return new Step3Plan(zPairs, extraEdges, naive, method);
		}
		throw new IllegalStateException("planStep3: no valid Z-pair placement found for n=" + n
				+ ", k=" + k + " after 50000 attempts");
	}

	/** Structural Z-pair rows {i1,i2,i3,i4} (0-based) for roles (i2=c, i4=cp); null if invalid. */
	private static int[] tryZRows(int n, int k, int c, int cp,
			java.util.Set<Integer> leftoverSet, java.util.Set<Integer> threeRows) {
		int delta = Math.floorMod(cp - c, n);
		if (delta < 1 || delta > k) return null;
		int i2 = c, i4 = cp;
		int i1 = Math.floorMod(c + k + 1, n);
		int i3 = Math.floorMod(cp - k - 1, n);
		if (java.util.Set.of(i1, i2, i3, i4).size() != 4) return null;
		if (leftoverSet.contains(i1) || leftoverSet.contains(i3)) return null;
		if (threeRows.contains(i1) || threeRows.contains(i3)) return null;
		for (int t : threeRows) {
			int d = Math.floorMod(i1 - t, n);
			if (Math.min(d, n - d) <= k) return null;
		}
		return new int[] { i1, i2, i3, i4 };
	}

	/**
	 * Color the cycle minus the 3-rows in alternating 1,2 arcs. No phase-chain
	 * constraint: a boundary where {@code m(t−1) = m(t+1)} yields a dist-2
	 * same-method pair whose only interior is the 3-row — servable by the
	 * TRUE-reusable (·,·,bridge-3) emitters (task #7), so no parity spare is
	 * needed. Returns the 1-based method array.
	 */
	private static int[] colorArcs(int n, java.util.Set<Integer> threeRows, boolean startPhase) {
		int[] method = new int[n + 1];
		for (int t : threeRows) method[t + 1] = 3;
		java.util.List<Integer> threes = new java.util.ArrayList<>(threeRows);
		java.util.Collections.sort(threes);
		int g3 = threes.size();
		boolean phase = startPhase;
		for (int a = 0; a < g3; a++) {
			int from = threes.get(a) + 1;                       // arc starts after this 3-row
			int to = threes.get((a + 1) % g3);                  // and ends before this one
			int len = Math.floorMod(to - from, n);
			boolean cur = phase;
			for (int off = 0; off < len; off++) {
				int v = Math.floorMod(from + off, n);
				method[v + 1] = cur ? 1 : 2;
				cur = !cur;
			}
		}
		return method;
	}

	/** Simulate the emission's bridge selection for every same-method pair. */
	private static boolean bridgesAvailable(int n, int k, int[] method,
			java.util.List<int[]> extraEdges) {
		for (int dist = 1; dist <= k; dist++) {
			for (int i = 1; i <= n; i++) {
				int j = ((i - 1 + dist) % n) + 1;
				if (!pairServable(n, method, i, j, dist)) return false;
			}
		}
		for (int[] e : extraEdges) {
			if (!pairServable(n, method, e[0] + 1, e[1] + 1, k + 1)) return false;
		}
		return true;
	}

	private static boolean pairServable(int n, int[] method, int i, int j, int dist) {
		int mij = method[i], mji = method[j];
		if (dist == 1) return mij != mji; // direct pair needs distinct methods
		if (mij != mji) return true;
		if (mij == 3) return false; // (3,3) same-method pair: underivable
		int want = (mij == 1) ? 2 : 1;
		for (int e = 1; e < dist; e++) {
			int cand = ((i - 1 + e) % n) + 1;
			// Opposite-method bridge, or method-3 via the true-reusable emitters.
			if (method[cand] == want || method[cand] == 3) return true;
		}
		return false;
	}

	/**
	 * Emit the virtual Lemma-2 cross pair for a Z-pair (repaired HK Step 3):
	 * virtual rows {@code α = a_{i2}+a_{i3}} (method 3; products E(α), F(α) =
	 * the band (1,2)-pair (i2,i3) cross products) and {@code β = a_{i1}−a_{i4}}
	 * (method mBeta; products = the band (i4, i1≡3) pair's virtuals). Three NEW
	 * products compute {@code Z(β,α_x)} and {@code Z(α,β_x)}, written into the
	 * missing cells' W rows:
	 * <pre>
	 *   W[i1,i2] += Z(β,α_x)      (wrapper adds  −y_{i1,i3} + y_{i4,i2} + y_{i4,i3})
	 *   W[i3,i4] += −Z(α,β_x)     (wrapper adds  +y_{i2,i1} − y_{i2,i4} + y_{i3,i1})
	 * </pre>
	 * following {@code Z(β,α_x) = y_{i1i2} + y_{i1i3} − y_{i4i2} − y_{i4i3}} and
	 * {@code Z(α,β_x) = y_{i2i1} − y_{i2i4} + y_{i3i1} − y_{i3i4}}.
	 */
	private static void emitVirtualCrossPair(Builder b, int n, ZPair z) {
		int i1 = z.i1(), i2 = z.i2(), i3 = z.i3(), i4 = z.i4();
		// Virtual row/col linear forms.
		LinearA alpha1 = LinearA.aEntry(i2, 1, b.p, b.n).add(LinearA.aEntry(i3, 1, b.p, b.n));
		LinearA alpha2 = LinearA.aEntry(i2, 2, b.p, b.n).add(LinearA.aEntry(i3, 2, b.p, b.n));
		LinearB alphaX1 = LinearB.xEntry(1, i2, b.p, b.n).add(LinearB.xEntry(1, i3, b.p, b.n));
		LinearB alphaX2 = LinearB.xEntry(2, i2, b.p, b.n).add(LinearB.xEntry(2, i3, b.p, b.n));
		LinearA beta1 = LinearA.aEntry(i1, 1, b.p, b.n).sub(LinearA.aEntry(i4, 1, b.p, b.n));
		LinearA beta2 = LinearA.aEntry(i1, 2, b.p, b.n).sub(LinearA.aEntry(i4, 2, b.p, b.n));
		LinearB betaX1 = LinearB.xEntry(1, i1, b.p, b.n).sub(LinearB.xEntry(1, i4, b.p, b.n));
		LinearB betaX2 = LinearB.xEntry(2, i1, b.p, b.n).sub(LinearB.xEntry(2, i4, b.p, b.n));

		// Reused virtual diagonals. α (method 3): E(α), F(α) from the (1,2) band
		// pair (i2,i3). β (method mBeta): from the (·,3) band pair (i4,i1):
		//   mBeta==2 → C(β), D(β) (pair case (1,3));  mBeta==1 → A(β), B(β) ((2,3)).
		int kE_alpha = b.findProduct(LinearA.of(alpha2), alphaX2);
		int kF_alpha = b.findProduct(LinearA.of(alpha1), alphaX1);

		int outBetaAlpha = b.outIdx(i1, i2); // receives +Z(β, α_x)
		int outAlphaBeta = b.outIdx(i3, i4); // receives −Z(α, β_x)

		// Lemma-2 virtual pair with roles (i := β, method mBeta; j := α, method 3).
		// Mirrors the Java (1,3)/(2,3) branches with a_diff = α − β, x_diff = α_x − β_x.
		LinearA d1 = alpha1.sub(beta1);
		LinearA d2 = alpha2.sub(beta2);
		LinearB dx1 = alphaX1.sub(betaX1);
		LinearB dx2 = alphaX2.sub(betaX2);

		if (z.mBeta() == 1) {
			// case (1,3): y_{βα} = A(β) − D(d) + F(α) − G(β, d, βx, dx)
			//             y_{αβ} = B(β) + C(d) + E(α) + G(β, d, βx, dx)
			int kA_beta = b.findProduct(LinearA.of(beta2), betaX1.add(betaX2));
			int kB_beta = b.findProduct(beta1.sub(beta2), betaX1);
			int kD_d = b.addProduct(LinearA.of(d1), dx1.add(dx2));
			int kC_d = b.addProduct(d1.sub(d2), dx2);
			int kG = b.addProduct(d1.add(beta2), betaX1.sub(dx2));

			b.addToOutput(outBetaAlpha, kA_beta, +1.0);
			b.addToOutput(outBetaAlpha, kD_d,    -1.0);
			b.addToOutput(outBetaAlpha, kF_alpha, +1.0);
			b.addToOutput(outBetaAlpha, kG,      -1.0);

			b.addToOutput(outAlphaBeta, kB_beta, -1.0);
			b.addToOutput(outAlphaBeta, kC_d,    -1.0);
			b.addToOutput(outAlphaBeta, kE_alpha, -1.0);
			b.addToOutput(outAlphaBeta, kG,      -1.0);
		} else {
			// case (2,3): y_{βα} = −A(d) + D(β) + E(α) + G(d, β, dx, βx)
			//             y_{αβ} = −B(d) − C(β) + F(α) − G(d, β, dx, βx)
			int kC_beta = b.findProduct(beta1.sub(beta2), betaX2);
			int kD_beta = b.findProduct(LinearA.of(beta1), betaX1.add(betaX2));
			int kA_d = b.addProduct(LinearA.of(d2), dx1.add(dx2));
			int kB_d = b.addProduct(d1.sub(d2), dx1);
			int kG = b.addProduct(beta1.add(d2), dx1.sub(betaX2));

			b.addToOutput(outBetaAlpha, kA_d,    -1.0);
			b.addToOutput(outBetaAlpha, kD_beta, +1.0);
			b.addToOutput(outBetaAlpha, kE_alpha, +1.0);
			b.addToOutput(outBetaAlpha, kG,      +1.0);

			b.addToOutput(outAlphaBeta, kB_d,    +1.0);
			b.addToOutput(outAlphaBeta, kC_beta, +1.0);
			b.addToOutput(outAlphaBeta, kF_alpha, -1.0);
			b.addToOutput(outAlphaBeta, kG,      +1.0);
		}
	}

	/** Shared band-pair emission: direct Lemma-2 pair when methods differ (or
	 *  dist 1), else same-method with an arc-interior opposite-method bridge. */
	private static void emitBandPair(Builder b, int[] method, int n, int i, int j, int dist) {
		int mij = method[i], mji = method[j];
		if (dist == 1 || mij != mji) {
			emitOffDiagonalPair(b, i, j, mij, mji);
			return;
		}
		// Prefer the opposite {1,2}-method bridge ((1,1,b2)/(2,2,b1) emitters);
		// fall back to a method-3 interior position via the TRUE-reusable
		// (·,·,bridge-3) emitters (task #7 — derived 2026-06-11; the (3,3,·)
		// same-method case alone remains underivable).
		int want = (mij == 1) ? 2 : 1;
		int bridgePos = -1;
		for (int e = 1; e < dist; e++) {
			int cand = ((i - 1 + e) % n) + 1;
			if (method[cand] == want) {
				bridgePos = cand;
				break;
			}
		}
		if (bridgePos < 0) {
			for (int e = 1; e < dist; e++) {
				int cand = ((i - 1 + e) % n) + 1;
				if (method[cand] == 3) {
					bridgePos = cand;
					break;
				}
			}
		}
		if (bridgePos < 0) {
			throw new IllegalStateException(String.format(
					"no method-%d (nor method-3) bridge in the arc interior of pair (%d,%d) dist=%d",
					want, i, j, dist));
		}
		emitSameMethodPair(b, i, j, mij, bridgePos, method[bridgePos]);
	}

	// ─────────────────────────────────────────────────────────────────────
	// Lemma 2 — diagonal + off-diagonal pair emitters.
	// ─────────────────────────────────────────────────────────────────────

	/** Emit the 2 products for diagonal {@code y_{ii}} via method m. */
	private static void emitDiagonal(Builder b, int i, int m) {
		// a_i = row i of A (length 2). Linear forms over A:
		LinearA a_i_1 = LinearA.aEntry(i, 1, b.p, b.n);
		LinearA a_i_2 = LinearA.aEntry(i, 2, b.p, b.n);
		LinearA a_i_sum = a_i_1.add(a_i_2); // a_{i,1} + a_{i,2}
		LinearA a_i_diff = a_i_1.sub(a_i_2); // a_{i,1} - a_{i,2}

		LinearB x_i_1 = LinearB.xEntry(1, i, b.p, b.n); // x_{1,i}
		LinearB x_i_2 = LinearB.xEntry(2, i, b.p, b.n); // x_{2,i}
		LinearB x_i_sum = x_i_1.add(x_i_2);

		int outRow_ii = b.outIdx(i, i);

		switch (m) {
			case 1 -> {
				// A(a_i, x_i) = a_{i,2} * (x_{1,i} + x_{2,i})
				int kA = b.addProduct(LinearA.of(a_i_2), x_i_sum);
				// B(a_i, x_i) = (a_{i,1} - a_{i,2}) * x_{1,i}
				int kB = b.addProduct(a_i_diff, x_i_1);
				b.addToOutput(outRow_ii, kA, 1.0);
				b.addToOutput(outRow_ii, kB, 1.0);
			}
			case 2 -> {
				// C(a_i, x_i) = (a_{i,1} - a_{i,2}) * x_{2,i}
				int kC = b.addProduct(a_i_diff, x_i_2);
				// D(a_i, x_i) = a_{i,1} * (x_{1,i} + x_{2,i})
				int kD = b.addProduct(LinearA.of(a_i_1), x_i_sum);
				b.addToOutput(outRow_ii, kC, -1.0);
				b.addToOutput(outRow_ii, kD, 1.0);
			}
			case 3 -> {
				// E(a_i, x_i) = a_{i,2} * x_{2,i}
				int kE = b.addProduct(LinearA.of(a_i_2), x_i_2);
				// F(a_i, x_i) = a_{i,1} * x_{1,i}
				int kF = b.addProduct(LinearA.of(a_i_1), x_i_1);
				b.addToOutput(outRow_ii, kE, 1.0);
				b.addToOutput(outRow_ii, kF, 1.0);
			}
			default -> throw new IllegalArgumentException("bad method " + m);
		}
	}

	/**
	 * Emit 3 new products for the pair {@code (y_{ij}, y_{ji})} when
	 * {@code y_{ii}} uses {@code mi} and {@code y_{jj}} uses {@code mj}
	 * (mi ≠ mj). Reuses the diagonal products via output-combination
	 * coefficients in W.
	 */
	private static void emitOffDiagonalPair(Builder b, int i, int j, int mi, int mj) {
		// Always emit such that mi < mj for symmetric handling.
		if (mi > mj) { int t=mi; mi=mj; mj=t; int ti=i; i=j; j=ti; }

		LinearA a_i_1 = LinearA.aEntry(i, 1, b.p, b.n);
		LinearA a_i_2 = LinearA.aEntry(i, 2, b.p, b.n);
		LinearA a_j_1 = LinearA.aEntry(j, 1, b.p, b.n);
		LinearA a_j_2 = LinearA.aEntry(j, 2, b.p, b.n);

		LinearB x_i_1 = LinearB.xEntry(1, i, b.p, b.n);
		LinearB x_i_2 = LinearB.xEntry(2, i, b.p, b.n);
		LinearB x_j_1 = LinearB.xEntry(1, j, b.p, b.n);
		LinearB x_j_2 = LinearB.xEntry(2, j, b.p, b.n);

		int outIJ = b.outIdx(i, j);
		int outJI = b.outIdx(j, i);

		// Note: diagonal product lookups happen lazily within each case branch,
		// because the diagonals only emit the products matching the method used
		// (method 1 → A, B; method 2 → C, D; method 3 → E, F).

		if (mi == 1 && mj == 2) {
			// y_{ij} = -B(a_i,x_i) - D(a_j,x_j) + F(a_i+a_j, x_i+x_j) - G(a_i, a_j, x_i, x_j)
			// y_{ji} = -A(a_i,x_i) + C(a_j,x_j) + E(a_i+a_j, x_i+x_j) + G(a_i, a_j, x_i, x_j)
			int kA_ii = b.findProduct(LinearA.of(a_i_2), x_i_1.add(x_i_2));
			int kB_ii = b.findProduct(a_i_1.sub(a_i_2), x_i_1);
			int kC_jj = b.findProduct(a_j_1.sub(a_j_2), x_j_2);
			int kD_jj = b.findProduct(LinearA.of(a_j_1), x_j_1.add(x_j_2));

			LinearA a_i_plus_j_1 = a_i_1.add(a_j_1);
			LinearA a_i_plus_j_2 = a_i_2.add(a_j_2);
			LinearB x_i_plus_j_1 = x_i_1.add(x_j_1);
			LinearB x_i_plus_j_2 = x_i_2.add(x_j_2);

			int kE = b.addProduct(LinearA.of(a_i_plus_j_2), x_i_plus_j_2);
			int kF = b.addProduct(LinearA.of(a_i_plus_j_1), x_i_plus_j_1);
			// G(a_i, a_j, x_i, x_j) = (a_{j,1} + a_{i,2}) * (x_{1,i} - x_{2,j})
			int kG = b.addProduct(a_j_1.add(a_i_2), x_i_1.sub(x_j_2));

			b.addToOutput(outIJ, kB_ii, -1.0);
			b.addToOutput(outIJ, kD_jj, -1.0);
			b.addToOutput(outIJ, kF,    +1.0);
			b.addToOutput(outIJ, kG,    -1.0);
			b.addToOutput(outJI, kA_ii, -1.0);
			b.addToOutput(outJI, kC_jj, +1.0);
			b.addToOutput(outJI, kE,    +1.0);
			b.addToOutput(outJI, kG,    +1.0);
		} else if (mi == 1 && mj == 3) {
			// y_{ij} = A(a_i,x_i) - D(-a_i+a_j, -x_i+x_j) + F(a_j, x_j) - G(a_i, -a_i+a_j, x_i, -x_i+x_j)
			// y_{ji} = B(a_i,x_i) + C(-a_i+a_j, -x_i+x_j) + E(a_j, x_j) + G(a_i, -a_i+a_j, x_i, -x_i+x_j)
			int kA_ii = b.findProduct(LinearA.of(a_i_2), x_i_1.add(x_i_2));
			int kB_ii = b.findProduct(a_i_1.sub(a_i_2), x_i_1);
			int kE_jj = b.findProduct(LinearA.of(a_j_2), x_j_2);
			int kF_jj = b.findProduct(LinearA.of(a_j_1), x_j_1);

			LinearA a_diff_1 = a_j_1.sub(a_i_1);
			LinearA a_diff_2 = a_j_2.sub(a_i_2);
			LinearB x_diff_1 = x_j_1.sub(x_i_1);
			LinearB x_diff_2 = x_j_2.sub(x_i_2);

			int kD_sub = b.addProduct(LinearA.of(a_diff_1), x_diff_1.add(x_diff_2));
			int kC_sub = b.addProduct(a_diff_1.sub(a_diff_2), x_diff_2);
			// G(a_i, -a_i+a_j, x_i, -x_i+x_j) = ((-a_{i,1}+a_{j,1}) + a_{i,2}) * (x_{1,i} - (-x_{2,i}+x_{2,j}))
			//                                  = (a_diff_1 + a_{i,2})            * (x_{1,i} - (x_diff_2))
			int kG = b.addProduct(a_diff_1.add(a_i_2), x_i_1.sub(x_diff_2));

			b.addToOutput(outIJ, kA_ii,  +1.0);
			b.addToOutput(outIJ, kD_sub, -1.0);
			b.addToOutput(outIJ, kF_jj,  +1.0);
			b.addToOutput(outIJ, kG,     -1.0);
			b.addToOutput(outJI, kB_ii,  +1.0);
			b.addToOutput(outJI, kC_sub, +1.0);
			b.addToOutput(outJI, kE_jj,  +1.0);
			b.addToOutput(outJI, kG,     +1.0);
		} else if (mi == 2 && mj == 3) {
			// y_{ij} = -A(-a_i+a_j, -x_i+x_j) + D(a_i, x_i) + E(a_j, x_j) + G(-a_i+a_j, a_i, -x_i+x_j, x_i)
			// y_{ji} = -B(-a_i+a_j, -x_i+x_j) - C(a_i, x_i) + F(a_j, x_j) - G(-a_i+a_j, a_i, -x_i+x_j, x_i)
			int kC_ii = b.findProduct(a_i_1.sub(a_i_2), x_i_2);
			int kD_ii = b.findProduct(LinearA.of(a_i_1), x_i_1.add(x_i_2));
			int kE_jj = b.findProduct(LinearA.of(a_j_2), x_j_2);
			int kF_jj = b.findProduct(LinearA.of(a_j_1), x_j_1);

			LinearA a_diff_1 = a_j_1.sub(a_i_1);
			LinearA a_diff_2 = a_j_2.sub(a_i_2);
			LinearB x_diff_1 = x_j_1.sub(x_i_1);
			LinearB x_diff_2 = x_j_2.sub(x_i_2);

			int kA_sub = b.addProduct(LinearA.of(a_diff_2), x_diff_1.add(x_diff_2));
			int kB_sub = b.addProduct(a_diff_1.sub(a_diff_2), x_diff_1);
			int kG = b.addProduct(a_i_1.add(a_diff_2), x_diff_1.sub(x_i_2));

			b.addToOutput(outIJ, kA_sub, -1.0);
			b.addToOutput(outIJ, kD_ii,  +1.0);
			b.addToOutput(outIJ, kE_jj,  +1.0);
			b.addToOutput(outIJ, kG,     +1.0);
			b.addToOutput(outJI, kB_sub, -1.0);
			b.addToOutput(outJI, kC_ii,  -1.0);
			b.addToOutput(outJI, kF_jj,  +1.0);
			b.addToOutput(outJI, kG,     -1.0);
		} else {
			throw new IllegalArgumentException("bad method pair: " + mi + ", " + mj);
		}
	}

	/**
	 * Page-10 "same-method" fallback. Both {@code y_{ii}} and
	 * {@code y_{jj}} use method {@code mij}, bridged by
	 * {@code y_{i+1,i+1}} using method {@code mBridge ≠ mij}.
	 *
	 * <p>Reuses the bridge-pair products E and F (emitted earlier when
	 * pair {@code (i, i+1)} was processed via Lemma 2) plus the
	 * diagonal A/B/C/D from {@code y_{jj}}. Emits 3 NEW products with
	 * substituted arguments. Total cost per same-method pair: 3 mults
	 * (vs 7 if computed independently).</p>
	 *
	 * <p>Currently implemented: case {@code (mij=1, mBridge=2)} — the
	 * explicit page-10 formula. The 5 other same-method/bridge
	 * combinations follow by symmetry (swapping the (A,B) ↔ (C,D) ↔
	 * (E,F) trios) and are tracked as TODO.</p>
	 */
	private static void emitSameMethodPair(Builder b, int i, int j, int mij, int mBridge) {
		emitSameMethodPair(b, i, j, mij, i + 1, mBridge);
	}

	/**
	 * Same-method pair with an EXPLICIT bridge position {@code bridgePos} (the
	 * linear callers pass {@code i+1}; the cyclic-banded caller passes a cyclic,
	 * method-selected position). Historically the position was recomputed here as
	 * {@code i+1} regardless of caller — in the cyclic builder this walked off the
	 * end of the index range (e.g. {@code a[n+1,·]} for the wrapped pair starting
	 * at {@code i=n}), the "expected product not found" failures of task #7.
	 */
	private static void emitSameMethodPair(Builder b, int i, int j, int mij, int bridgePos, int mBridge) {
		int ipp = bridgePos;
		if (mij == 1 && mBridge == 2) {
			emitSameMethodPair_11_bridge2(b, i, j, ipp);
			return;
		}
		if (mij == 2 && mBridge == 1) {
			emitSameMethodPair_22_bridge1(b, i, j, ipp);
			return;
		}
		if (mij == 1 && mBridge == 3) {
			emitSameMethodPair_11_bridge3(b, i, j, ipp);
			return;
		}
		if (mij == 2 && mBridge == 3) {
			emitSameMethodPair_22_bridge3(b, i, j, ipp);
			return;
		}
		// (3,3,bridge1) and (3,3,bridge2): can't arise in Case-1 odd-n (only one row gets
		// method 3 in the alternating coloring, so same-method-3 pairs are impossible).
		// Only relevant for Case-2 even-n via Lemma 3 sequences — sympy enumeration
		// returned 0 solutions with our standard atom catalog, suggesting these cases
		// need an extended set of candidate products. Tracked as future work.
		throw new UnsupportedOperationException(
				"same-method pair (" + i + "," + j + ") with mij=" + mij + ", bridge="
						+ mBridge + ": (3,3,*) cases not yet derivable from page-10 "
						+ "atoms (sympy enumeration returned 0 solutions; need extended candidates)");
	}

	/**
	 * Page-10 same-method case: both {@code y_{ii}}, {@code y_{jj}}
	 * use method 1 ({@code A+B}); bridge {@code y_{i+1,i+1}} uses
	 * method 2 ({@code -C+D}).
	 *
	 * <p>Formulas (after sympy verification — both pass):</p>
	 * <pre>
	 * y_{ij}+y_{i+1,j} = +B(a_j, x_j)
	 *                   +C(-a_j+a_i+a_{i+1}, -x_j+x_i+x_{i+1})
	 *                   +E( a_i+a_{i+1},      x_i+x_{i+1})    ← shared from pair (i, i+1)
	 *                   +G( a_j, -a_j+a_i+a_{i+1}, x_j, -x_j+x_i+x_{i+1})
	 *
	 * y_{ji}+y_{j,i+1} = +A(a_j, x_j)
	 *                   -D(-a_j+a_i+a_{i+1}, -x_j+x_i+x_{i+1})
	 *                   +F( a_i+a_{i+1},      x_i+x_{i+1})    ← shared from pair (i, i+1)
	 *                   -G( a_j, -a_j+a_i+a_{i+1}, x_j, -x_j+x_i+x_{i+1})
	 * </pre>
	 *
	 * <p>Then recover {@code y_{ij} = (y_{ij}+y_{i+1,j}) - y_{i+1,j}}
	 * (y_{i+1,j} already computed via Lemma 2 for pair (i+1, j)),
	 * similarly {@code y_{ji} = (y_{ji}+y_{j,i+1}) - y_{j,i+1}}.</p>
	 */
	private static void emitSameMethodPair_11_bridge2(Builder b, int i, int j, int ipp) {
		LinearA a_i_1 = LinearA.aEntry(i, 1, b.p, b.n);
		LinearA a_i_2 = LinearA.aEntry(i, 2, b.p, b.n);
		LinearA a_ipp_1 = LinearA.aEntry(ipp, 1, b.p, b.n);
		LinearA a_ipp_2 = LinearA.aEntry(ipp, 2, b.p, b.n);
		LinearA a_j_1 = LinearA.aEntry(j, 1, b.p, b.n);
		LinearA a_j_2 = LinearA.aEntry(j, 2, b.p, b.n);

		LinearB x_i_1 = LinearB.xEntry(1, i, b.p, b.n);
		LinearB x_i_2 = LinearB.xEntry(2, i, b.p, b.n);
		LinearB x_ipp_1 = LinearB.xEntry(1, ipp, b.p, b.n);
		LinearB x_ipp_2 = LinearB.xEntry(2, ipp, b.p, b.n);
		LinearB x_j_1 = LinearB.xEntry(1, j, b.p, b.n);
		LinearB x_j_2 = LinearB.xEntry(2, j, b.p, b.n);

		// Substituted args used in C(sub), D(sub), G(sub).
		LinearA sub_a_1 = a_i_1.add(a_ipp_1).sub(a_j_1);   // -a_j + a_i + a_{i+1}
		LinearA sub_a_2 = a_i_2.add(a_ipp_2).sub(a_j_2);
		LinearB sub_x_1 = x_i_1.add(x_ipp_1).sub(x_j_1);
		LinearB sub_x_2 = x_i_2.add(x_ipp_2).sub(x_j_2);

		// 3 NEW products.
		int kC_sub = b.addProduct(sub_a_1.sub(sub_a_2), sub_x_2);
		int kD_sub = b.addProduct(LinearA.of(sub_a_1), sub_x_1.add(sub_x_2));
		// G(a_j, sub, x_j, sub_x): arg1=a_j (sub2=a_{j,2}), arg2=sub_a (sub1=sub_a_1),
		// arg3=x_j (sub1=x_{1,j}), arg4=sub_x (sub2=sub_x_2)
		int kG = b.addProduct(sub_a_1.add(a_j_2), x_j_1.sub(sub_x_2));

		// Look up SHARED products. Already emitted by:
		// - y_{jj} method 1: A(a_j,x_j), B(a_j,x_j)
		// - pair (i, i+1) Lemma 2 case (1,2): E(a_i+a_{i+1}, x_i+x_{i+1}), F(...)
		int kA_jj = b.findProduct(LinearA.of(a_j_2), x_j_1.add(x_j_2));
		int kB_jj = b.findProduct(a_j_1.sub(a_j_2), x_j_1);
		int kE_adj = b.findProduct(LinearA.of(a_i_2.add(a_ipp_2)), x_i_2.add(x_ipp_2));
		int kF_adj = b.findProduct(LinearA.of(a_i_1.add(a_ipp_1)), x_i_1.add(x_ipp_1));

		// Build the SUM rows in W: outIdx(i, j) gets sum_{ij}, outIdx(j, i) gets sum_{ji}.
		// Then SUBTRACT y_{i+1, j} / y_{j, i+1} (computed by Lemma 2 for pair (i+1, j))
		// by NEGATING those W rows into ours.
		int outIJ = b.outIdx(i, j);
		int outJI = b.outIdx(j, i);
		int outIpp_J = b.outIdx(ipp, j);
		int outJ_Ipp = b.outIdx(j, ipp);

		// sum_{ij}: +B(a_j,x_j) + C(sub) + E(a_i+a_{i+1}, x_i+x_{i+1}) + G(sub)
		b.addToOutput(outIJ, kB_jj,  +1.0);
		b.addToOutput(outIJ, kC_sub, +1.0);
		b.addToOutput(outIJ, kE_adj, +1.0);
		b.addToOutput(outIJ, kG,     +1.0);

		// sum_{ji}: +A(a_j,x_j) - D(sub) + F(a_i+a_{i+1}, x_i+x_{i+1}) - G(sub)
		b.addToOutput(outJI, kA_jj,  +1.0);
		b.addToOutput(outJI, kD_sub, -1.0);
		b.addToOutput(outJI, kF_adj, +1.0);
		b.addToOutput(outJI, kG,     -1.0);

		// Now SUBTRACT y_{i+1, j} from outIJ — copy W[outIpp_J] with negative sign.
		Map<Integer, Double> ippJ = b.wAccum.get((long) outIpp_J);
		if (ippJ != null) {
			for (var e : ippJ.entrySet()) {
				b.addToOutput(outIJ, e.getKey(), -e.getValue());
			}
		}
		Map<Integer, Double> jIpp = b.wAccum.get((long) outJ_Ipp);
		if (jIpp != null) {
			for (var e : jIpp.entrySet()) {
				b.addToOutput(outJI, e.getKey(), -e.getValue());
			}
		}
	}

	/**
	 * Page-10 same-method case, derived by symbolic enumeration:
	 * both {@code y_{ii}}, {@code y_{jj}} use method 2 ({@code -C+D});
	 * bridge {@code y_{i+1,i+1}} uses method 1 ({@code A+B}).
	 *
	 * <p>Verified formulas (both PASS sympy):</p>
	 * <pre>
	 * y_{ij}+y_{i+1,j} = −C(a_j, x_j)
	 *                   +F(a_i+a_{i+1}, x_i+x_{i+1})         ← shared from pair (i, i+1)
	 *                   −B(sub_a_1, sub_a_2, sub_x_1)        ← NEW (substituted args)
	 *                   −G(sub_a_2, a_{j,1}, sub_x_1, x_{2,j}) ← NEW
	 *
	 * y_{ji}+y_{j,i+1} = +D(a_j, x_j)
	 *                   +E(a_i+a_{i+1}, x_i+x_{i+1})         ← shared from pair (i, i+1)
	 *                   −A(sub_a_2, sub_x_1, sub_x_2)        ← NEW
	 *                   +G(sub_a_2, a_{j,1}, sub_x_1, x_{2,j}) ← NEW (same as above)
	 * </pre>
	 *
	 * <p>3 NEW products (A_sub, B_sub, G), same count as (1,1,bridge2).
	 * Mirror-structure under method swap (A↔D, B↔C, E↔F) with G's
	 * argument order swapped.</p>
	 */
	private static void emitSameMethodPair_22_bridge1(Builder b, int i, int j, int ipp) {
		LinearA a_i_1 = LinearA.aEntry(i, 1, b.p, b.n);
		LinearA a_i_2 = LinearA.aEntry(i, 2, b.p, b.n);
		LinearA a_ipp_1 = LinearA.aEntry(ipp, 1, b.p, b.n);
		LinearA a_ipp_2 = LinearA.aEntry(ipp, 2, b.p, b.n);
		LinearA a_j_1 = LinearA.aEntry(j, 1, b.p, b.n);
		LinearA a_j_2 = LinearA.aEntry(j, 2, b.p, b.n);

		LinearB x_i_1 = LinearB.xEntry(1, i, b.p, b.n);
		LinearB x_i_2 = LinearB.xEntry(2, i, b.p, b.n);
		LinearB x_ipp_1 = LinearB.xEntry(1, ipp, b.p, b.n);
		LinearB x_ipp_2 = LinearB.xEntry(2, ipp, b.p, b.n);
		LinearB x_j_1 = LinearB.xEntry(1, j, b.p, b.n);
		LinearB x_j_2 = LinearB.xEntry(2, j, b.p, b.n);

		LinearA sub_a_1 = a_i_1.add(a_ipp_1).sub(a_j_1);
		LinearA sub_a_2 = a_i_2.add(a_ipp_2).sub(a_j_2);
		LinearB sub_x_1 = x_i_1.add(x_ipp_1).sub(x_j_1);
		LinearB sub_x_2 = x_i_2.add(x_ipp_2).sub(x_j_2);

		// 3 NEW products.
		int kA_sub = b.addProduct(LinearA.of(sub_a_2), sub_x_1.add(sub_x_2));
		int kB_sub = b.addProduct(sub_a_1.sub(sub_a_2), sub_x_1);
		// G_v2: arg1 = sub_a (use sub_a_2), arg2 = a_j (use a_{j,1}), arg3 = sub_x (use sub_x_1), arg4 = x_j (use x_{2,j})
		// G = (arg2.sub1 + arg1.sub2) * (arg3.sub1 - arg4.sub2) = (a_{j,1} + sub_a_2) * (sub_x_1 - x_{2,j})
		int kG = b.addProduct(a_j_1.add(sub_a_2), sub_x_1.sub(x_j_2));

		// Look up SHARED products from y_{jj} method 2 and pair (i, i+1) Lemma 2.
		int kC_jj = b.findProduct(a_j_1.sub(a_j_2), x_j_2);
		int kD_jj = b.findProduct(LinearA.of(a_j_1), x_j_1.add(x_j_2));
		int kE_adj = b.findProduct(LinearA.of(a_i_2.add(a_ipp_2)), x_i_2.add(x_ipp_2));
		int kF_adj = b.findProduct(LinearA.of(a_i_1.add(a_ipp_1)), x_i_1.add(x_ipp_1));

		int outIJ = b.outIdx(i, j);
		int outJI = b.outIdx(j, i);
		int outIpp_J = b.outIdx(ipp, j);
		int outJ_Ipp = b.outIdx(j, ipp);

		// sum_{ij}: -C_jj + F_iIpp - B_sub - G_v2
		b.addToOutput(outIJ, kC_jj,  -1.0);
		b.addToOutput(outIJ, kF_adj, +1.0);
		b.addToOutput(outIJ, kB_sub, -1.0);
		b.addToOutput(outIJ, kG,     -1.0);

		// sum_{ji}: +D_jj + E_iIpp - A_sub + G_v2
		b.addToOutput(outJI, kD_jj,  +1.0);
		b.addToOutput(outJI, kE_adj, +1.0);
		b.addToOutput(outJI, kA_sub, -1.0);
		b.addToOutput(outJI, kG,     +1.0);

		// Subtract y_{i+1, j} and y_{j, i+1} (already computed via Lemma 2).
		Map<Integer, Double> ippJ = b.wAccum.get((long) outIpp_J);
		if (ippJ != null) {
			for (var e : ippJ.entrySet()) {
				b.addToOutput(outIJ, e.getKey(), -e.getValue());
			}
		}
		Map<Integer, Double> jIpp = b.wAccum.get((long) outJ_Ipp);
		if (jIpp != null) {
			for (var e : jIpp.entrySet()) {
				b.addToOutput(outJI, e.getKey(), -e.getValue());
			}
		}
	}

	/**
	 * OPERATIONAL (1,1,bridge-3) same-method pair (task #7, 2026-06-11): both
	 * {@code y_ii}, {@code y_jj} method 1; bridge row {@code bp} method 3. The
	 * historical emitter relied on E/F sum-products that never exist when the
	 * bridge is genuinely method-3; this identity was derived over the TRUE
	 * 12-product reusable set (diagonals of i, j, bp + all products of the
	 * Lemma-2 pairs (i,bp) and (bp,j)) by
	 * {@code references/hopcroftkerr1971/sympy/derive_bridge_true_reusables.py}:
	 * <pre>
	 * y_ij = +A(i) +B(i) +F(bp) − D(d_ib) − G_ib + D(d_bj) − F(sub) + G_new
	 * y_ji = +A(i) +B(i) +E(bp) + C(d_ib) + G_ib − C(d_bj) − E(sub) − G_new
	 * </pre>
	 * with {@code d_ib = a_bp − a_i}, {@code d_bj = a_bp − a_j},
	 * {@code sub = a_i + a_bp − a_j},
	 * {@code G_new = (a_bp1 − a_j1 + a_i2)·(x_1i − x_2bp + x_2j)}.
	 */
	private static void emitSameMethodPair_11_bridge3(Builder b, int i, int j, int bp) {
		LinearA a_i_1 = LinearA.aEntry(i, 1, b.p, b.n);
		LinearA a_i_2 = LinearA.aEntry(i, 2, b.p, b.n);
		LinearA a_b_1 = LinearA.aEntry(bp, 1, b.p, b.n);
		LinearA a_b_2 = LinearA.aEntry(bp, 2, b.p, b.n);
		LinearA a_j_1 = LinearA.aEntry(j, 1, b.p, b.n);
		LinearA a_j_2 = LinearA.aEntry(j, 2, b.p, b.n);
		LinearB x_i_1 = LinearB.xEntry(1, i, b.p, b.n);
		LinearB x_i_2 = LinearB.xEntry(2, i, b.p, b.n);
		LinearB x_b_1 = LinearB.xEntry(1, bp, b.p, b.n);
		LinearB x_b_2 = LinearB.xEntry(2, bp, b.p, b.n);
		LinearB x_j_1 = LinearB.xEntry(1, j, b.p, b.n);
		LinearB x_j_2 = LinearB.xEntry(2, j, b.p, b.n);

		// Reused: diagonals of i (method 1) and bp (method 3).
		int kA_i = b.findProduct(LinearA.of(a_i_2), x_i_1.add(x_i_2));
		int kB_i = b.findProduct(a_i_1.sub(a_i_2), x_i_1);
		int kE_b = b.findProduct(LinearA.of(a_b_2), x_b_2);
		int kF_b = b.findProduct(LinearA.of(a_b_1), x_b_1);
		// Reused: pair (i,bp) methods (1,3) → D(d_ib), C(d_ib), G_ib.
		LinearA dib_1 = a_b_1.sub(a_i_1);
		LinearA dib_2 = a_b_2.sub(a_i_2);
		LinearB dxib_1 = x_b_1.sub(x_i_1);
		LinearB dxib_2 = x_b_2.sub(x_i_2);
		int kD_ib = b.findProduct(LinearA.of(dib_1), dxib_1.add(dxib_2));
		int kC_ib = b.findProduct(dib_1.sub(dib_2), dxib_2);
		int kG_ib = b.findProduct(dib_1.add(a_i_2), x_i_1.sub(dxib_2));
		// Reused: pair (bp,j) methods (3,1) → normalized (1,3) on (j,bp): D(d_bj), C(d_bj).
		LinearA dbj_1 = a_b_1.sub(a_j_1);
		LinearA dbj_2 = a_b_2.sub(a_j_2);
		LinearB dxbj_1 = x_b_1.sub(x_j_1);
		LinearB dxbj_2 = x_b_2.sub(x_j_2);
		int kD_bj = b.findProduct(LinearA.of(dbj_1), dxbj_1.add(dxbj_2));
		int kC_bj = b.findProduct(dbj_1.sub(dbj_2), dxbj_2);

		// 3 NEW products.
		LinearA sub_1 = a_i_1.add(a_b_1).sub(a_j_1);
		LinearA sub_2 = a_i_2.add(a_b_2).sub(a_j_2);
		LinearB subx_1 = x_i_1.add(x_b_1).sub(x_j_1);
		LinearB subx_2 = x_i_2.add(x_b_2).sub(x_j_2);
		int kE_sub = b.addProduct(LinearA.of(sub_2), subx_2);
		int kF_sub = b.addProduct(LinearA.of(sub_1), subx_1);
		int kG_new = b.addProduct(a_b_1.sub(a_j_1).add(a_i_2), x_i_1.sub(x_b_2).add(x_j_2));

		int outIJ = b.outIdx(i, j);
		int outJI = b.outIdx(j, i);
		b.addToOutput(outIJ, kA_i, +1.0);
		b.addToOutput(outIJ, kB_i, +1.0);
		b.addToOutput(outIJ, kF_b, +1.0);
		b.addToOutput(outIJ, kD_ib, -1.0);
		b.addToOutput(outIJ, kG_ib, -1.0);
		b.addToOutput(outIJ, kD_bj, +1.0);
		b.addToOutput(outIJ, kF_sub, -1.0);
		b.addToOutput(outIJ, kG_new, +1.0);
		b.addToOutput(outJI, kA_i, +1.0);
		b.addToOutput(outJI, kB_i, +1.0);
		b.addToOutput(outJI, kE_b, +1.0);
		b.addToOutput(outJI, kC_ib, +1.0);
		b.addToOutput(outJI, kG_ib, +1.0);
		b.addToOutput(outJI, kC_bj, -1.0);
		b.addToOutput(outJI, kE_sub, -1.0);
		b.addToOutput(outJI, kG_new, -1.0);
	}

	/**
	 * OPERATIONAL (2,2,bridge-3) same-method pair (task #7, 2026-06-11) — see
	 * {@link #emitSameMethodPair_11_bridge3}. Identity (true reusables):
	 * <pre>
	 * y_ij = −C(i) +D(i) +E(bp) − A(d_ib) + G_ib + A(d_bj) − E(sub) − G_new
	 * y_ji = −C(i) +D(i) +F(bp) − B(d_ib) − G_ib + B(d_bj) − F(sub) + G_new
	 * </pre>
	 * with {@code G_new = (a_i1 + a_bp2 − a_j2)·(x_1bp − x_1j − x_2i)}.
	 */
	private static void emitSameMethodPair_22_bridge3(Builder b, int i, int j, int bp) {
		LinearA a_i_1 = LinearA.aEntry(i, 1, b.p, b.n);
		LinearA a_i_2 = LinearA.aEntry(i, 2, b.p, b.n);
		LinearA a_b_1 = LinearA.aEntry(bp, 1, b.p, b.n);
		LinearA a_b_2 = LinearA.aEntry(bp, 2, b.p, b.n);
		LinearA a_j_1 = LinearA.aEntry(j, 1, b.p, b.n);
		LinearA a_j_2 = LinearA.aEntry(j, 2, b.p, b.n);
		LinearB x_i_1 = LinearB.xEntry(1, i, b.p, b.n);
		LinearB x_i_2 = LinearB.xEntry(2, i, b.p, b.n);
		LinearB x_b_1 = LinearB.xEntry(1, bp, b.p, b.n);
		LinearB x_b_2 = LinearB.xEntry(2, bp, b.p, b.n);
		LinearB x_j_1 = LinearB.xEntry(1, j, b.p, b.n);
		LinearB x_j_2 = LinearB.xEntry(2, j, b.p, b.n);

		// Reused: diagonals of i (method 2) and bp (method 3).
		int kC_i = b.findProduct(a_i_1.sub(a_i_2), x_i_2);
		int kD_i = b.findProduct(LinearA.of(a_i_1), x_i_1.add(x_i_2));
		int kE_b = b.findProduct(LinearA.of(a_b_2), x_b_2);
		int kF_b = b.findProduct(LinearA.of(a_b_1), x_b_1);
		// Reused: pair (i,bp) methods (2,3) → A(d_ib), B(d_ib), G_ib.
		LinearA dib_1 = a_b_1.sub(a_i_1);
		LinearA dib_2 = a_b_2.sub(a_i_2);
		LinearB dxib_1 = x_b_1.sub(x_i_1);
		LinearB dxib_2 = x_b_2.sub(x_i_2);
		int kA_ib = b.findProduct(LinearA.of(dib_2), dxib_1.add(dxib_2));
		int kB_ib = b.findProduct(dib_1.sub(dib_2), dxib_1);
		int kG_ib = b.findProduct(a_i_1.add(dib_2), dxib_1.sub(x_i_2));
		// Reused: pair (bp,j) methods (3,2) → normalized (2,3) on (j,bp): A(d_bj), B(d_bj).
		LinearA dbj_1 = a_b_1.sub(a_j_1);
		LinearA dbj_2 = a_b_2.sub(a_j_2);
		LinearB dxbj_1 = x_b_1.sub(x_j_1);
		LinearB dxbj_2 = x_b_2.sub(x_j_2);
		int kA_bj = b.findProduct(LinearA.of(dbj_2), dxbj_1.add(dxbj_2));
		int kB_bj = b.findProduct(dbj_1.sub(dbj_2), dxbj_1);

		// 3 NEW products.
		LinearA sub_1 = a_i_1.add(a_b_1).sub(a_j_1);
		LinearA sub_2 = a_i_2.add(a_b_2).sub(a_j_2);
		LinearB subx_1 = x_i_1.add(x_b_1).sub(x_j_1);
		LinearB subx_2 = x_i_2.add(x_b_2).sub(x_j_2);
		int kE_sub = b.addProduct(LinearA.of(sub_2), subx_2);
		int kF_sub = b.addProduct(LinearA.of(sub_1), subx_1);
		int kG_new = b.addProduct(a_i_1.add(a_b_2).sub(a_j_2), x_b_1.sub(x_j_1).sub(x_i_2));

		int outIJ = b.outIdx(i, j);
		int outJI = b.outIdx(j, i);
		b.addToOutput(outIJ, kC_i, -1.0);
		b.addToOutput(outIJ, kD_i, +1.0);
		b.addToOutput(outIJ, kE_b, +1.0);
		b.addToOutput(outIJ, kA_ib, -1.0);
		b.addToOutput(outIJ, kG_ib, +1.0);
		b.addToOutput(outIJ, kA_bj, +1.0);
		b.addToOutput(outIJ, kE_sub, -1.0);
		b.addToOutput(outIJ, kG_new, -1.0);
		b.addToOutput(outJI, kC_i, -1.0);
		b.addToOutput(outJI, kD_i, +1.0);
		b.addToOutput(outJI, kF_b, +1.0);
		b.addToOutput(outJI, kB_ib, -1.0);
		b.addToOutput(outJI, kG_ib, -1.0);
		b.addToOutput(outJI, kB_bj, +1.0);
		b.addToOutput(outJI, kF_sub, -1.0);
		b.addToOutput(outJI, kG_new, +1.0);
	}

	/** Subtract bridge product W rows: W[targetOut] -= W[bridgeOut]. */
	private static void subtractBridge(Builder b, int targetOut, int bridgeOut) {
		Map<Integer, Double> br = b.wAccum.get((long) bridgeOut);
		if (br == null) return;
		for (var e : br.entrySet()) {
			b.addToOutput(targetOut, e.getKey(), -e.getValue());
		}
	}

	// ─────────────────────────────────────────────────────────────────────
	// Linear forms over A entries and X entries.
	// ─────────────────────────────────────────────────────────────────────

	/**
	 * Linear form over the {@code p × 2} A-matrix: a map from
	 * {@code (i, j)} → coefficient. Two LinearA are equal iff they
	 * map to the same nonzero coefficient set, used as a key for
	 * sharing identical products.
	 */
	private static final class LinearA {
		private final Map<Long, Double> coeffs; // key = i*10 + j
		private LinearA(Map<Long, Double> c) { coeffs = c; }

		static LinearA aEntry(int i, int j, int p, int n) {
			Map<Long, Double> c = new HashMap<>();
			c.put(key(i, j), 1.0);
			return new LinearA(c);
		}
		static LinearA of(LinearA other) { return other; }
		static long key(int i, int j) { return ((long) i << 4) | j; }
		static int unkeyI(long k) { return (int) (k >> 4); }
		static int unkeyJ(long k) { return (int) (k & 0xF); }

		LinearA add(LinearA o) { return combine(o, +1.0); }
		LinearA sub(LinearA o) { return combine(o, -1.0); }
		private LinearA combine(LinearA o, double sign) {
			Map<Long, Double> out = new HashMap<>(coeffs);
			for (var e : o.coeffs.entrySet()) {
				out.merge(e.getKey(), sign * e.getValue(), Double::sum);
			}
			out.entrySet().removeIf(e -> Math.abs(e.getValue()) < 1e-12);
			return new LinearA(out);
		}

		@Override
		public int hashCode() { return coeffs.hashCode(); }
		@Override
		public boolean equals(Object o) {
			return o instanceof LinearA l && l.coeffs.equals(coeffs);
		}

		@Override
		public String toString() {
			StringBuilder sb = new StringBuilder();
			coeffs.entrySet().stream()
					.sorted(Map.Entry.comparingByKey())
					.forEach(e -> sb.append(e.getValue() > 0 ? "+" : "")
							.append(e.getValue() % 1 == 0 ? String.valueOf(e.getValue().longValue()) : e.getValue())
							.append("·a[").append(unkeyI(e.getKey())).append(',').append(unkeyJ(e.getKey())).append("] "));
			return sb.toString().trim();
		}
	}

	/**
	 * Linear form over the {@code 2 × n} X-matrix. Same shape as
	 * LinearA but indexed by {@code (row, col)} of X.
	 */
	private static final class LinearB {
		private final Map<Long, Double> coeffs;
		private LinearB(Map<Long, Double> c) { coeffs = c; }

		static LinearB xEntry(int row, int col, int p, int n) {
			Map<Long, Double> c = new HashMap<>();
			c.put(key(row, col), 1.0);
			return new LinearB(c);
		}
		static long key(int row, int col) { return ((long) row << 8) | col; }
		static int unkeyRow(long k) { return (int) (k >> 8); }
		static int unkeyCol(long k) { return (int) (k & 0xFF); }

		LinearB add(LinearB o) { return combine(o, +1.0); }
		LinearB sub(LinearB o) { return combine(o, -1.0); }
		private LinearB combine(LinearB o, double sign) {
			Map<Long, Double> out = new HashMap<>(coeffs);
			for (var e : o.coeffs.entrySet()) {
				out.merge(e.getKey(), sign * e.getValue(), Double::sum);
			}
			out.entrySet().removeIf(e -> Math.abs(e.getValue()) < 1e-12);
			return new LinearB(out);
		}

		@Override
		public int hashCode() { return coeffs.hashCode(); }
		@Override
		public boolean equals(Object o) {
			return o instanceof LinearB l && l.coeffs.equals(coeffs);
		}

		@Override
		public String toString() {
			StringBuilder sb = new StringBuilder();
			coeffs.entrySet().stream()
					.sorted(Map.Entry.comparingByKey())
					.forEach(e -> sb.append(e.getValue() > 0 ? "+" : "")
							.append(e.getValue() % 1 == 0 ? String.valueOf(e.getValue().longValue()) : e.getValue())
							.append("·x[").append(unkeyRow(e.getKey())).append(',').append(unkeyCol(e.getKey())).append("] "));
			return sb.toString().trim();
		}
	}

	// ─────────────────────────────────────────────────────────────────────
	// Accumulator.
	// ─────────────────────────────────────────────────────────────────────

	private static final class Builder {
		final int p;       // A is p × 2
		final int twoM;    // = 2 (middle dim)
		final int n;       // X is 2 × n
		final List<LinearA> uForms = new ArrayList<>();
		final List<LinearB> vForms = new ArrayList<>();
		final Map<Long, Map<Integer, Double>> wAccum = new HashMap<>(); // outRow → (productIdx → coef)
		// Sharing: index products by (u, v) pair so identical bilinear forms collapse.
		final Map<UVKey, Integer> productIndex = new HashMap<>();

		Builder(int p, int twoM, int n) {
			if (twoM != 2) throw new IllegalArgumentException("middle dim must be 2, got " + twoM);
			this.p = p; this.twoM = twoM; this.n = n;
		}

		int outIdx(int i, int l) {
			// Output Y is p×n, row-major: index = (i-1)*n + (l-1)
			return (i - 1) * n + (l - 1);
		}

		int addProduct(LinearA u, LinearB v) {
			UVKey key = new UVKey(u, v);
			Integer existing = productIndex.get(key);
			if (existing != null) return existing;
			int idx = uForms.size();
			uForms.add(u);
			vForms.add(v);
			productIndex.put(key, idx);
			return idx;
		}

		int findProduct(LinearA u, LinearB v) {
			Integer i = productIndex.get(new UVKey(u, v));
			if (i == null) throw new IllegalStateException("expected product not found: " + u + " · " + v);
			return i;
		}

		void addToOutput(int outRow, int productIdx, double coef) {
			wAccum.computeIfAbsent((long) outRow, k -> new HashMap<>())
					.merge(productIdx, coef, Double::sum);
		}

		NonCubicBilinearAlgorithm build() {
			int r = uForms.size();
			double[][] U = new double[p * twoM][r]; // A flat = i*2 + j (0-indexed)
			double[][] V = new double[twoM * n][r]; // X flat = (row-1)*n + (col-1)
			double[][] W = new double[p * n][r];

			for (int k = 0; k < r; k++) {
				for (var e : uForms.get(k).coeffs.entrySet()) {
					int i = LinearA.unkeyI(e.getKey()); // 1-based
					int j = LinearA.unkeyJ(e.getKey()); // 1-based
					U[(i - 1) * twoM + (j - 1)][k] = e.getValue();
				}
				for (var e : vForms.get(k).coeffs.entrySet()) {
					int row = LinearB.unkeyRow(e.getKey()); // 1-based
					int col = LinearB.unkeyCol(e.getKey()); // 1-based
					V[(row - 1) * n + (col - 1)][k] = e.getValue();
				}
			}
			for (var e : wAccum.entrySet()) {
				int outRow = e.getKey().intValue();
				for (var p : e.getValue().entrySet()) {
					W[outRow][p.getKey()] = p.getValue();
				}
			}
			return new NonCubicBilinearAlgorithm(p, twoM, n, U, V, W);
		}
	}

	private record UVKey(LinearA u, LinearB v) {}
}
