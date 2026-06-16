package eu.solven.matmul.search.flip;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import eu.solven.matmul.NonCubicBilinearAlgorithm;

/**
 * Mutable, exact-integer bilinear scheme for flip-graph walks (Kauers–Moosbauer
 * 2022 formalism; Perminov's ternary variant, arXiv:2511.20317). Each product
 * {@code l} holds three dense {@code int[]} factor vectors — {@code u} (n·m),
 * {@code v} (m·p), {@code w} (n·p) — i.e. the COLUMNS of the catalog's
 * {@code U/V/W}, stored product-major so flips mutate two small vectors.
 *
 * <p>Integer-only by design: flips, reductions and splits below all preserve
 * integrality, so a walk seeded from a {@code Z} scheme stays {@code Z}-exact at
 * every vertex — no field drift, no epsilon. Seeds with non-integer coefficients
 * are rejected at {@link #of(NonCubicBilinearAlgorithm)} (rational-coefficient
 * walks are a follow-up; would need exact rational vectors).</p>
 *
 * <p><strong>Sign-only sharing (load-bearing restriction):</strong> moves fire on
 * pairs whose shared vector is equal up to sign ({@code λ ∈ {+1,−1}}), not up to
 * an arbitrary scalar. For ternary ({@code ZT}) schemes this is complete —
 * proportional ternary vectors are always ±equal — and it keeps every move
 * integral. General integer-λ flips are only integral in one variant direction,
 * so they are deferred rather than half-supported.</p>
 */
public final class FlipScheme {

	/** Tensor slot whose vectors two products share (the flip pivot). */
	public enum Slot { U, V, W }

	/** Format axis for the meta-moves ({@link #extendAxis}, {@link #projectAxis}). */
	public enum Axis { N, M, P }

	public final int n;
	public final int m;
	public final int p;

	// Parallel per-product lists; removal is swap-with-last (order is not
	// semantically meaningful — a scheme is a SET of rank-one terms).
	private final List<int[]> us = new ArrayList<>();
	private final List<int[]> vs = new ArrayList<>();
	private final List<int[]> ws = new ArrayList<>();

	private FlipScheme(int n, int m, int p) {
		this.n = n;
		this.m = m;
		this.p = p;
	}

	/** @throws IllegalArgumentException on any non-integer coefficient. */
	public static FlipScheme of(NonCubicBilinearAlgorithm a) {
		FlipScheme s = new FlipScheme(a.n, a.m, a.p);
		double[][] du = a.denseU();
		double[][] dv = a.denseV();
		double[][] dw = a.denseW();
		for (int l = 0; l < a.r; l++) {
			s.us.add(column(du, l));
			s.vs.add(column(dv, l));
			s.ws.add(column(dw, l));
		}
		return s;
	}

	private static int[] column(double[][] dense, int col) {
		int[] out = new int[dense.length];
		for (int i = 0; i < dense.length; i++) {
			double v = dense[i][col];
			long iv = Math.round(v);
			if (Math.abs(v - iv) > 1e-9 || Math.abs(iv) > Integer.MAX_VALUE) {
				throw new IllegalArgumentException(
						"non-integer coefficient " + v + " — integer/ZT seeds only");
			}
			out[i] = (int) iv;
		}
		return out;
	}

	public int rank() {
		return us.size();
	}

	public int[] u(int l) { return us.get(l); }
	public int[] v(int l) { return vs.get(l); }
	public int[] w(int l) { return ws.get(l); }

	public int[] vec(Slot s, int l) {
		return switch (s) { case U -> us.get(l); case V -> vs.get(l); case W -> ws.get(l); };
	}

	public FlipScheme copy() {
		FlipScheme c = new FlipScheme(n, m, p);
		for (int l = 0; l < rank(); l++) {
			c.us.add(us.get(l).clone());
			c.vs.add(vs.get(l).clone());
			c.ws.add(ws.get(l).clone());
		}
		return c;
	}

	public NonCubicBilinearAlgorithm toAlgorithm() {
		int r = rank();
		double[][] du = new double[n * m][r];
		double[][] dv = new double[m * p][r];
		double[][] dw = new double[n * p][r];
		for (int l = 0; l < r; l++) {
			fill(du, us.get(l), l);
			fill(dv, vs.get(l), l);
			fill(dw, ws.get(l), l);
		}
		return new NonCubicBilinearAlgorithm(n, m, p, du, dv, dw);
	}

	private static void fill(double[][] dense, int[] vec, int col) {
		for (int i = 0; i < vec.length; i++) {
			dense[i][col] = vec[i];
		}
	}

	// ── sign-proportionality ────────────────────────────────────────────────

	/** {@code +1} if {@code b == a}, {@code −1} if {@code b == −a}, else {@code 0}
	 *  (zero vectors return 0 — they are reduction fodder, not flip pivots). */
	public static int signRatio(int[] a, int[] b) {
		int sign = 0;
		for (int i = 0; i < a.length; i++) {
			if (a[i] == 0 && b[i] == 0) {
				continue;
			}
			if (a[i] == b[i]) {
				if (sign < 0) return 0;
				sign = 1;
			} else if (a[i] == -b[i]) {
				if (sign > 0) return 0;
				sign = -1;
			} else {
				return 0;
			}
		}
		return sign;
	}

	public static boolean isZero(int[] a) {
		for (int x : a) {
			if (x != 0) return false;
		}
		return true;
	}

	/**
	 * Groups of ≥2 products whose {@code slot} vectors are equal up to sign
	 * (zero vectors excluded). These are exactly the flip pivots — and, for
	 * ternary schemes, exactly the slot's buds.
	 */
	public List<int[]> signClasses(Slot slot) {
		Map<VecKey, List<Integer>> byKey = new HashMap<>();
		for (int l = 0; l < rank(); l++) {
			int[] vec = vec(slot, l);
			if (isZero(vec)) {
				continue;
			}
			byKey.computeIfAbsent(new VecKey(signNormalize(vec)), k -> new ArrayList<>()).add(l);
		}
		List<int[]> out = new ArrayList<>();
		for (List<Integer> g : byKey.values()) {
			if (g.size() >= 2) {
				out.add(g.stream().mapToInt(Integer::intValue).toArray());
			}
		}
		return out;
	}

	/** Copy of {@code vec} negated if its first nonzero entry is negative. */
	static int[] signNormalize(int[] vec) {
		for (int x : vec) {
			if (x != 0) {
				if (x > 0) return vec;
				int[] neg = new int[vec.length];
				for (int i = 0; i < vec.length; i++) {
					neg[i] = -vec[i];
				}
				return neg;
			}
		}
		return vec;
	}

	/** Value-semantics wrapper so int[] vectors can key a HashMap. */
	record VecKey(int[] vec) {
		@Override
		public boolean equals(Object o) {
			return o instanceof VecKey k && java.util.Arrays.equals(vec, k.vec);
		}

		@Override
		public int hashCode() {
			return java.util.Arrays.hashCode(vec);
		}
	}

	// ── moves ───────────────────────────────────────────────────────────────

	/**
	 * Flip on products {@code (i,j)} sharing {@code slot} up to sign
	 * ({@code vec(slot,j) == sign · vec(slot,i)}, caller-supplied). With the two
	 * non-shared slots in canonical order {@code (X,Y)} (U→(V,W), V→(U,W),
	 * W→(U,V)):
	 * <ul>
	 *   <li>variant A: {@code X_i += sign·X_j;  Y_j −= Y_i}</li>
	 *   <li>variant B: {@code Y_i += sign·Y_j;  X_j −= X_i}</li>
	 * </ul>
	 * Both rewrite the two rank-one terms without changing their sum, so the
	 * represented tensor is invariant by construction. Throws if the pair does
	 * not actually share the slot at {@code sign} (a silent no-op here would be
	 * a silent search corruption).
	 */
	public void flip(Slot slot, int i, int j, int sign, boolean variantB) {
		if (i == j) {
			throw new IllegalArgumentException("flip needs two distinct products, got " + i);
		}
		if (sign == 0 || signRatio(vec(slot, i), vec(slot, j)) != sign) {
			throw new IllegalArgumentException(
					"products " + i + "," + j + " do not share slot " + slot + " at sign " + sign);
		}
		Slot recv = recvSlot(slot, variantB);
		Slot give = giveSlot(slot, variantB);
		addInto(vec(recv, i), vec(recv, j), sign);
		addInto(vec(give, j), vec(give, i), -1);
	}

	/**
	 * {@link #flip} guarded by the alphabet cap: applies the flip, and if any
	 * coefficient of the two modified vectors exceeds {@code cap} (in absolute
	 * value), reverts it exactly — flips are invertible because the shared-slot
	 * vectors and the donor vectors are untouched. Returns whether the flip
	 * stuck. {@code cap ≤ 0} = unbounded (never reverts).
	 */
	public boolean flipWithinCap(Slot slot, int i, int j, int sign, boolean variantB, int cap) {
		flip(slot, i, j, sign, variantB);
		if (cap <= 0) {
			return true;
		}
		Slot recv = recvSlot(slot, variantB);
		Slot give = giveSlot(slot, variantB);
		if (maxAbs(vec(recv, i)) <= cap && maxAbs(vec(give, j)) <= cap) {
			return true;
		}
		addInto(vec(recv, i), vec(recv, j), -sign);
		addInto(vec(give, j), vec(give, i), 1);
		return false;
	}

	private static Slot recvSlot(Slot shared, boolean variantB) {
		Slot x = switch (shared) { case U -> Slot.V; case V, W -> Slot.U; };
		Slot y = (shared == Slot.W) ? Slot.V : Slot.W;
		return variantB ? y : x;
	}

	private static Slot giveSlot(Slot shared, boolean variantB) {
		return recvSlot(shared, !variantB);
	}

	private static void addInto(int[] target, int[] src, int factor) {
		for (int k = 0; k < target.length; k++) {
			target[k] += factor * src[k];
		}
	}

	/**
	 * One reduction pass: drop products with a zero vector, and merge any pair
	 * sharing TWO slots up to sign ({@code u_j = s₁·u_i, v_j = s₂·v_i ⇒
	 * w_i += s₁s₂·w_j}, drop {@code j} — same for the other slot pairs).
	 * Returns whether anything changed; {@link #reduce()} loops to fixpoint.
	 */
	public boolean reduceOnce() {
		boolean changed = dropZeroProducts();
		return mergeOnce() || changed;
	}

	/** Drop every product with a zero vector (they contribute nothing). Always
	 *  safe for any objective: rank drops, no structure score can decrease. */
	public boolean dropZeroProducts() {
		boolean changed = false;
		for (int l = rank() - 1; l >= 0; l--) {
			if (isZero(us.get(l)) || isZero(vs.get(l)) || isZero(ws.get(l))) {
				removeProduct(l);
				changed = true;
			}
		}
		return changed;
	}

	/**
	 * Merge ONE pair sharing two slots up to sign, if any. Exposed separately
	 * from {@link #dropZeroProducts} because a merge is NOT neutral for
	 * structure objectives — it consumes exactly the doubly-proportional pairs
	 * a bud walk tries to build — so walks may want it cost-gated rather than
	 * automatic.
	 */
	public boolean mergeOnce() {
		for (int i = 0; i < rank(); i++) {
			for (int j = i + 1; j < rank(); j++) {
				int su = signRatio(us.get(i), us.get(j));
				int sv = signRatio(vs.get(i), vs.get(j));
				int sw = signRatio(ws.get(i), ws.get(j));
				if (su != 0 && sv != 0) {
					addInto(ws.get(i), ws.get(j), su * sv);
				} else if (su != 0 && sw != 0) {
					addInto(vs.get(i), vs.get(j), su * sw);
				} else if (sv != 0 && sw != 0) {
					addInto(us.get(i), us.get(j), sv * sw);
				} else {
					continue;
				}
				removeProduct(j);
				return true;
			}
		}
		return false;
	}

	/** Reduce to fixpoint (zero-drops + two-slot merges). Rank only decreases. */
	public void reduce() {
		while (reduceOnce()) {
			// loop
		}
	}

	// ── meta-moves (Kauers–Wood meta flip graph, arXiv:2510.19787) ──────────

	/**
	 * Extend the format by one along {@code axis} — the meta-move that grows
	 * ⟨n,m,p⟩ to a neighbouring format. Existing vectors are re-indexed into the
	 * larger flatten; the new slice is covered by NAIVE products (one elementary
	 * product per entry it interacts with), so the result is exact by
	 * construction. Rank grows by {@code m·p} (N), {@code n·p} (M) or
	 * {@code n·m} (P). Returns a NEW scheme; {@code this} is untouched.
	 */
	public FlipScheme extendAxis(Axis axis) {
		int n2 = n + (axis == Axis.N ? 1 : 0);
		int m2 = m + (axis == Axis.M ? 1 : 0);
		int p2 = p + (axis == Axis.P ? 1 : 0);
		FlipScheme out = new FlipScheme(n2, m2, p2);
		for (int l = 0; l < rank(); l++) {
			out.us.add(grow(us.get(l), n, m, n2, m2));
			out.vs.add(grow(vs.get(l), m, p, m2, p2));
			out.ws.add(grow(ws.get(l), n, p, n2, p2));
		}
		switch (axis) {
			case N -> {
				// New output row: C[n,k] = Σ_j A[n,j]·B[j,k]
				for (int j = 0; j < m; j++) {
					for (int k = 0; k < p; k++) {
						out.addProduct(unit(n2 * m2, n * m2 + j),
								unit(m2 * p2, j * p2 + k),
								unit(n2 * p2, n * p2 + k));
					}
				}
			}
			case M -> {
				// New inner index: C[i,k] += A[i,m]·B[m,k]
				for (int i = 0; i < n; i++) {
					for (int k = 0; k < p; k++) {
						out.addProduct(unit(n2 * m2, i * m2 + m),
								unit(m2 * p2, m * p2 + k),
								unit(n2 * p2, i * p2 + k));
					}
				}
			}
			case P -> {
				// New output column: C[i,p] = Σ_j A[i,j]·B[j,p]
				for (int i = 0; i < n; i++) {
					for (int j = 0; j < m; j++) {
						out.addProduct(unit(n2 * m2, i * m2 + j),
								unit(m2 * p2, j * p2 + p),
								unit(n2 * p2, i * p2 + p));
					}
				}
			}
		}
		return out;
	}

	/**
	 * Project the format down by dropping index {@code drop} on {@code axis} —
	 * the meta-move that shrinks ⟨n,m,p⟩. Vectors are re-indexed into the
	 * smaller flatten and products whose support dies are dead-code-eliminated
	 * (the {@code R → R − μ} operator the projection-margin objective feeds).
	 * Exact by construction. Returns a NEW scheme; {@code this} is untouched.
	 */
	public FlipScheme projectAxis(Axis axis, int drop) {
		int dim = switch (axis) { case N -> n; case M -> m; case P -> p; };
		if (dim <= 1) {
			throw new IllegalArgumentException("cannot project axis " + axis + " below 1");
		}
		if (drop < 0 || drop >= dim) {
			throw new IllegalArgumentException("drop index " + drop + " out of range for " + axis);
		}
		int n2 = n - (axis == Axis.N ? 1 : 0);
		int m2 = m - (axis == Axis.M ? 1 : 0);
		int p2 = p - (axis == Axis.P ? 1 : 0);
		FlipScheme out = new FlipScheme(n2, m2, p2);
		for (int l = 0; l < rank(); l++) {
			int[] u2 = switch (axis) {
				case N -> dropRow(us.get(l), n, m, drop);
				case M -> dropCol(us.get(l), n, m, drop);
				case P -> us.get(l).clone();
			};
			int[] v2 = switch (axis) {
				case N -> vs.get(l).clone();
				case M -> dropRow(vs.get(l), m, p, drop);
				case P -> dropCol(vs.get(l), m, p, drop);
			};
			int[] w2 = switch (axis) {
				case N -> dropRow(ws.get(l), n, p, drop);
				case M -> ws.get(l).clone();
				case P -> dropCol(ws.get(l), n, p, drop);
			};
			out.addProduct(u2, v2, w2);
		}
		out.dropZeroProducts();
		return out;
	}

	private void addProduct(int[] u, int[] v, int[] w) {
		us.add(u);
		vs.add(v);
		ws.add(w);
	}

	/** Re-flatten a (rows×cols) vector into a (rows2×cols2) flatten, rows2≥rows, cols2≥cols. */
	private static int[] grow(int[] vec, int rows, int cols, int rows2, int cols2) {
		int[] out = new int[rows2 * cols2];
		for (int r = 0; r < rows; r++) {
			for (int c = 0; c < cols; c++) {
				out[r * cols2 + c] = vec[r * cols + c];
			}
		}
		return out;
	}

	private static int[] dropRow(int[] vec, int rows, int cols, int d) {
		int[] out = new int[(rows - 1) * cols];
		for (int r = 0; r < rows; r++) {
			if (r == d) {
				continue;
			}
			int r2 = r < d ? r : r - 1;
			System.arraycopy(vec, r * cols, out, r2 * cols, cols);
		}
		return out;
	}

	private static int[] dropCol(int[] vec, int rows, int cols, int d) {
		int[] out = new int[rows * (cols - 1)];
		for (int r = 0; r < rows; r++) {
			for (int c = 0; c < cols; c++) {
				if (c == d) {
					continue;
				}
				out[r * (cols - 1) + (c < d ? c : c - 1)] = vec[r * cols + c];
			}
		}
		return out;
	}

	private static int[] unit(int len, int idx) {
		int[] out = new int[len];
		out[idx] = 1;
		return out;
	}

	private void removeProduct(int l) {
		int last = rank() - 1;
		us.set(l, us.get(last));
		vs.set(l, vs.get(last));
		ws.set(l, ws.get(last));
		us.remove(last);
		vs.remove(last);
		ws.remove(last);
	}

	/**
	 * Plus transition (Moosbauer–Poole): split one product into two by
	 * partitioning the nonzero support of one of its vectors — rank+1, tensor
	 * unchanged, and alphabet-safe (each part is a mask of the original, so no
	 * coefficient grows). Returns {@code false} if no product has a vector with
	 * ≥2 nonzeros (nothing splittable).
	 */
	public boolean split(Random rng) {
		int r = rank();
		int start = rng.nextInt(r);
		Slot[] slots = Slot.values();
		int slotStart = rng.nextInt(3);
		for (int dl = 0; dl < r; dl++) {
			int l = (start + dl) % r;
			for (int ds = 0; ds < 3; ds++) {
				Slot slot = slots[(slotStart + ds) % 3];
				int[] vec = vec(slot, l);
				List<Integer> nz = new ArrayList<>();
				for (int k = 0; k < vec.length; k++) {
					if (vec[k] != 0) nz.add(k);
				}
				if (nz.size() < 2) {
					continue;
				}
				// Random non-trivial subset of the support moves to the new product.
				int[] part = new int[vec.length];
				int moved = 0;
				for (int k : nz) {
					if (rng.nextBoolean()) {
						part[k] = vec[k];
						moved++;
					}
				}
				if (moved == 0) {
					int k = nz.get(rng.nextInt(nz.size()));
					part[k] = vec[k];
					moved = 1;
				}
				if (moved == nz.size()) {
					int k = nz.get(rng.nextInt(nz.size()));
					part[k] = 0;
				}
				for (int k = 0; k < vec.length; k++) {
					if (part[k] != 0) vec[k] = 0;
				}
				int[] nu = (slot == Slot.U) ? part : us.get(l).clone();
				int[] nv = (slot == Slot.V) ? part : vs.get(l).clone();
				int[] nw = (slot == Slot.W) ? part : ws.get(l).clone();
				us.add(nu);
				vs.add(nv);
				ws.add(nw);
				return true;
			}
		}
		return false;
	}

	/**
	 * EXACT integer Brent-equation check: the represented tensor equals the
	 * ⟨n,m,p⟩ matmul tensor, computed in {@code long} arithmetic — immune to the
	 * floating-point cancellation that makes {@code Verifier} unreliable once an
	 * unbounded walk grows large coefficients. O(r·(nm)(mp)(np)) — test/boundary
	 * use, not per-step. Overflows (silently) only past |coef|³·r ≈ 2⁶³, far
	 * beyond any sane walk alphabet.
	 */
	public boolean isExactIntTensor() {
		int r = rank();
		for (int a = 0; a < n * m; a++) {
			int i = a / m;
			int j = a % m;
			for (int b = 0; b < m * p; b++) {
				int j2 = b / p;
				int k = b % p;
				for (int c = 0; c < n * p; c++) {
					int i2 = c / p;
					int k2 = c % p;
					long sum = 0;
					for (int l = 0; l < r; l++) {
						sum += (long) us.get(l)[a] * vs.get(l)[b] * ws.get(l)[c];
					}
					long expected = (i == i2 && j == j2 && k == k2) ? 1 : 0;
					if (sum != expected) {
						return false;
					}
				}
			}
		}
		return true;
	}

	/** Max |coefficient| over all vectors — the alphabet check ({@code ≤1} = ZT). */
	public int maxAbsCoefficient() {
		int max = 0;
		for (int l = 0; l < rank(); l++) {
			max = Math.max(max, maxAbs(us.get(l)));
			max = Math.max(max, maxAbs(vs.get(l)));
			max = Math.max(max, maxAbs(ws.get(l)));
		}
		return max;
	}

	private static int maxAbs(int[] vec) {
		int max = 0;
		for (int x : vec) {
			max = Math.max(max, Math.abs(x));
		}
		return max;
	}
}
