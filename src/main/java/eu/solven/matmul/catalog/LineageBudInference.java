package eu.solven.matmul.catalog;

import eu.solven.matmul.recombination.Recombination;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Function;

import eu.solven.matmul.NonCubicBilinearAlgorithm;

/**
 * Infer a scheme's bud structure from its construction {@link Lineage} —
 * <em>without expanding the (possibly huge) composite scheme</em>. This
 * promotes buds from a validate-tier quantity (needs explicit U/V/W) to a
 * detect-tier one (ref-level), just like rank: the only thing that ever needs
 * expansion is the small {@link Lineage.Atom} at the leaves; everything
 * composed on top is propagated by the per-operator rules below.
 *
 * <p>The bud structure is carried as three independent class-size multisets
 * (U, V, W) — see {@link SerendipitousBudProduct#independentClassSizes} — which
 * is the composition-stable notion (it propagates exactly through Kronecker and
 * tensor-symmetry relabelling).</p>
 *
 * <h3>Per-operator rules &amp; certainty</h3>
 * <table>
 *   <caption>How each lineage op acts on buds</caption>
 *   <tr><th>op</th><th>rule</th><th>certainty</th></tr>
 *   <tr><td>{@code Atom}</td><td>expand (small) → independent class sizes</td><td>EXACT</td></tr>
 *   <tr><td>{@code KronProduct/KronChain}</td><td>per type, Cartesian product of
 *       class sizes ({@code a·b}); cancellation-free</td><td>EXACT</td></tr>
 *   <tr><td>{@code Transpose}</td><td>relabel U/V/W partitions per the perm</td><td>EXACT</td></tr>
 *   <tr><td>{@code AxisFlip / AxisPermute / Dce}</td><td>row reorder / zero-term
 *       removal preserves proportionality classes → passthrough</td><td>EXACT</td></tr>
 *   <tr><td>{@code ConcatCols / ConcatRows / SumInner}</td><td>per type, multiset
 *       <em>union</em> of operand classes — exact on the tiled axis, a lower bound
 *       on the other two (cross-operand merges not predictable from sizes)</td>
 *       <td>STRUCTURAL_ESTIMATE</td></tr>
 *   <tr><td>{@code Recombination*, AugmentSquareDiscard, SerendipitousProduct,
 *       Project, DisjointSum}</td><td>rank-reducing / zero-sparing /
 *       allocation-dependent — not reliably inferable from sizes alone</td>
 *       <td>UNKNOWN (expansion required)</td></tr>
 * </table>
 *
 * <p>Certainty combines as the lattice minimum over the subtree
 * (UNKNOWN &lt; STRUCTURAL_ESTIMATE &lt; EXACT): one estimate or unknown leaf
 * downgrades the whole. This is the optimality-discipline label — a caller must
 * not present a {@code STRUCTURAL_ESTIMATE}/{@code UNKNOWN} profile as the
 * certified bud structure.</p>
 */
public final class LineageBudInference {

	private LineageBudInference() {}

	/** Honesty tier of an inferred profile (optimality discipline). */
	public enum Certainty {
		/** Proven from the lineage by cancellation-free / relabelling rules. */
		EXACT,
		/** A structural bound (e.g. concat union ignores cross-operand merges). */
		STRUCTURAL_ESTIMATE,
		/** Not inferable from the lineage; expanding the scheme is required. */
		UNKNOWN;

		/** Lattice minimum: the weakest of the two. */
		static Certainty min(Certainty a, Certainty b) {
			if (a == UNKNOWN || b == UNKNOWN) return UNKNOWN;
			if (a == STRUCTURAL_ESTIMATE || b == STRUCTURAL_ESTIMATE) return STRUCTURAL_ESTIMATE;
			return EXACT;
		}
	}

	/**
	 * Inferred bud structure. {@code uClasses/vClasses/wClasses} are the
	 * (descending) class-size multisets per factor; each sums to {@code rank}
	 * (when {@code certainty != UNKNOWN}). A class of size ≥ 2 is a bud.
	 */
	public record Profile(int rank, int[] uClasses, int[] vClasses, int[] wClasses,
			Certainty certainty) {

		static final Profile UNKNOWN = new Profile(-1, new int[0], new int[0], new int[0],
				Certainty.UNKNOWN);

		public boolean known() {
			return certainty != Certainty.UNKNOWN;
		}

		public int uBuds() { return countBuds(uClasses); }
		public int vBuds() { return countBuds(vClasses); }
		public int wBuds() { return countBuds(wClasses); }

		public boolean hasBuds() {
			return uBuds() + vBuds() + wBuds() > 0;
		}

		private static int countBuds(int[] classes) {
			int n = 0;
			for (int s : classes) if (s >= 2) n++;
			return n;
		}

		/** Render as {@code "k×U⟨1,1,s⟩ + …"}, matching {@link SerendipitousBudProduct} tags. */
		public String summary() {
			java.util.TreeMap<String, Integer> groups = new java.util.TreeMap<>();
			for (int s : uClasses) if (s >= 2) groups.merge("U⟨1,1," + s + "⟩", 1, Integer::sum);
			for (int s : vClasses) if (s >= 2) groups.merge("V⟨" + s + ",1,1⟩", 1, Integer::sum);
			for (int s : wClasses) if (s >= 2) groups.merge("W⟨1," + s + ",1⟩", 1, Integer::sum);
			StringBuilder sb = new StringBuilder();
			groups.forEach((tag, c) -> {
				if (sb.length() > 0) sb.append(" + ");
				sb.append(c).append("×").append(tag);
			});
			return sb.toString();
		}
	}

	/** Reference profile from an already-expanded scheme (the certified ground truth). */
	public static Profile fromExpanded(NonCubicBilinearAlgorithm a) {
		int[][] cs = SerendipitousBudProduct.independentClassSizes(a);
		return new Profile(a.r, cs[0], cs[1], cs[2], Certainty.EXACT);
	}

	/**
	 * Infer the bud profile of {@code node} by propagating up the lineage.
	 * {@code atomResolver} maps an {@link Lineage.Atom} ref to its expanded
	 * (small) scheme; returning {@code null} for an unresolvable atom yields an
	 * {@code UNKNOWN} profile. The resolver is the <em>only</em> place expansion
	 * happens.
	 */
	public static Profile infer(Lineage.Node node,
			Function<String, NonCubicBilinearAlgorithm> atomResolver) {
		return infer(node, atomResolver, sz -> Profile.UNKNOWN);
	}

	/**
	 * As {@link #infer(Lineage.Node, Function)} but with a {@code leafShape}
	 * resolver — {@code (n,m,p) → Profile} for the catalog-best scheme at a
	 * sub-shape — which enables the {@code RecombinationN} estimate (the only op
	 * that needs leaf buds <em>by shape</em>, since recombination leaves are
	 * looked up by sub-shape, not by atom ref). Returning a {@code null}/UNKNOWN
	 * leaf profile leaves the recombination UNKNOWN.
	 */
	public static Profile infer(Lineage.Node node,
			Function<String, NonCubicBilinearAlgorithm> atomResolver,
			Function<int[], Profile> leafShape) {
		switch (node) {
			case Lineage.Atom a -> {
				NonCubicBilinearAlgorithm alg = atomResolver.apply(a.ref());
				return alg == null ? Profile.UNKNOWN : fromExpanded(alg);
			}
			case Lineage.KronProduct kp -> {
				return kron(infer(kp.outer(), atomResolver, leafShape),
						infer(kp.inner(), atomResolver, leafShape));
			}
			case Lineage.KronChain kc -> {
				Profile acc = null;
				for (Lineage.Node f : kc.factors()) {
					Profile p = infer(f, atomResolver, leafShape);
					acc = (acc == null) ? p : kron(acc, p);
				}
				return acc == null ? Profile.UNKNOWN : acc;
			}
			case Lineage.Transpose t -> {
				return transpose(infer(t.child(), atomResolver, leafShape), t.perm());
			}
			// OrientAs is an S₃ axis relabel of the child — passthrough (it preserves
			// the proportionality classes of nonzero terms, like AxisPermute below).
			case Lineage.OrientAs o -> {
				return infer(o.child(), atomResolver, leafShape);
			}
			// PeeledViaTa: a TA cross-fusion assembly — bud structure not inferred
			// (conservative; UNKNOWN leaves the recombination bud profile UNKNOWN).
			case Lineage.PeeledViaTa t -> {
				return Profile.UNKNOWN;
			}
			// Generic TA-fused recombination — bud structure not inferred (conservative).
			case Lineage.RecombinationTaN r -> {
				return Profile.UNKNOWN;
			}
			// Row reorder / flip / dead-zero-term removal all preserve the
			// proportionality classes of the surviving nonzero terms → passthrough.
			case Lineage.AxisFlip af -> {
				return infer(af.child(), atomResolver, leafShape);
			}
			case Lineage.AxisPermute ap -> {
				return infer(ap.child(), atomResolver, leafShape);
			}
			case Lineage.Dce d -> {
				return infer(d.child(), atomResolver, leafShape);
			}
			// Direct-sum-like: union of class multisets — exact on the tiled axis,
			// a lower bound elsewhere (cross-operand merges unpredictable) → estimate.
			case Lineage.ConcatCols c -> {
				return concat(infer(c.left(), atomResolver, leafShape),
						infer(c.right(), atomResolver, leafShape));
			}
			case Lineage.ConcatRows c -> {
				return concat(infer(c.top(), atomResolver, leafShape),
						infer(c.bottom(), atomResolver, leafShape));
			}
			case Lineage.SumInner s -> {
				return concat(infer(s.left(), atomResolver, leafShape),
						infer(s.right(), atomResolver, leafShape));
			}
			// Recombination: re-derive per-base-term sub-shapes from the base +
			// allocation, then lift each base-factor class by its leaf's classes
			// (g·u). Needs the base expanded (small atom) and a leaf-by-shape
			// resolver. Structural lower bound on buds (misses coincidental merges).
			case Lineage.RecombinationN r -> {
				NonCubicBilinearAlgorithm baseAlg =
						(r.base() instanceof Lineage.Atom a) ? atomResolver.apply(a.ref()) : null;
				if (baseAlg == null) return Profile.UNKNOWN;
				return inferRecombination(baseAlg, r.allocA(), r.allocB(), r.allocC(), leafShape);
			}
			case Lineage.RecombinationWithPairN ignored -> { return Profile.UNKNOWN; }
			case Lineage.AugmentSquareDiscard ignored -> { return Profile.UNKNOWN; }
			case Lineage.SerendipitousProduct ignored -> { return Profile.UNKNOWN; }
			case Lineage.Project ignored -> { return Profile.UNKNOWN; }
			case Lineage.DisjointSum ignored -> { return Profile.UNKNOWN; }
		}
	}

	// ── propagation rules ──

	/** Kronecker: per type, the class sizes of A⊗B are {a·b : a∈A, b∈B}. Exact. */
	private static Profile kron(Profile a, Profile b) {
		if (!a.known() || !b.known()) return Profile.UNKNOWN;
		return new Profile(a.rank() * b.rank(),
				cartesian(a.uClasses(), b.uClasses()),
				cartesian(a.vClasses(), b.vClasses()),
				cartesian(a.wClasses(), b.wClasses()),
				Certainty.min(a.certainty(), b.certainty()));
	}

	/**
	 * Recombination bud estimate. The composite's term {@code (kBase, kSub)} is
	 * base-term {@code kBase}'s block-spread combined with leaf-term {@code kSub};
	 * two result U-columns are proportional iff the base-U columns are (same
	 * base-U-class), the leaves are the same (same sub-shape), and the leaf-U
	 * columns are. So per type the result classes are, grouped by (base class,
	 * sub-shape) of size {@code g}: {@code g·u} for each leaf class {@code u}.
	 *
	 * <p>Sub-shapes are re-derived from {@code base + alloc} (independent of the
	 * SOTA resolver — a dummy {@code →0} is used). Exact given exact leaf
	 * profiles, except it cannot see <em>coincidental</em> cross-class merges, so
	 * it is a lower bound → {@code STRUCTURAL_ESTIMATE}.</p>
	 */
	static Profile inferRecombination(NonCubicBilinearAlgorithm baseAlg,
			int[] allocA, int[] allocB, int[] allocC, Function<int[], Profile> leafShape) {
		if (baseAlg == null || leafShape == null) return Profile.UNKNOWN;
		Recombination.Result rec;
		try {
			rec = Recombination.recombineWithAllocation(baseAlg, (n, m, p) -> 0,
					allocA, allocB, allocC);
		} catch (RuntimeException e) {
			return Profile.UNKNOWN;
		}
		int[][] sub = rec.smallMatrixSizes;
		int baseR = baseAlg.r;
		int[][] ids = SerendipitousBudProduct.independentClassIds(baseAlg);

		// Resolve each distinct non-degenerate sub-shape's leaf profile once.
		Map<String, Profile> leafCache = new HashMap<>();
		Certainty cert = Certainty.STRUCTURAL_ESTIMATE;
		int rank = 0;
		for (int k = 0; k < baseR; k++) {
			int[] sz = sub[k];
			if (sz[0] == 0 || sz[1] == 0 || sz[2] == 0) continue;  // peeled-away block
			String sk = sz[0] + "x" + sz[1] + "x" + sz[2];
			Profile lp = leafCache.get(sk);
			if (lp == null) {
				lp = leafShape.apply(sz);
				if (lp == null || !lp.known()) return Profile.UNKNOWN;
				leafCache.put(sk, lp);
			}
			cert = Certainty.min(cert, lp.certainty());
			rank += lp.rank();
		}
		return new Profile(rank,
				mergeRecombType(baseR, sub, ids[0], leafCache, 0),
				mergeRecombType(baseR, sub, ids[1], leafCache, 1),
				mergeRecombType(baseR, sub, ids[2], leafCache, 2),
				cert);
	}

	/** Per type: group base terms by (base-class, sub-shape) → emit {@code g·u}. */
	private static int[] mergeRecombType(int baseR, int[][] sub, int[] classId,
			Map<String, Profile> leafCache, int type) {
		Map<String, Integer> count = new LinkedHashMap<>();
		Map<String, String> shapeKeyOf = new HashMap<>();
		for (int k = 0; k < baseR; k++) {
			int[] sz = sub[k];
			if (sz[0] == 0 || sz[1] == 0 || sz[2] == 0) continue;
			String sk = sz[0] + "x" + sz[1] + "x" + sz[2];
			String key = classId[k] + "|" + sk;
			count.merge(key, 1, Integer::sum);
			shapeKeyOf.put(key, sk);
		}
		List<Integer> out = new ArrayList<>();
		for (Map.Entry<String, Integer> e : count.entrySet()) {
			Profile lp = leafCache.get(shapeKeyOf.get(e.getKey()));
			int[] classes = switch (type) {
				case 0 -> lp.uClasses();
				case 1 -> lp.vClasses();
				default -> lp.wClasses();
			};
			int g = e.getValue();
			for (int c : classes) out.add(g * c);
		}
		return descending(out.stream().mapToInt(Integer::intValue).toArray());
	}

	/** Direct-sum-like union of class multisets (structural lower bound). */
	private static Profile concat(Profile a, Profile b) {
		if (!a.known() || !b.known()) return Profile.UNKNOWN;
		return new Profile(a.rank() + b.rank(),
				union(a.uClasses(), b.uClasses()),
				union(a.vClasses(), b.vClasses()),
				union(a.wClasses(), b.wClasses()),
				Certainty.min(Certainty.STRUCTURAL_ESTIMATE,
						Certainty.min(a.certainty(), b.certainty())));
	}

	/**
	 * Tensor-symmetry relabel. {@code perm} is {@code "ABC->XYZ"}: the new
	 * A/B/C factor takes the partition of the OLD factor named by X/Y/Z. Rank
	 * and the partitions themselves are unchanged.
	 */
	private static Profile transpose(Profile c, String perm) {
		if (!c.known()) return Profile.UNKNOWN;
		int arrow = perm.indexOf("->");
		if (arrow < 0 || perm.length() < arrow + 5) {
			return c;  // unrecognised perm string — leave classes as-is
		}
		String dst = perm.substring(arrow + 2);
		int[] u = classesFor(c, dst.charAt(0));
		int[] v = classesFor(c, dst.charAt(1));
		int[] w = classesFor(c, dst.charAt(2));
		return new Profile(c.rank(), u, v, w, c.certainty());
	}

	private static int[] classesFor(Profile c, char letter) {
		return switch (letter) {
			case 'A' -> c.uClasses();
			case 'B' -> c.vClasses();
			case 'C' -> c.wClasses();
			default -> c.uClasses();
		};
	}

	private static int[] cartesian(int[] x, int[] y) {
		int[] out = new int[x.length * y.length];
		int k = 0;
		for (int a : x) for (int b : y) out[k++] = a * b;
		return descending(out);
	}

	private static int[] union(int[] x, int[] y) {
		int[] out = Arrays.copyOf(x, x.length + y.length);
		System.arraycopy(y, 0, out, x.length, y.length);
		return descending(out);
	}

	private static int[] descending(int[] a) {
		Arrays.sort(a);
		for (int i = 0, j = a.length - 1; i < j; i++, j--) {
			int t = a[i]; a[i] = a[j]; a[j] = t;
		}
		return a;
	}
}
