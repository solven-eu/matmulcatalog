package eu.solven.matmul.catalog;

import eu.solven.matmul.recombination.Recombination;

import java.util.ArrayList;
import java.util.List;

import eu.solven.matmul.NonCubicBilinearAlgorithm;

/**
 * Serendipitous product via BUD decomposition (#159). Originates as
 * <strong>Smith 2002 §9.3 "serendipitous equalities"</strong> (Warren D. Smith,
 * "Fast matrix multiplication formulae") — eq. (69):
 * {@code Rk⟨3a,3b,3c⟩ ≤ 19·Rk⟨a,b,c⟩ + 2·Rk⟨2a,b,c⟩}, from the
 * Johnson–McLoughlin {@code ⟨3,3,3⟩=23} scheme whose 2 proportional
 * coefficient-vector pairs (= buds) let 2 of the 23 blocks be realised as
 * doubled-axis {@code ⟨2a,b,c⟩} blocks instead (→ {@code ⟨9,9,9⟩≤527<529}).
 * Perminov's draft (Def 2.9–2.12, §2.6) and Sedoglavic generalise eq. (69) to
 * an ARBITRARY bud structure, which is what this class computes. A scheme
 * decomposes into elementary matmul tensors by its buds —
 * rank-one terms sharing a {@code u} (or {@code v}/{@code w}) vector up to
 * scaling. The product with a second scheme {@code ⟨n₂,m₂,p₂⟩} realizes each
 * enlarged elementary block at its best known rank:
 * {@code r_s = Σ Sᵢ·R(⟨Nᵢn₂,Mᵢm₂,Pᵢp₂⟩)}.
 *
 * <p>This is the SAME objective the recombination / block-split search solves —
 * {@code min Σ R(sub-shapes)} over a decomposition, against a SOTA rank oracle
 * {@code R(·)}. There the decomposition is a freely-chosen block ALLOCATION
 * (branch-and-bound searched → optimal-within-scope); here it is the base's BUD
 * STRUCTURE, which is non-unique (alternative groupings exist). We currently take
 * the deterministic GREEDY decomposition, so {@code r_s} and {@code bud_score}
 * are UPPER BOUNDS — an optimal-bud-structure search (the dual of the allocation
 * search) would tighten them. See {@code references/SERENDIPITOUS_PARTIAL_PRODUCT.md}
 * (the saving {@code k·rB − Σ R(...)} framing).
 *
 * <p>Covers single-type {@code U}-buds ({@code ⟨1,1,k⟩}), {@code V}-buds
 * ({@code ⟨k,1,1⟩}) and {@code W}-buds ({@code ⟨1,k,1⟩}) plus trivial terms.
 * Combined buds (§2.6.4, e.g. {@code ⟨2,1,2⟩}) are a follow-up. The assembled
 * scheme is a standard flatten {@code (U,V,W)}; the caller verifies it with
 * {@link eu.solven.matmul.verifiers.Verifier#isExactNonCubic} (the oracle).</p>
 */
public final class SerendipitousBudProduct {

	private SerendipitousBudProduct() {}

	public enum BudType { U, V, W }

	/** A bud: a group of ≥2 term indices sharing one factor vector (up to scaling). */
	public record Bud(BudType type, int[] terms) {}

	/** Full bud decomposition: typed buds + leftover trivial terms. */
	public record BudDecomposition(List<Bud> buds, int[] trivial) {
		/** U-buds only (back-compat with the original probe). */
		public List<int[]> uBuds() {
			List<int[]> out = new ArrayList<>();
			for (Bud b : buds) if (b.type() == BudType.U) out.add(b.terms());
			return out;
		}
	}

	/** Compact, JSON-friendly summary of a scheme's (greedy U→V→W) bud structure. */
	public record BudSummary(boolean hasBuds, int uBuds, int vBuds, int wBuds, int trivial,
			String summary) {}

	/**
	 * Summarise the bud structure: one canonical greedy decomposition rendered as
	 * e.g. {@code "4×U⟨1,1,2⟩ + 12×⟨1,1,1⟩"}. Cheap (O(r²·dim)); intended for the
	 * catalog manifest / per-scheme JSON. NOT unique — alternative groupings exist
	 * (see {@code references/SERENDIPITOUS_PARTIAL_PRODUCT.md}); this is the
	 * deterministic greedy one.
	 */
	public static BudSummary summarise(NonCubicBilinearAlgorithm a) {
		BudDecomposition dec = findBuds(a);
		java.util.TreeMap<String, Integer> groups = new java.util.TreeMap<>();
		int u = 0, v = 0, w = 0;
		for (Bud b : dec.buds()) {
			int k = b.terms().length;
			String tag = switch (b.type()) {
				case U -> "U⟨1,1," + k + "⟩";
				case V -> "V⟨" + k + ",1,1⟩";
				case W -> "W⟨1," + k + ",1⟩";
			};
			groups.merge(tag, 1, Integer::sum);
			switch (b.type()) { case U -> u++; case V -> v++; case W -> w++; }
		}
		StringBuilder sb = new StringBuilder();
		groups.forEach((shape, c) -> {
			if (sb.length() > 0) sb.append(" + ");
			sb.append(c).append("×").append(shape);
		});
		if (dec.trivial().length > 0) {
			if (sb.length() > 0) sb.append(" + ");
			sb.append(dec.trivial().length).append("×⟨1,1,1⟩");
		}
		return new BudSummary(!dec.buds().isEmpty(), u, v, w, dec.trivial().length, sb.toString());
	}

	/** Default greedy type ordering for {@link #findBuds(NonCubicBilinearAlgorithm)}. */
	public static final BudType[] DEFAULT_ORDER = { BudType.U, BudType.V, BudType.W };

	/**
	 * All 6 bud-type orderings. The greedy is <strong>order-sensitive</strong> —
	 * a term that belongs to both a U-class and a V-class is consumed by whichever
	 * type is processed first — so the U→V→W default can hide a larger, cheaper bud
	 * of a later type (e.g. a size-3 V-bud masked into a size-2 U-bud + leftovers,
	 * which is exactly the ⟨8,9,9⟩=430 = ⟨4,3,3⟩⊗⟨2,3,3⟩ + ⟨6,3,3⟩ case). Callers
	 * after the cheapest serendipitous cost should try them all and keep the min.
	 */
	public static final BudType[][] ALL_ORDERINGS = {
			{ BudType.U, BudType.V, BudType.W }, { BudType.U, BudType.W, BudType.V },
			{ BudType.V, BudType.U, BudType.W }, { BudType.V, BudType.W, BudType.U },
			{ BudType.W, BudType.U, BudType.V }, { BudType.W, BudType.V, BudType.U } };

	/** Greedy U→V→W decomposition (back-compat default). */
	public static BudDecomposition findBuds(NonCubicBilinearAlgorithm a) {
		return findBuds(a, DEFAULT_ORDER);
	}

	/** Greedy bud decomposition under an explicit type ordering (see {@link #ALL_ORDERINGS}). */
	public static BudDecomposition findBuds(NonCubicBilinearAlgorithm a, BudType[] order) {
		double[][] srcU = a.denseU();
		double[][] srcV = a.denseV();
		double[][] srcW = a.denseW();
		boolean[] used = new boolean[a.r];
		List<Bud> buds = new ArrayList<>();
		for (BudType t : order) {
			double[][] src = switch (t) { case U -> srcU; case V -> srcV; case W -> srcW; };
			groupBy(a, src, t, used, buds);
		}
		List<Integer> trivial = new ArrayList<>();
		for (int l = 0; l < a.r; l++) if (!used[l]) trivial.add(l);
		return new BudDecomposition(buds, trivial.stream().mapToInt(Integer::intValue).toArray());
	}

	/**
	 * Independent per-factor class-size distributions: for each of U, V, W
	 * <em>separately</em>, partition the r terms by proportional direction and
	 * return the multiset of class sizes (descending, summing to r). Unlike the
	 * greedy disjoint {@link #findBuds} — which assigns each term to at most one
	 * bud (U first, then V, then W) — these are three independent partitions, so
	 * one term may simultaneously belong to a U-class, a V-class and a W-class.
	 *
	 * <p>This is the <strong>composition-stable</strong> notion of bud structure
	 * used by {@link LineageBudInference}: it is the one that propagates exactly
	 * through Kronecker (Cartesian product of class sizes) and tensor-symmetry
	 * relabelling. It generally reports more/larger buds than the human-display
	 * greedy {@link #summarise}. Returns {@code {uClassSizes, vClassSizes,
	 * wClassSizes}}.</p>
	 */
	public static int[][] independentClassSizes(NonCubicBilinearAlgorithm a) {
		double[][] srcU = a.denseU();
		double[][] srcV = a.denseV();
		double[][] srcW = a.denseW();
		return new int[][] { classSizes(srcU, a.r), classSizes(srcV, a.r), classSizes(srcW, a.r) };
	}

	/**
	 * Independent per-factor class <em>IDs</em>: {@code {uIds, vIds, wIds}}, each
	 * length {@code r}, where {@code uIds[l]} is the U-direction class index of
	 * term {@code l} (terms with proportional U columns share an id). Companion to
	 * {@link #independentClassSizes} used by recombination bud-inference to group
	 * base terms by their base-factor class. IDs are assigned in first-seen order.
	 */
	public static int[][] independentClassIds(NonCubicBilinearAlgorithm a) {
		double[][] srcU = a.denseU();
		double[][] srcV = a.denseV();
		double[][] srcW = a.denseW();
		return new int[][] { classIds(srcU, a.r), classIds(srcV, a.r), classIds(srcW, a.r) };
	}

	private static int[] classIds(double[][] factor, int r) {
		java.util.LinkedHashMap<String, Integer> idOf = new java.util.LinkedHashMap<>();
		int[] ids = new int[r];
		for (int l = 0; l < r; l++) {
			String key = java.util.Arrays.toString(canonicalDirection(column(factor, l)));
			ids[l] = idOf.computeIfAbsent(key, k -> idOf.size());
		}
		return ids;
	}

	private static int[] classSizes(double[][] factor, int r) {
		java.util.LinkedHashMap<String, Integer> byDir = new java.util.LinkedHashMap<>();
		for (int l = 0; l < r; l++) {
			byDir.merge(java.util.Arrays.toString(canonicalDirection(column(factor, l))), 1, Integer::sum);
		}
		int[] sizes = byDir.values().stream().mapToInt(Integer::intValue).toArray();
		java.util.Arrays.sort(sizes);
		for (int i = 0, j = sizes.length - 1; i < j; i++, j--) {
			int t = sizes[i]; sizes[i] = sizes[j]; sizes[j] = t;  // descending
		}
		return sizes;
	}

	/** Back-compat: U-buds + trivial. */
	public static BudDecomposition findUBuds(NonCubicBilinearAlgorithm a) {
		double[][] srcU = a.denseU();
		boolean[] used = new boolean[a.r];
		List<Bud> buds = new ArrayList<>();
		groupBy(a, srcU, BudType.U, used, buds);
		List<Integer> trivial = new ArrayList<>();
		for (int l = 0; l < a.r; l++) if (!used[l]) trivial.add(l);
		return new BudDecomposition(buds, trivial.stream().mapToInt(Integer::intValue).toArray());
	}

	private static void groupBy(NonCubicBilinearAlgorithm a, double[][] factor, BudType type,
			boolean[] used, List<Bud> buds) {
		// Hash by canonical direction → O(r·dim) instead of O(r²·dim).
		java.util.LinkedHashMap<String, List<Integer>> byDir = new java.util.LinkedHashMap<>();
		for (int l = 0; l < a.r; l++) {
			if (used[l]) continue;
			String key = java.util.Arrays.toString(canonicalDirection(column(factor, l)));
			byDir.computeIfAbsent(key, k -> new ArrayList<>()).add(l);
		}
		for (List<Integer> grp : byDir.values()) {
			if (grp.size() >= 2) {
				int[] terms = grp.stream().mapToInt(Integer::intValue).toArray();
				for (int t : terms) used[t] = true;
				buds.add(new Bud(type, terms));
			}
		}
	}

	/**
	 * Predicted rank of the serendipitous product {@code T1 ⊗ ⟨n2,m2,p2⟩}
	 * <em>without building it</em> — each bud of size {@code k} contributes the
	 * rank of its enlarged inner block ({@code U}: grows p, {@code V}: grows n,
	 * {@code W}: grows m), and each lone term contributes {@code R(⟨n2,m2,p2⟩)}.
	 * Returns {@link Long#MAX_VALUE}/4 if any required enlarged inner rank is
	 * unknown (so callers treat it as "not computable"). This is the cost
	 * primitive for the bud-base factory.
	 */
	public static long serendipitousCost(NonCubicBilinearAlgorithm t1, FieldAwareLookup lookup,
			int n2, int m2, int p2) {
		long best = Long.MAX_VALUE / 4;
		for (BudType[] order : ALL_ORDERINGS) {
			best = Math.min(best, costOf(findBuds(t1, order), lookup, n2, m2, p2));
		}
		return best;
	}

	/** Predicted cost of a specific decomposition; {@code Long.MAX_VALUE/4} if any
	 *  enlarged inner rank is unknown. */
	public static long costOf(BudDecomposition dec, FieldAwareLookup lookup, int n2, int m2, int p2) {
		final long UNKNOWN = Long.MAX_VALUE / 4;
		long inner = lookup.findRank(n2, m2, p2);
		if (inner >= Recombination.SotaResolver.UNKNOWN_RANK) return UNKNOWN;
		long cost = (long) dec.trivial().length * inner;
		for (Bud b : dec.buds()) {
			int k = b.terms().length;
			long r = switch (b.type()) {
				case U -> lookup.findRank(n2, m2, k * p2);
				case V -> lookup.findRank(k * n2, m2, p2);
				case W -> lookup.findRank(n2, k * m2, p2);
			};
			if (r >= Recombination.SotaResolver.UNKNOWN_RANK) return UNKNOWN;
			cost += r;
		}
		return cost;
	}

	/** Serendipitous product {@code T1 ⊗ ⟨n2,m2,p2⟩} using all bud types. */
	public static NonCubicBilinearAlgorithm productViaBuds(
			NonCubicBilinearAlgorithm t1, FieldAwareLookup lookup, int n2, int m2, int p2) {
		return productViaBudsTyped(t1, lookup, n2, m2, p2, java.util.EnumSet.allOf(BudType.class));
	}

	/**
	 * Serendipitous product {@code T1 ⊗ ⟨n2,m2,p2⟩} built under the CHEAPEST bud-type
	 * ordering — the {@code min} over {@link #ALL_ORDERINGS}, exactly as
	 * {@code SerendipitousSearch.bestFor} selects. The greedy decomposition is
	 * order-sensitive (a term shared between a U- and a V-class goes to whichever type
	 * is processed first), so the {@link #DEFAULT_ORDER} that {@link #productViaBuds}
	 * uses can build at a higher rank than the search predicted. This is what made the
	 * {@code SerendipitousProduct} replay drift off the search's rank (e.g. ⟨14,16,25⟩
	 * predicted 3297, replayed 3310) → the write-guard discarded the win. Trying every
	 * ordering and keeping the lowest-rank build makes replay reproduce the search's
	 * choice deterministically. [[lineage replay must be bit-exact]]
	 */
	public static NonCubicBilinearAlgorithm productViaBudsBest(
			NonCubicBilinearAlgorithm t1, FieldAwareLookup lookup, int n2, int m2, int p2) {
		return productViaBudsBest(t1, InnerResolver.of(lookup), n2, m2, p2);
	}

	/** {@link #productViaBudsBest} against an explicit {@link InnerResolver} —
	 *  stub-capable callers (LineageReplayer, RecursiveMaterialiser) pass a
	 *  replaying resolver so stub-only fusion targets stay buildable. */
	public static NonCubicBilinearAlgorithm productViaBudsBest(
			NonCubicBilinearAlgorithm t1, InnerResolver resolver, int n2, int m2, int p2) {
		java.util.Set<BudType> allow = java.util.EnumSet.allOf(BudType.class);
		NonCubicBilinearAlgorithm best = null;
		for (BudType[] order : ALL_ORDERINGS) {
			BudDecomposition dec = findBuds(t1, order);
			if (dec.buds().isEmpty()) {
				continue; // no buds under this ordering → nothing to beat the others with
			}
			try {
				NonCubicBilinearAlgorithm built =
						productFromDecomposition(t1, dec, resolver, n2, m2, p2, allow);
				if (best == null || built.r < best.r) {
					best = built;
				}
			} catch (RuntimeException e) {
				// This ordering's fusion target is unavailable — skip to the next.
			}
		}
		// No ordering yielded a buildable bud decomposition → default (term-by-term Kron).
		return best != null ? best : productFromDecomposition(
				t1, findBuds(t1), resolver, n2, m2, p2, java.util.EnumSet.noneOf(BudType.class));
	}

	/**
	 * Serendipitous product fusing only the bud types in {@code allow}; buds of a
	 * disallowed type are realised term-by-term (plain Kronecker, no fusion).
	 * Used to bisect which bud-block construction is correct (verification debug)
	 * and to fall back when a type's fusion target is unknown/unbeneficial.
	 */
	public static NonCubicBilinearAlgorithm productViaBudsTyped(
			NonCubicBilinearAlgorithm t1, FieldAwareLookup lookup, int n2, int m2, int p2,
			java.util.Set<BudType> allow) {
		return productFromDecomposition(t1, findBuds(t1), lookup, n2, m2, p2, allow);
	}

	/**
	 * Resolves an explicit (buildable, ORIENTED) scheme at {@code ⟨n,m,p⟩} — the
	 * build-time ingredient fetch for serendipitous products. The default
	 * ({@link #of}) is {@code findWithSource}, which SKIPS lineage-only stubs; a
	 * stub-capable caller (RecursiveMaterialiser, LineageReplayer) passes a
	 * replaying resolver so a stub-only enlarged shape (e.g. the ⟨4,4,20⟩=230
	 * ConcatCols stub that prices ⟨20,28,28⟩=8434) is still buildable. Without
	 * this hook the search silently dropped every candidate whose fusion target
	 * had no dense file — predict saw the stub's rank via findRank, build threw.
	 */
	public interface InnerResolver {
		java.util.Optional<NonCubicBilinearAlgorithm> find(int n, int m, int p);

		/** Dense-file-only resolver (no stub replay) — the historical behaviour. */
		static InnerResolver of(FieldAwareLookup lookup) {
			return (n, m, p) -> lookup.findWithSource(n, m, p).map(FieldAwareLookup.WithSource::alg);
		}
	}

	/**
	 * Build the serendipitous product from a <strong>precomputed</strong>
	 * decomposition, so the built scheme matches the ordering whose cost was
	 * predicted (the greedy ordering changes which buds are chosen).
	 */
	public static NonCubicBilinearAlgorithm productFromDecomposition(
			NonCubicBilinearAlgorithm t1, BudDecomposition dec, FieldAwareLookup lookup,
			int n2, int m2, int p2, java.util.Set<BudType> allow) {
		return productFromDecomposition(t1, dec, InnerResolver.of(lookup), n2, m2, p2, allow);
	}

	/** {@link #productFromDecomposition} against an explicit {@link InnerResolver}. */
	public static NonCubicBilinearAlgorithm productFromDecomposition(
			NonCubicBilinearAlgorithm t1, BudDecomposition dec, InnerResolver resolver,
			int n2, int m2, int p2, java.util.Set<BudType> allow) {
		NonCubicBilinearAlgorithm s2 = resolver.find(n2, m2, p2).orElseThrow(
				() -> new IllegalStateException("no buildable ⟨" + n2 + "," + m2 + "," + p2 + "⟩"));
		List<NonCubicBilinearAlgorithm> parts = new ArrayList<>();
		for (int i : dec.trivial()) parts.add(Compose.kroneckerGeneral(rankOne(t1, i), s2));
		for (Bud bud : dec.buds()) {
			if (!allow.contains(bud.type())) {
				// no fusion: realise each bud term as a plain ⟨n2,m2,p2⟩ copy.
				for (int term : bud.terms()) {
					parts.add(Compose.kroneckerGeneral(rankOne(t1, term), s2));
				}
				continue;
			}
			int k = bud.terms().length;
			int en = bud.type() == BudType.V ? k * n2 : n2;
			int em = bud.type() == BudType.W ? k * m2 : m2;
			int ep = bud.type() == BudType.U ? k * p2 : p2;
			NonCubicBilinearAlgorithm s3 = resolver.find(en, em, ep).orElseThrow(
					() -> new IllegalStateException("no buildable enlarged ⟨" + en + "," + em + ","
							+ ep + "⟩ (stub-only? pass a replaying InnerResolver)"));
			parts.add(buildBudBlock(t1, bud, s3, n2, m2, p2));
		}
		return concatColumns(parts, t1.n * n2, t1.m * m2, t1.p * p2);
	}

	/** Back-compat alias (U-buds + trivial path is subsumed by productViaBuds). */
	public static NonCubicBilinearAlgorithm productViaUBuds(
			NonCubicBilinearAlgorithm t1, FieldAwareLookup lookup, int n2, int m2, int p2) {
		return productViaBuds(t1, lookup, n2, m2, p2);
	}

	// ── bud-block construction (Perminov §2.6.2–2.6.3), our index convention ──
	private static NonCubicBilinearAlgorithm buildBudBlock(
			NonCubicBilinearAlgorithm t1, Bud bud, NonCubicBilinearAlgorithm s3,
			int n2, int m2, int p2) {
		double[][] srcU = t1.denseU();
		double[][] srcV = t1.denseV();
		double[][] srcW = t1.denseW();
		int n1 = t1.n, m1 = t1.m, p1 = t1.p, k = bud.terms().length;
		int N = n1 * n2, M = m1 * m2, P = p1 * p2, r3 = s3.r;
		double[][] U = new double[N * M][r3], V = new double[M * P][r3], W = new double[N * P][r3];

		// Shared vector = the bud's common factor (rescale terms so it is exact).
		double[][] shared = switch (bud.type()) { case U -> srcU; case V -> srcV; case W -> srcW; };
		double[] sbar = column(shared, bud.terms()[0]);
		double[] scale = new double[k];
		for (int l = 0; l < k; l++) scale[l] = proportionFactor(sbar, column(shared, bud.terms()[l]));

		for (int j = 0; j < r3; j++) {
			for (int l = 0; l < k; l++) {
				int term = bud.terms()[l];
				double sc = scale[l];
				double[] uL = column(srcU, term), vL = column(srcV, term), wL = column(srcW, term);
				switch (bud.type()) {
					case U -> { // shared u; S3=⟨n2,m2,k·p2⟩; split p; sum v,w; scale→v
						if (l == 0) kronU(U, sbar, n1, m1, s3, j, n2, m2, M, /*off*/ 0, /*span*/ m2, true);
						addV(V, mul(vL, sc), m1, p1, s3, j, m2, p2, k, l, P, BudType.U);
						addW(W, wL, n1, p1, s3, j, n2, p2, k, l, P, BudType.U);
					}
					case V -> { // shared v; S3=⟨k·n2,m2,p2⟩; split n; sum u,w; scale→u
						if (l == 0) kronV(V, sbar, m1, p1, s3, j, m2, p2, P);
						addU(U, mul(uL, sc), n1, m1, s3, j, n2, m2, k, l, M, BudType.V);
						addW(W, wL, n1, p1, s3, j, n2, p2, k, l, P, BudType.V);
					}
					case W -> { // shared w; S3=⟨n2,k·m2,p2⟩; split m; sum u,v; scale→u
						if (l == 0) kronW(W, sbar, n1, p1, s3, j, n2, p2, P);
						addU(U, mul(uL, sc), n1, m1, s3, j, n2, m2, k, l, M, BudType.W);
						addV(V, vL, m1, p1, s3, j, m2, p2, k, l, P, BudType.W);
					}
				}
			}
		}
		return new NonCubicBilinearAlgorithm(N, M, P, U, V, W);
	}

	// shared-factor Kronecker (the common vector ⊗ S3's column, full not split)
	private static void kronU(double[][] U, double[] ubar, int n1, int m1, NonCubicBilinearAlgorithm s3,
			int j, int n2, int m2, int M, int off, int span, boolean ignore) {
		double[][] srcU = s3.denseU();
		for (int i1 = 0; i1 < n1; i1++) for (int j1 = 0; j1 < m1; j1++) {
			double uu = ubar[i1 * m1 + j1]; if (uu == 0) continue;
			for (int i2 = 0; i2 < n2; i2++) for (int j2 = 0; j2 < m2; j2++) {
				double v = uu * srcU[i2 * m2 + j2][j];
				if (v != 0) U[(i1 * n2 + i2) * M + (j1 * m2 + j2)][j] = v;
			}
		}
	}

	private static void kronV(double[][] V, double[] vbar, int m1, int p1, NonCubicBilinearAlgorithm s3,
			int j, int m2, int p2, int P) {
		double[][] srcV = s3.denseV();
		for (int j1 = 0; j1 < m1; j1++) for (int k1 = 0; k1 < p1; k1++) {
			double vv = vbar[j1 * p1 + k1]; if (vv == 0) continue;
			for (int j2 = 0; j2 < m2; j2++) for (int k2 = 0; k2 < p2; k2++) {
				double v = vv * srcV[j2 * p2 + k2][j];
				if (v != 0) V[(j1 * m2 + j2) * P + (k1 * p2 + k2)][j] = v;
			}
		}
	}

	private static void kronW(double[][] W, double[] wbar, int n1, int p1, NonCubicBilinearAlgorithm s3,
			int j, int n2, int p2, int P) {
		double[][] srcW = s3.denseW();
		for (int i1 = 0; i1 < n1; i1++) for (int k1 = 0; k1 < p1; k1++) {
			double ww = wbar[i1 * p1 + k1]; if (ww == 0) continue;
			for (int i2 = 0; i2 < n2; i2++) for (int k2 = 0; k2 < p2; k2++) {
				double v = ww * srcW[i2 * p2 + k2][j];
				if (v != 0) W[(i1 * n2 + i2) * P + (k1 * p2 + k2)][j] = v;
			}
		}
	}

	// summed-factor contributions: per-bud-term ⊗ the l-th split block of S3
	private static void addU(double[][] U, double[] uL, int n1, int m1, NonCubicBilinearAlgorithm s3,
			int j, int n2, int m2, int k, int l, int M, BudType type) {
		// U index in S3 depends on which axis is split: V-bud splits n (S3 U is (k·n2)×m2);
		// W-bud splits m (S3 U is n2×(k·m2)).
		double[][] srcU = s3.denseU();
		for (int i1 = 0; i1 < n1; i1++) for (int j1 = 0; j1 < m1; j1++) {
			double uu = uL[i1 * m1 + j1]; if (uu == 0) continue;
			for (int i2 = 0; i2 < n2; i2++) for (int j2 = 0; j2 < m2; j2++) {
				double s = (type == BudType.V)
						? srcU[(l * n2 + i2) * m2 + j2][j]        // split n
						: srcU[i2 * (k * m2) + (l * m2 + j2)][j]; // split m (W-bud)
				if (s != 0) U[(i1 * n2 + i2) * M + (j1 * m2 + j2)][j] += uu * s;
			}
		}
	}

	private static void addV(double[][] V, double[] vL, int m1, int p1, NonCubicBilinearAlgorithm s3,
			int j, int m2, int p2, int k, int l, int P, BudType type) {
		// V index: U-bud splits p (S3 V is m2×(k·p2)); W-bud splits m (S3 V is (k·m2)×p2).
		double[][] srcV = s3.denseV();
		for (int j1 = 0; j1 < m1; j1++) for (int k1 = 0; k1 < p1; k1++) {
			double vv = vL[j1 * p1 + k1]; if (vv == 0) continue;
			for (int j2 = 0; j2 < m2; j2++) for (int k2 = 0; k2 < p2; k2++) {
				double s = (type == BudType.U)
						? srcV[j2 * (k * p2) + (l * p2 + k2)][j]  // split p
						: srcV[(l * m2 + j2) * p2 + k2][j];       // split m (W-bud)
				if (s != 0) V[(j1 * m2 + j2) * P + (k1 * p2 + k2)][j] += vv * s;
			}
		}
	}

	private static void addW(double[][] W, double[] wL, int n1, int p1, NonCubicBilinearAlgorithm s3,
			int j, int n2, int p2, int k, int l, int P, BudType type) {
		// W index: U-bud splits p (S3 W is n2×(k·p2)); V-bud splits n (S3 W is (k·n2)×p2).
		double[][] srcW = s3.denseW();
		for (int i1 = 0; i1 < n1; i1++) for (int k1 = 0; k1 < p1; k1++) {
			double ww = wL[i1 * p1 + k1]; if (ww == 0) continue;
			for (int i2 = 0; i2 < n2; i2++) for (int k2 = 0; k2 < p2; k2++) {
				double s = (type == BudType.U)
						? srcW[i2 * (k * p2) + (l * p2 + k2)][j]  // split p
						: srcW[(l * n2 + i2) * p2 + k2][j];       // split n (V-bud)
				if (s != 0) W[(i1 * n2 + i2) * P + (k1 * p2 + k2)][j] += ww * s;
			}
		}
	}

	// ── small helpers ──
	private static double[] mul(double[] v, double s) {
		double[] o = new double[v.length];
		for (int i = 0; i < v.length; i++) o[i] = v[i] * s;
		return o;
	}

	private static NonCubicBilinearAlgorithm rankOne(NonCubicBilinearAlgorithm a, int col) {
		double[][] srcU = a.denseU();
		double[][] srcV = a.denseV();
		double[][] srcW = a.denseW();
		double[][] U = new double[a.dimU()][1], V = new double[a.dimV()][1], W = new double[a.dimW()][1];
		for (int i = 0; i < a.dimU(); i++) U[i][0] = srcU[i][col];
		for (int i = 0; i < a.dimV(); i++) V[i][0] = srcV[i][col];
		for (int i = 0; i < a.dimW(); i++) W[i][0] = srcW[i][col];
		return new NonCubicBilinearAlgorithm(a.n, a.m, a.p, U, V, W);
	}

	private static NonCubicBilinearAlgorithm concatColumns(
			List<NonCubicBilinearAlgorithm> parts, int N, int M, int P) {
		int r = 0;
		for (NonCubicBilinearAlgorithm p : parts) r += p.r;
		double[][] U = new double[N * M][r], V = new double[M * P][r], W = new double[N * P][r];
		int off = 0;
		for (NonCubicBilinearAlgorithm p : parts) {
			double[][] srcU = p.denseU();
			double[][] srcV = p.denseV();
			double[][] srcW = p.denseW();
			for (int row = 0; row < N * M; row++) for (int c = 0; c < p.r; c++) U[row][off + c] = srcU[row][c];
			for (int row = 0; row < M * P; row++) for (int c = 0; c < p.r; c++) V[row][off + c] = srcV[row][c];
			for (int row = 0; row < N * P; row++) for (int c = 0; c < p.r; c++) W[row][off + c] = srcW[row][c];
			off += p.r;
		}
		return new NonCubicBilinearAlgorithm(N, M, P, U, V, W);
	}

	private static double[] column(double[][] M, int col) {
		double[] c = new double[M.length];
		for (int i = 0; i < M.length; i++) c[i] = M[i][col];
		return c;
	}

	private static int[] canonicalDirection(double[] v) {
		double s = 0;
		for (double x : v) if (x != 0) { s = x; break; }
		int[] key = new int[v.length];
		if (s == 0) return key;
		for (int i = 0; i < v.length; i++) key[i] = (int) Math.round(v[i] / s * 1_000_000.0);
		return key;
	}

	private static double proportionFactor(double[] base, double[] target) {
		for (int i = 0; i < base.length; i++) if (base[i] != 0) return target[i] / base[i];
		return 1.0;
	}
}
