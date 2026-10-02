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
 * ({@code ⟨k,1,1⟩}) and {@code W}-buds ({@code ⟨1,k,1⟩}), trivial terms, and
 * <b>combined two-axis buds</b> (§2.6.4 — {@link GridBud}: {@code ⟨1,b,c⟩},
 * {@code ⟨a,1,c⟩}, {@code ⟨a,b,1⟩}, e.g. the {@code ⟨1,2,2⟩} of a ⟨3,3,3⟩:23 that
 * fuses four ⟨6,3,3⟩ copies into one ⟨6,6,6⟩:153 → ⟨9,9,18⟩ ≤ 913). Three-axis
 * buds ({@code a,b,c} all ≥ 2) are not searched. The assembled scheme is a
 * standard flatten {@code (U,V,W)}; the caller verifies it with
 * {@link eu.solven.matmul.verifiers.Verifier#isExactNonCubic} (the oracle).</p>
 *
 * <p><b>Optimality tier: bound.</b> The decomposition is the cheapest of a fixed
 * structural family ({@link #candidateDecompositions}: grid strategies × type
 * orderings, each greedy) — not an optimum over all bud partitions.</p>
 */
public final class SerendipitousBudProduct {

	private SerendipitousBudProduct() {}

	public enum BudType { U, V, W }

	/** A bud: a group of ≥2 term indices sharing one factor vector (up to scaling). */
	public record Bud(BudType type, int[] terms) {}

	/**
	 * A combined bud (Perminov §2.6.4): terms that together form an embedded
	 * {@code ⟨a,b,c⟩} matmul tensor with exactly one of {@code a,b,c} equal to 1 —
	 * a full grid over the proportionality classes of two factors:
	 * <ul>
	 *   <li>{@code ⟨1,b,c⟩} — {@code b} U-classes × {@code c} W-classes (V free);</li>
	 *   <li>{@code ⟨a,1,c⟩} — {@code a} U-classes × {@code c} V-classes (W free);</li>
	 *   <li>{@code ⟨a,b,1⟩} — {@code a} W-classes × {@code b} V-classes (U free).</li>
	 * </ul>
	 * {@code terms} is row-major over {@code (i,j,k)}: index {@code (i·b + j)·c + k}.
	 * Against an inner {@code ⟨n₂,m₂,p₂⟩} the {@code a·b·c} copies fuse into ONE
	 * {@code ⟨a·n₂, b·m₂, c·p₂⟩}. The single-type buds are the degenerate grids with
	 * two unit dims; a grid beats its own row/column split whenever
	 * {@code R(⟨an₂,bm₂,cp₂⟩)} is below the sum of the single-axis fusions.
	 *
	 * <p>No scale condition is needed: each term is
	 * {@code (α·Û) ⊗ V ⊗ (γ·Ŵ)} for its two class representatives, and
	 * {@code α·γ} is absorbed by the free factor.</p>
	 */
	public record GridBud(int a, int b, int c, int[] terms) {
		public GridBud {
			if ((a == 1 ? 1 : 0) + (b == 1 ? 1 : 0) + (c == 1 ? 1 : 0) != 1 || terms.length != a * b * c) {
				throw new IllegalArgumentException("a combined bud is ⟨1,b,c⟩ / ⟨a,1,c⟩ / ⟨a,b,1⟩ with a·b·c terms, got ⟨"
						+ a + "," + b + "," + c + "⟩ with " + terms.length);
			}
		}
	}

	/** Which combined buds a decomposition extracts before the single-type greedy. */
	public enum GridStrategy {
		/** None — single-type buds only (the pre-2026-09-30 behaviour). */
		NONE,
		/** Every grid type, largest first. */
		ALL,
		/**
		 * Every grid type, but only CLOSED grids: each row-class and column-class lies
		 * entirely inside the grid, so fusing it strands no class-mate.
		 */
		CLOSED,
		/**
		 * Every grid type, but no grid may SPLIT a U-class (each U-class it touches lies
		 * entirely inside it). Where a grid would take two of a three-term U-class, the
		 * class is left whole for the single-type greedy — ⟨3,10,14⟩:312 ⊗ˢ ⟨3,3,2⟩ keeps
		 * its four {@code ⟨1,1,3⟩} next to 69 {@code ⟨1,2,2⟩} (4244) instead of 71 grids
		 * and four broken triples (4252).
		 */
		KEEP_U,
		/** As {@link #KEEP_U}, for V-classes. */
		KEEP_V,
		/** As {@link #KEEP_U}, for W-classes. */
		KEEP_W,
		/**
		 * Every grid type, each grid cut down to its 2×2 CORE: the two rows and two
		 * columns whose four terms are least entangled with the free factor (fewest
		 * terms that also sit in a class of the factor the grid leaves free). A 2×5 grid
		 * built from the third members of eight U-triples costs those triples more than
		 * it saves; its core over the two triple-free columns does not —
		 * ⟨2,5,13⟩:104 ⊗ˢ ⟨2,5,2⟩ = 1657 with one {@code ⟨2,2,1⟩}, 1659 without.
		 */
		CORE,
		/** Only {@code ⟨1,b,c⟩} grids (U-classes × W-classes). */
		UW,
		/** Only {@code ⟨a,1,c⟩} grids (U-classes × V-classes). */
		UV,
		/** Only {@code ⟨a,b,1⟩} grids (W-classes × V-classes). */
		VW
	}

	/** Rank oracle for pricing a decomposition; {@code ≥ UNKNOWN_RANK} = unknown. */
	@FunctionalInterface
	public interface RankOracle {
		long rank(int n, int m, int p);

		static RankOracle of(FieldAwareLookup lookup) {
			return lookup::findRank;
		}

		/** Prices by what the resolver can actually BUILD (the replay-side oracle). */
		static RankOracle of(InnerResolver resolver) {
			java.util.Map<String, Long> memo = new java.util.HashMap<>();
			return (n, m, p) -> {
				String key = n + "x" + m + "x" + p;
				Long known = memo.get(key);
				if (known != null) {
					return known;
				}
				long r;
				try {
					r = resolver.find(n, m, p).map(a -> (long) a.r).orElse(Long.MAX_VALUE / 4);
				} catch (RuntimeException e) {
					r = Long.MAX_VALUE / 4;
				}
				memo.put(key, r);
				return r;
			};
		}
	}

	/** Full bud decomposition: typed buds + combined buds + leftover trivial terms. */
	public record BudDecomposition(List<Bud> buds, List<GridBud> grids, int[] trivial) {
		public BudDecomposition(List<Bud> buds, int[] trivial) {
			this(buds, List.of(), trivial);
		}

		/** Any fusion at all (single-type or combined)? */
		public boolean hasBuds() {
			return !buds.isEmpty() || !grids.isEmpty();
		}

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
		return decompose(a.r, independentClassIds(a), List.of(), order);
	}

	/**
	 * The structural family of decompositions a serendipitous product is chosen from:
	 * every {@link GridStrategy} that extracts at least one combined bud (plus
	 * {@link GridStrategy#NONE}, first — its six greedy members are the whole
	 * pre-2026-09-30 family) × every type ordering of {@link #ALL_ORDERINGS} for the
	 * terms left over, each ordering both as the class-at-a-time greedy and as the
	 * fewest-options-first greedy. It depends on the base ALONE — not on any rank oracle — so
	 * the search ({@code SerendipitousSearch.bestFor}, pricing with {@code findRank})
	 * and the replay ({@link #productViaBudsBest}, pricing with what it can build)
	 * minimise over the same set, and a better inner can only lower a replayed rank.
	 */
	public static List<BudDecomposition> candidateDecompositions(NonCubicBilinearAlgorithm a) {
		int[][] ids = independentClassIds(a);
		List<BudDecomposition> out = new ArrayList<>();
		java.util.Set<String> seenGridSets = new java.util.HashSet<>();
		for (GridStrategy strategy : GridStrategy.values()) {
			List<GridBud> grids = selectGrids(a.r, ids, strategy);
			if (strategy != GridStrategy.NONE && grids.isEmpty()) {
				continue;
			}
			StringBuilder sig = new StringBuilder();
			for (GridBud g : grids) {
				sig.append(java.util.Arrays.toString(g.terms())).append('|').append(g.a()).append(',').append(g.b());
			}
			if (!seenGridSets.add(sig.toString())) {
				continue; // e.g. ALL picked exactly the UW grids
			}
			for (BudType[] order : ALL_ORDERINGS) {
				out.add(decompose(a.r, ids, grids, order));
			}
			for (BudType[] order : ALL_ORDERINGS) {
				out.add(decomposeFewestOptionsFirst(a.r, ids, grids, order));
			}
		}
		return out;
	}

	/** Combined buds first (given), then the greedy single-type grouping of the rest. */
	private static BudDecomposition decompose(int r, int[][] ids, List<GridBud> grids, BudType[] order) {
		boolean[] used = new boolean[r];
		for (GridBud g : grids) for (int t : g.terms()) used[t] = true;
		List<Bud> buds = new ArrayList<>();
		for (BudType t : order) {
			groupBy(r, ids[t.ordinal()], t, used, buds);
		}
		List<Integer> trivial = new ArrayList<>();
		for (int l = 0; l < r; l++) if (!used[l]) trivial.add(l);
		return new BudDecomposition(buds, grids, trivial.stream().mapToInt(Integer::intValue).toArray());
	}

	/**
	 * Single-type grouping that serves the most constrained terms first (Karp–Sipser's
	 * matching heuristic, lifted from pairs to classes). The type-by-type greedy of
	 * {@link #decompose} commits a whole class at once, so a three-term U-class
	 * {@code {a,b,c}} swallows {@code c} even when {@code c} is the ONLY partner of some
	 * {@code d} on another factor — one bud where {@code {a,b}} + {@code {c,d}} gives
	 * two. No type ordering repairs that when the pattern occurs on two factors at
	 * once: ⟨4,5,9⟩:132 ⊗ˢ ⟨4,3,3⟩ stays at 3748 under all six orderings, 3744 here.
	 *
	 * <p>Loop, until no unassigned term has a class with a partner left:</p>
	 * <ol>
	 *   <li>a term with exactly ONE class that still holds another unassigned term
	 *       opens a bud there, with the class-mate that has the fewest other options
	 *       (pairs before joins: {@code d} claims {@code c} before {@code c} is
	 *       absorbed);</li>
	 *   <li>else every term with a class that is already a bud joins it;</li>
	 *   <li>else the class with the most unassigned terms opens with its two
	 *       least-flexible members.</li>
	 * </ol>
	 * {@code order} breaks ties between factors. Structural — no rank oracle.
	 */
	private static BudDecomposition decomposeFewestOptionsFirst(int r, int[][] ids, List<GridBud> grids,
			BudType[] order) {
		boolean[] done = new boolean[r];
		for (GridBud g : grids) for (int t : g.terms()) done[t] = true;
		int types = BudType.values().length;
		int[][] free = new int[types][];        // unassigned terms per class
		int[][] opened = new int[types][];      // terms already committed per class
		for (int ty = 0; ty < types; ty++) {
			int max = 0;
			for (int t = 0; t < r; t++) max = Math.max(max, ids[ty][t] + 1);
			free[ty] = new int[max];
			opened[ty] = new int[max];
			for (int t = 0; t < r; t++) if (!done[t]) free[ty][ids[ty][t]]++;
		}
		int[] assignedType = new int[r];
		java.util.Arrays.fill(assignedType, -1);
		while (true) {
			// 1. a term with exactly one usable class — the one with the fewest possible
			//    partners first (a term whose only class is a pair before a term whose only
			//    class is a triple, or the triple's pick could strand the pair's partner).
			int seed = -1, seedType = -1, seedChoices = Integer.MAX_VALUE;
			for (int t = 0; t < r; t++) {
				if (done[t]) continue;
				int usable = 0, last = -1;
				for (BudType ty : order) {
					if (free[ty.ordinal()][ids[ty.ordinal()][t]] >= 2) { usable++; last = ty.ordinal(); }
				}
				if (usable == 1 && free[last][ids[last][t]] < seedChoices) {
					seed = t;
					seedType = last;
					seedChoices = free[last][ids[last][t]];
				}
			}
			if (seed < 0) {
				// 2. join an open bud.
				boolean progressed = false;
				for (int t = 0; t < r; t++) {
					if (done[t]) continue;
					for (BudType ty : order) {
						if (opened[ty.ordinal()][ids[ty.ordinal()][t]] > 0) {
							commit(t, ty.ordinal(), ids, done, free, opened, assignedType);
							progressed = true;
							break;
						}
					}
				}
				if (progressed) continue;
			}
			// 3. else the fullest class, its least flexible member.
			if (seed < 0) {
				int best = 1;
				for (BudType ty : order) {
					for (int t = 0; t < r; t++) {
						if (done[t]) continue;
						int size = free[ty.ordinal()][ids[ty.ordinal()][t]];
						if (size > best) { best = size; seedType = ty.ordinal(); seed = t; }
					}
				}
				if (seed < 0) break; // nobody has a partner left
				seed = leastFlexible(seedType, ids[seedType][seed], -1, r, ids, done, free);
			}
			int mate = leastFlexible(seedType, ids[seedType][seed], seed, r, ids, done, free);
			commit(seed, seedType, ids, done, free, opened, assignedType);
			commit(mate, seedType, ids, done, free, opened, assignedType);
		}
		// Emit buds in the same shape as the greedy: per type in `order`, classes first-seen.
		List<Bud> buds = new ArrayList<>();
		for (BudType ty : order) {
			java.util.LinkedHashMap<Integer, List<Integer>> byClass = new java.util.LinkedHashMap<>();
			for (int t = 0; t < r; t++) {
				if (assignedType[t] == ty.ordinal()) {
					byClass.computeIfAbsent(ids[ty.ordinal()][t], k -> new ArrayList<>()).add(t);
				}
			}
			for (List<Integer> grp : byClass.values()) {
				buds.add(new Bud(ty, grp.stream().mapToInt(Integer::intValue).toArray()));
			}
		}
		List<Integer> trivial = new ArrayList<>();
		for (int t = 0; t < r; t++) if (assignedType[t] < 0 && !inGrid(grids, t)) trivial.add(t);
		return new BudDecomposition(buds, grids, trivial.stream().mapToInt(Integer::intValue).toArray());
	}

	private static boolean inGrid(List<GridBud> grids, int term) {
		for (GridBud g : grids) for (int t : g.terms()) if (t == term) return true;
		return false;
	}

	private static void commit(int t, int type, int[][] ids, boolean[] done, int[][] free, int[][] opened,
			int[] assignedType) {
		done[t] = true;
		assignedType[t] = type;
		for (int ty = 0; ty < free.length; ty++) free[ty][ids[ty][t]]--;
		opened[type][ids[type][t]]++;
	}

	/** The unassigned member of a class (other than {@code except}) with the fewest usable classes; lowest index on ties. */
	private static int leastFlexible(int type, int classId, int except, int r, int[][] ids, boolean[] done,
			int[][] free) {
		int best = -1, bestOptions = Integer.MAX_VALUE;
		for (int t = 0; t < r; t++) {
			if (done[t] || t == except || ids[type][t] != classId) continue;
			int options = 0;
			for (int ty = 0; ty < free.length; ty++) if (free[ty][ids[ty][t]] >= 2) options++;
			if (options < bestOptions) { bestOptions = options; best = t; }
		}
		return best;
	}

	/**
	 * Greedy, largest-first, disjoint selection of combined buds. Candidates: for every
	 * pair of row-classes sharing ≥ 2 column-classes, the grid over their common
	 * columns, and that grid extended by every further row-class containing all those
	 * columns. One term per cell (a second term in the same cell stays available to the
	 * single-type greedy). Deterministic: ties break on grid type, then on the terms.
	 */
	static List<GridBud> selectGrids(int r, int[][] ids, GridStrategy strategy) {
		if (strategy == GridStrategy.NONE) {
			return List.of();
		}
		List<GridBud> candidates = new ArrayList<>();
		boolean keepU = strategy == GridStrategy.CLOSED || strategy == GridStrategy.KEEP_U;
		boolean keepV = strategy == GridStrategy.CLOSED || strategy == GridStrategy.KEEP_V;
		boolean keepW = strategy == GridStrategy.CLOSED || strategy == GridStrategy.KEEP_W;
		boolean core = strategy == GridStrategy.CORE;
		boolean every = strategy == GridStrategy.ALL || keepU || keepV || keepW || core;
		// (rows, cols, free factor) per grid type: UW = (U, W, V), UV = (U, V, W), VW = (W, V, U).
		if (every || strategy == GridStrategy.UW) {
			gridCandidates(r, ids[0], ids[2], core ? ids[1] : null, GridStrategy.UW, keepU, keepW, candidates);
		}
		if (every || strategy == GridStrategy.UV) {
			gridCandidates(r, ids[0], ids[1], core ? ids[2] : null, GridStrategy.UV, keepU, keepV, candidates);
		}
		if (every || strategy == GridStrategy.VW) {
			gridCandidates(r, ids[2], ids[1], core ? ids[0] : null, GridStrategy.VW, keepW, keepV, candidates);
		}
		candidates.sort((x, y) -> {
			if (x.terms().length != y.terms().length) return Integer.compare(y.terms().length, x.terms().length);
			int tx = x.a() == 1 ? 0 : x.b() == 1 ? 1 : 2, ty = y.a() == 1 ? 0 : y.b() == 1 ? 1 : 2;
			if (tx != ty) return Integer.compare(tx, ty);
			return java.util.Arrays.compare(x.terms(), y.terms());
		});
		boolean[] used = new boolean[r];
		List<GridBud> picked = new ArrayList<>();
		next: for (GridBud g : candidates) {
			for (int t : g.terms()) if (used[t]) continue next;
			for (int t : g.terms()) used[t] = true;
			picked.add(g);
		}
		return picked;
	}

	/** Enumerate the grids of one type: {@code rows} / {@code cols} are per-term class ids. */
	private static void gridCandidates(int r, int[] rows, int[] cols, int[] coreFree, GridStrategy type,
			boolean rowsWhole, boolean colsWhole, List<GridBud> out) {
		// row-class → (col-class → first term in that cell), both in first-seen order.
		java.util.LinkedHashMap<Integer, java.util.LinkedHashMap<Integer, Integer>> cells = new java.util.LinkedHashMap<>();
		java.util.Map<Integer, Integer> rowSize = new java.util.HashMap<>(), colSize = new java.util.HashMap<>();
		for (int t = 0; t < r; t++) {
			cells.computeIfAbsent(rows[t], k -> new java.util.LinkedHashMap<>()).putIfAbsent(cols[t], t);
			rowSize.merge(rows[t], 1, Integer::sum);
			colSize.merge(cols[t], 1, Integer::sum);
		}
		java.util.function.BiPredicate<List<Integer>, List<Integer>> admissible = (gridRows, gridCols) -> {
			if (rowsWhole) for (Integer row : gridRows) if (rowSize.get(row) != gridCols.size()) return false;
			if (colsWhole) for (Integer col : gridCols) if (colSize.get(col) != gridRows.size()) return false;
			return true;
		};
		List<Integer> rowIds = new ArrayList<>();
		for (var e : cells.entrySet()) if (e.getValue().size() >= 2) rowIds.add(e.getKey());
		java.util.Set<String> seen = new java.util.HashSet<>();
		for (int x = 0; x < rowIds.size(); x++) {
			for (int y = x + 1; y < rowIds.size(); y++) {
				var cx = cells.get(rowIds.get(x));
				var cy = cells.get(rowIds.get(y));
				List<Integer> common = new ArrayList<>();
				for (Integer c : cx.keySet()) if (cy.containsKey(c)) common.add(c);
				if (common.size() < 2) continue;
				List<Integer> pair = List.of(rowIds.get(x), rowIds.get(y));
				List<Integer> extended = new ArrayList<>(pair);
				for (int z = 0; z < rowIds.size(); z++) {
					if (z == x || z == y) continue;
					if (cells.get(rowIds.get(z)).keySet().containsAll(common)) extended.add(rowIds.get(z));
				}
				if (coreFree != null && (extended.size() > 2 || common.size() > 2)) {
					// CORE: a biclique larger than 2×2 is replaced by its least-entangled 2×2.
					addCore(r, cells, extended, common, coreFree, type, seen, out);
					continue;
				}
				if (admissible.test(pair, common)) addGrid(cells, pair, common, type, seen, out);
				if (extended.size() > 2 && admissible.test(extended, common)) {
					addGrid(cells, extended, common, type, seen, out);
				}
			}
		}
	}

	/** The 2×2 sub-grid of {@code gridRows × gridCols} with the fewest terms that also
	 *  belong to a (≥ 2-term) class of the free factor; first found on ties. */
	private static void addCore(int r, java.util.Map<Integer, java.util.LinkedHashMap<Integer, Integer>> cells,
			List<Integer> gridRows, List<Integer> gridCols, int[] freeIds, GridStrategy type,
			java.util.Set<String> seen, List<GridBud> out) {
		java.util.Map<Integer, Integer> freeSize = new java.util.HashMap<>();
		for (int t = 0; t < r; t++) freeSize.merge(freeIds[t], 1, Integer::sum);
		List<Integer> rs = new ArrayList<>(gridRows);
		java.util.Collections.sort(rs);
		int best = Integer.MAX_VALUE;
		List<Integer> bestRows = null, bestCols = null;
		for (int x = 0; x < rs.size(); x++) for (int y = x + 1; y < rs.size(); y++) {
			for (int u = 0; u < gridCols.size(); u++) for (int v = u + 1; v < gridCols.size(); v++) {
				int entangled = 0;
				for (int row : new int[] { rs.get(x), rs.get(y) }) for (int col : new int[] { gridCols.get(u), gridCols.get(v) }) {
					if (freeSize.get(freeIds[cells.get(row).get(col)]) >= 2) entangled++;
				}
				if (entangled < best) {
					best = entangled;
					bestRows = List.of(rs.get(x), rs.get(y));
					bestCols = List.of(gridCols.get(u), gridCols.get(v));
				}
			}
		}
		addGrid(cells, bestRows, bestCols, type, seen, out);
	}

	private static void addGrid(java.util.Map<Integer, java.util.LinkedHashMap<Integer, Integer>> cells,
			List<Integer> gridRows, List<Integer> gridCols, GridStrategy type, java.util.Set<String> seen,
			List<GridBud> out) {
		List<Integer> rs = new ArrayList<>(gridRows);
		java.util.Collections.sort(rs);
		if (!seen.add(rs + "×" + gridCols)) return;
		int nr = rs.size(), nc = gridCols.size();
		int[] terms = new int[nr * nc];
		// Row-major (i,j,k) with the unit axis dropped — see GridBud:
		//   UW ⟨1,b,c⟩: rows = U-classes (j), cols = W-classes (k) → j·c + k
		//   UV ⟨a,1,c⟩: rows = U-classes (i), cols = V-classes (k) → i·c + k
		//   VW ⟨a,b,1⟩: rows = W-classes (i), cols = V-classes (j) → i·b + j
		for (int x = 0; x < nr; x++) for (int y = 0; y < nc; y++) {
			terms[x * nc + y] = cells.get(rs.get(x)).get(gridCols.get(y));
		}
		out.add(switch (type) {
			case UW -> new GridBud(1, nr, nc, terms);
			case UV -> new GridBud(nr, 1, nc, terms);
			case VW -> new GridBud(nr, nc, 1, terms);
			default -> throw new IllegalArgumentException("not a grid type: " + type);
		});
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
		boolean[] used = new boolean[a.r];
		List<Bud> buds = new ArrayList<>();
		groupBy(a.r, classIds(a.denseU(), a.r), BudType.U, used, buds);
		List<Integer> trivial = new ArrayList<>();
		for (int l = 0; l < a.r; l++) if (!used[l]) trivial.add(l);
		return new BudDecomposition(buds, trivial.stream().mapToInt(Integer::intValue).toArray());
	}

	/** Group the still-unused terms by their class id on one factor ({@code ids} from
	 *  {@link #independentClassIds}: equal id ⇔ proportional columns). */
	private static void groupBy(int r, int[] ids, BudType type, boolean[] used, List<Bud> buds) {
		java.util.LinkedHashMap<Integer, List<Integer>> byDir = new java.util.LinkedHashMap<>();
		for (int l = 0; l < r; l++) {
			if (used[l]) continue;
			byDir.computeIfAbsent(ids[l], k -> new ArrayList<>()).add(l);
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
		for (BudDecomposition dec : candidateDecompositions(t1)) {
			best = Math.min(best, costOf(dec, lookup, n2, m2, p2));
		}
		return best;
	}

	/** Predicted cost of a specific decomposition; {@code Long.MAX_VALUE/4} if any
	 *  enlarged inner rank is unknown. */
	public static long costOf(BudDecomposition dec, FieldAwareLookup lookup, int n2, int m2, int p2) {
		return costOf(dec, RankOracle.of(lookup), n2, m2, p2);
	}

	/** {@link #costOf(BudDecomposition, FieldAwareLookup, int, int, int)} against any oracle. */
	public static long costOf(BudDecomposition dec, RankOracle oracle, int n2, int m2, int p2) {
		final long UNKNOWN = Long.MAX_VALUE / 4;
		long inner = oracle.rank(n2, m2, p2);
		if (inner >= Recombination.SotaResolver.UNKNOWN_RANK) return UNKNOWN;
		long cost = (long) dec.trivial().length * inner;
		for (Bud b : dec.buds()) {
			int k = b.terms().length;
			long r = switch (b.type()) {
				case U -> oracle.rank(n2, m2, k * p2);
				case V -> oracle.rank(k * n2, m2, p2);
				case W -> oracle.rank(n2, k * m2, p2);
			};
			if (r >= Recombination.SotaResolver.UNKNOWN_RANK) return UNKNOWN;
			cost += r;
		}
		for (GridBud g : dec.grids()) {
			long r = oracle.rank(g.a() * n2, g.b() * m2, g.c() * p2);
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
		// A built product's rank IS its priced cost (the parts' ranks add up), so price
		// every candidate with what the resolver can build and build cheapest-first —
		// one build instead of one per candidate. Stable sort: ties keep family order.
		RankOracle oracle = RankOracle.of(resolver);
		List<BudDecomposition> cands = new ArrayList<>();
		java.util.Map<BudDecomposition, Long> price = new java.util.IdentityHashMap<>();
		for (BudDecomposition dec : candidateDecompositions(t1)) {
			if (!dec.hasBuds()) {
				continue; // nothing fused → nothing to beat the plain product with
			}
			long cost = costOf(dec, oracle, n2, m2, p2);
			if (cost >= Long.MAX_VALUE / 4) {
				continue; // a fusion target is unavailable
			}
			cands.add(dec);
			price.put(dec, cost);
		}
		cands.sort(java.util.Comparator.comparingLong(price::get));
		for (BudDecomposition dec : cands) {
			try {
				return productFromDecomposition(t1, dec, resolver, n2, m2, p2, allow);
			} catch (RuntimeException e) {
				// Priced but not buildable after all — fall through to the next cheapest.
			}
		}
		// No candidate yielded a buildable bud decomposition → default (term-by-term Kron).
		return productFromDecomposition(
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
		for (GridBud grid : dec.grids()) {
			if (allow.isEmpty()) {
				for (int term : grid.terms()) {
					parts.add(Compose.kroneckerGeneral(rankOne(t1, term), s2));
				}
				continue;
			}
			int en = grid.a() * n2, em = grid.b() * m2, ep = grid.c() * p2;
			NonCubicBilinearAlgorithm s3 = resolver.find(en, em, ep).orElseThrow(
					() -> new IllegalStateException("no buildable enlarged ⟨" + en + "," + em + ","
							+ ep + "⟩ for a combined bud (stub-only? pass a replaying InnerResolver)"));
			parts.add(buildGridBlock(t1, grid, s3, n2, m2, p2));
		}
		return concatColumns(parts, t1.n * n2, t1.m * m2, t1.p * p2);
	}

	/**
	 * Fused block of a combined bud {@code ⟨a,b,c⟩} against the inner {@code ⟨n₂,m₂,p₂⟩},
	 * realised by one {@code S3 = ⟨a·n₂, b·m₂, c·p₂⟩}. With class representatives
	 * {@code Û[i][j]} (base A-coefficients, {@code n₁×m₁}), {@code V̂[j][k]}, {@code Ŵ[i][k]}
	 * such that the bud's terms sum to {@code Σ_{ijk} Û_ij ⊗ V̂_jk ⊗ Ŵ_ik}, product
	 * {@code q} of the block is
	 * <pre>
	 *   U_q = Σ_{i,j} Û_ij ⊗ S3.U_q[block (i,j)]     (block = the n₂×m₂ sub-matrix)
	 *   V_q = Σ_{j,k} V̂_jk ⊗ S3.V_q[block (j,k)]
	 *   W_q = Σ_{i,k} Ŵ_ik ⊗ S3.W_q[block (i,k)]
	 * </pre>
	 * Summing over {@code q}, S3's matmul identity kills every mismatched block triple
	 * and leaves {@code ⟨n₂,m₂,p₂⟩} on the matched ones — i.e. exactly
	 * {@code (Σ Û_ij⊗V̂_jk⊗Ŵ_ik) ⊗ ⟨n₂,m₂,p₂⟩}. The single-type blocks of
	 * {@link #buildBudBlock} are this formula with two unit dims.
	 */
	private static NonCubicBilinearAlgorithm buildGridBlock(
			NonCubicBilinearAlgorithm t1, GridBud grid, NonCubicBilinearAlgorithm s3,
			int n2, int m2, int p2) {
		double[][] srcU = t1.denseU(), srcV = t1.denseV(), srcW = t1.denseW();
		int a = grid.a(), b = grid.b(), c = grid.c();
		int n1 = t1.n, m1 = t1.m, p1 = t1.p;
		double[][][] uHat = new double[a][b][], vHat = new double[b][c][], wHat = new double[a][c][];
		for (int i = 0; i < a; i++) for (int j = 0; j < b; j++) for (int k = 0; k < c; k++) {
			int term = grid.terms()[(i * b + j) * c + k];
			double[] u = column(srcU, term), v = column(srcV, term), w = column(srcW, term);
			if (a == 1) {
				// U-class j, W-class k; V is free and absorbs both scales.
				if (uHat[0][j] == null) uHat[0][j] = u;
				if (wHat[0][k] == null) wHat[0][k] = w;
				vHat[j][k] = mul(v, proportionFactor(uHat[0][j], u) * proportionFactor(wHat[0][k], w));
			} else if (b == 1) {
				// U-class i, V-class k; W is free.
				if (uHat[i][0] == null) uHat[i][0] = u;
				if (vHat[0][k] == null) vHat[0][k] = v;
				wHat[i][k] = mul(w, proportionFactor(uHat[i][0], u) * proportionFactor(vHat[0][k], v));
			} else {
				// W-class i, V-class j; U is free.
				if (wHat[i][0] == null) wHat[i][0] = w;
				if (vHat[j][0] == null) vHat[j][0] = v;
				uHat[i][j] = mul(u, proportionFactor(wHat[i][0], w) * proportionFactor(vHat[j][0], v));
			}
		}
		int N = n1 * n2, M = m1 * m2, P = p1 * p2, r3 = s3.r;
		int m3 = b * m2, p3 = c * p2;
		double[][] s3U = s3.denseU(), s3V = s3.denseV(), s3W = s3.denseW();
		double[][] U = new double[N * M][r3], V = new double[M * P][r3], W = new double[N * P][r3];
		for (int i = 0; i < a; i++) for (int j = 0; j < b; j++) {
			double[] hat = uHat[i][j];
			for (int i1 = 0; i1 < n1; i1++) for (int j1 = 0; j1 < m1; j1++) {
				double coef = hat[i1 * m1 + j1]; if (coef == 0) continue;
				for (int i2 = 0; i2 < n2; i2++) for (int j2 = 0; j2 < m2; j2++) {
					double[] src = s3U[(i * n2 + i2) * m3 + (j * m2 + j2)];
					double[] dst = U[(i1 * n2 + i2) * M + (j1 * m2 + j2)];
					for (int q = 0; q < r3; q++) if (src[q] != 0) dst[q] += coef * src[q];
				}
			}
		}
		for (int j = 0; j < b; j++) for (int k = 0; k < c; k++) {
			double[] hat = vHat[j][k];
			for (int j1 = 0; j1 < m1; j1++) for (int k1 = 0; k1 < p1; k1++) {
				double coef = hat[j1 * p1 + k1]; if (coef == 0) continue;
				for (int j2 = 0; j2 < m2; j2++) for (int k2 = 0; k2 < p2; k2++) {
					double[] src = s3V[(j * m2 + j2) * p3 + (k * p2 + k2)];
					double[] dst = V[(j1 * m2 + j2) * P + (k1 * p2 + k2)];
					for (int q = 0; q < r3; q++) if (src[q] != 0) dst[q] += coef * src[q];
				}
			}
		}
		for (int i = 0; i < a; i++) for (int k = 0; k < c; k++) {
			double[] hat = wHat[i][k];
			for (int i1 = 0; i1 < n1; i1++) for (int k1 = 0; k1 < p1; k1++) {
				double coef = hat[i1 * p1 + k1]; if (coef == 0) continue;
				for (int i2 = 0; i2 < n2; i2++) for (int k2 = 0; k2 < p2; k2++) {
					double[] src = s3W[(i * n2 + i2) * p3 + (k * p2 + k2)];
					double[] dst = W[(i1 * n2 + i2) * P + (k1 * p2 + k2)];
					for (int q = 0; q < r3; q++) if (src[q] != 0) dst[q] += coef * src[q];
				}
			}
		}
		return new NonCubicBilinearAlgorithm(N, M, P, U, V, W);
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
