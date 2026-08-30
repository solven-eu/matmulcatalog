package eu.solven.matmul.catalog;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Optional;

import eu.solven.matmul.algebra.Algebra;
import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.KnownAlgorithm.Optimality;

/**
 * Curated list of known fast bilinear matrix-multiplication algorithms,
 * ordered chronologically per {@code (format, algebra)}. Mirrors and is
 * cross-checked with {@code SMALL_MATMUL_CATALOG.md}.
 *
 * <p>**This is a starting set.** Pull requests welcome to extend
 * coverage, verify ranks against {@code [drons]}'s structured database,
 * and add the AlphaTensor 2022 non-cubic improvements not yet captured
 * here.</p>
 *
 * <p><strong>Constructions only.</strong> Every row must be an algorithm
 * that attains its {@code rank}. A rank-only claim with no factor matrices
 * belongs in {@code docs/cited-bounds.json} (CLAUDE.md, "Scheme registration"
 * rule 2); an impossibility result belongs in {@code docs/lower-bounds.json}.
 * Each row carries an {@link KnownAlgorithm.Optimality} tier, and
 * {@code TestBoundsVsLowerBounds} rejects any row that falls below a
 * published floor — the guard that would have caught the commutative
 * {@code ⟨2,2,2⟩ = 6} row (attributed, for years, to the two 1971 papers
 * that prove 7 is minimal).</p>
 */
public final class KnownAlgorithmCatalog {

	private static final List<KnownAlgorithm> ALL = new ArrayList<>();

	private static final Algebra NC_R = Algebra.nonCommutative(Field.R);
	private static final Algebra NC_F2 = Algebra.nonCommutative(Field.F2);
	private static final Algebra NC_C = Algebra.nonCommutative(Field.C);
	private static final Algebra CMT_R = Algebra.commutative(Field.R);

	// Optimality tiers (CLAUDE.md optimality discipline). PROVEN means a matching
	// floor is published in docs/lower-bounds.json; TestBoundsVsLowerBounds asserts
	// both directions, so neither an over-claim nor a stale BOUND can survive.
	private static final Optimality PROVEN = Optimality.PROVEN_OPTIMAL;
	private static final Optimality BOUND = Optimality.BOUND;

	static {
		// ───────── ⟨2,2,2⟩ ─────────
		// Strassen's algorithm uses ±1 coefficients — field-agnostic; explicitly
		// catalogue under every algebra we track so RecursiveComposition can compose.
		add(2, 2, 2, NC_R, 7, PROVEN, 1969, "Strassen",
				"https://doi.org/10.1007/BF02165411",
				"Original recursive matmul; non-commutative, valid for matrix entries.");
		add(2, 2, 2, NC_F2, 7, PROVEN, 1969, "Strassen",
				"https://doi.org/10.1007/BF02165411",
				"Strassen's ±1 coefficients reduce mod 2 to a valid F₂ scheme.");
		add(2, 2, 2, NC_C, 7, PROVEN, 1969, "Strassen",
				"https://doi.org/10.1007/BF02165411",
				"Strassen's algorithm works unchanged over C.");
		add(2, 2, 2, CMT_R, 7, PROVEN, 1969, "Strassen",
				"https://doi.org/10.1007/BF02165411",
				"Also optimal over commutative scalars; Winograd proved the matching lower bound in 1971.");

		// ───────── ⟨3,3,3⟩ ─────────
		add(3, 3, 3, NC_R, 23, BOUND, 1976, "Laderman",
				"https://doi.org/10.1090/S0002-9904-1976-13988-2",
				"Best known since 1976; no improvement in 49 years. LB is 19 (Bläser 2003).");
		add(3, 3, 3, NC_F2, 23, BOUND, 1976, "Laderman",
				"https://doi.org/10.1090/S0002-9904-1976-13988-2",
				"Laderman's scheme reduces mod 2; AlphaTensor 2022 did NOT improve ⟨3,3,3⟩ over F₂. "
				+ "LB over F2 tightened to 20 by Wang (arXiv:2603.07280, 2026).");

		// ───────── ⟨2,2,3⟩ ─────────
		add(2, 2, 3, NC_R, 11, PROVEN, 1971, "Hopcroft–Kerr",
				"https://doi.org/10.1137/0120004", "Tight over fields with char 0.");

		// ───────── ⟨2,3,3⟩ ─────────
		add(2, 3, 3, NC_R, 15, PROVEN, 1973, "Hopcroft–Kerr / Pan",
				null, "Tight at 15 over R, Q.");

		// ───────── ⟨4,4,4⟩ ─────────
		add(4, 4, 4, NC_R, 49, BOUND, 1969, "Strassen² (recursive)",
				"https://doi.org/10.1007/BF02165411",
				"7 × 7 from two-level Strassen; still best known over R as of 2025.");
		add(4, 4, 4, NC_F2, 47, BOUND, 2022, "AlphaTensor",
				"https://www.nature.com/articles/s41586-022-05172-4",
				"RL-discovered; first improvement over Strassen² in this setting.");
		add(4, 4, 4, NC_C, 48, BOUND, 2025, "AlphaEvolve",
				"https://arxiv.org/abs/2506.13131",
				"First improvement on Strassen² over C 'after 56 years' (Novikov et al. 2025).");

		// ───────── ⟨5,5,5⟩ ─────────
		add(5, 5, 5, NC_F2, 96, BOUND, 2022, "AlphaTensor",
				"https://www.nature.com/articles/s41586-022-05172-4",
				"Improves Smirnov's prior best by RL search.");

		// ───────── ⟨4,5,5⟩ (AlphaTensor non-cubic improvement) ─────────
		add(4, 5, 5, NC_F2, 76, BOUND, 2022, "AlphaTensor",
				"https://www.nature.com/articles/s41586-022-05172-4",
				"Down from 80; one of the AlphaTensor 2022 highlight improvements.");
	}

	private static void add(int n, int m, int p, Algebra a, int rank, Optimality optimality,
			int year, String source, String link, String notes) {
		ALL.add(new KnownAlgorithm(n, m, p, a, rank, optimality, year, source, link, notes));
	}

	private KnownAlgorithmCatalog() {}

	public static List<KnownAlgorithm> all() {
		return new ArrayList<>(ALL);
	}

	/** All entries for a given {@code (format, algebra)} ordered by ascending year (oldest first). */
	public static List<KnownAlgorithm> historyFor(int n, int m, int p, Algebra algebra) {
		List<KnownAlgorithm> hits = new ArrayList<>();
		for (KnownAlgorithm a : ALL) {
			if (a.n == n && a.m == m && a.p == p && a.algebra.equals(algebra)) hits.add(a);
		}
		hits.sort(Comparator.comparingInt((KnownAlgorithm a) -> a.year));
		return hits;
	}

	/** Best (minimum) known rank for the given {@code (format, algebra)}, if any. */
	public static Optional<KnownAlgorithm> bestKnown(int n, int m, int p, Algebra algebra) {
		return historyFor(n, m, p, algebra).stream()
				.min(Comparator.comparingInt(a -> a.rank));
	}
}
