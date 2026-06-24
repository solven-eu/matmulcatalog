package eu.solven.matmul.search;

import eu.solven.matmul.recombination.Recombination;

import java.util.List;

import eu.solven.matmul.catalog.FieldAwareLookup;

/**
 * Hand-extracted catalog of Schönhage τ-theorem / disjoint-sum identities
 * published by FMM-Lille, Smirnov, Sedoglavic, et al. Each identity
 * expresses a target shape {@code ⟨n,m,p⟩} as a sum (with optional
 * "borrow-and-correct" coefficients) of smaller catalog atoms, with the
 * predicted total rank that arises from naively summing the sub-atom
 * ranks possibly with a τ-discount.
 *
 * <h2>Two identity families</h2>
 *
 * <p><strong>1. Plain disjoint sum.</strong> {@code rank(target) ≤ Σ rank(Tᵢ)}
 * over a multiset of sub-tensors {@code Tᵢ} whose joint computation embeds
 * the target tensor. The classical Schönhage τ-theorem covers when the
 * embedding exists; the FMM-Lille recipes are concrete instances.</p>
 *
 * <p>Example: <strong>⟨11,11,11⟩=873</strong>:
 * {@code R(⟨11,11,11⟩) ≤ R(⟨6,6,6⟩=153) + 3·R(⟨5,5,6⟩=110) + 3·R(⟨5,6,6⟩=130) = 873}
 * (FMM-Lille; underlying construction Smirnov / Sedoglavic family).</p>
 *
 * <p><strong>2. Serendipitous tensor product (Pan-style borrow).</strong>
 * {@code rank(⟨n₁n₂, m₁m₂, p₁p₂⟩) ≤ (rA − k)·rB + c·rCaux} — borrow {@code k}
 * sub-products from {@code A}, replace them with {@code c} copies of an
 * auxiliary shape {@code Caux}. Strictly more general than naïve Kronecker.</p>
 *
 * <p>Example: <strong>⟨6,8,9⟩=296</strong>:
 * {@code (R(⟨2,3,4⟩=20) − 8) · R(⟨3,3,2⟩=15) + 4·R(⟨3,3,4⟩=29) = 12·15 + 116 = 296}
 * (FMM-Lille; vs naïve Kronecker 20·15=300).</p>
 *
 * <h2>How to read entries</h2>
 *
 * <p>Each {@link Identity} records:</p>
 * <ul>
 *   <li>The target shape and predicted rank.</li>
 *   <li>The sub-shapes referenced (each lookup-able in the catalog via
 *       {@link FieldAwareLookup#findRank}).</li>
 *   <li>The arithmetic operator (plain sum vs borrow-and-correct).</li>
 *   <li>The source attribution.</li>
 * </ul>
 *
 * <p>An automatically-runnable {@link #verify} method checks each identity
 * against the live catalog — for plain sums, it asserts that summing the
 * catalog ranks of the referenced sub-atoms reproduces the predicted total.
 * For borrow-and-correct, it applies the {@code (rA−k)·rB + c·rC} formula.
 * Identities that fail verification (because the catalog has changed) get
 * flagged at next regen.</p>
 *
 * <p>This class is the staging ground for {@code DisjointSumSearch} +
 * {@code SerendipitousTensorSearch} (#159 / #102): once a critical mass of
 * identities is enumerated here, the search wraps them as candidates in
 * {@code BlockSplitSearch.findBestStrategy}.</p>
 */
public final class KnownTauIdentities {

	private KnownTauIdentities() {}

	/** A single τ-identity. Use static factory methods to construct. */
	public sealed interface Identity {
		int[] target();
		long predictedRank();
		String attribution();
		String description();
	}

	/**
	 * Plain disjoint-sum identity: {@code R(target) ≤ Σᵢ multiplicityᵢ · R(subShapeᵢ)}.
	 */
	public record DisjointSum(
			int[] target,
			List<SubTerm> terms,
			long predictedRank,
			String attribution,
			String description) implements Identity {}

	/**
	 * A term in a disjoint sum: {@code multiplicity · cost(n,m,p)}, where
	 * {@code cost} is either a catalog SOTA lookup ({@link Kind#ATOM_RANK})
	 * or Pan 1980's pair-fusion cost {@code abc + ab + bc + ca}
	 * ({@link Kind#PAN_PAIR_COST}).
	 *
	 * <p>{@link Kind#PAN_PAIR_COST} is required when an identity sums
	 * SAME-shape cyclically-related cubic products that get pair-fused —
	 * the cost is the Pan TA pair_cost, NOT a separate catalog atom.
	 * Example: Sedoglavic doubling ⟨14⟩³=1719 = ⟨7⟩³=249 + 3·pair_cost(7,7,7).</p>
	 */
	public record SubTerm(int n, int m, int p, int multiplicity, Kind kind) {
		public SubTerm(int n, int m, int p, int multiplicity) {
			this(n, m, p, multiplicity, Kind.ATOM_RANK);
		}
		public int[] shape() { return new int[] { n, m, p }; }
		public enum Kind { ATOM_RANK, PAN_PAIR_COST }
	}

	/**
	 * Serendipitous tensor product (Schönhage borrow-and-correct):
	 * {@code R(⟨nA·nB, mA·mB, pA·pB⟩) ≤ (rA − k) · rB + c · R(Caux)}.
	 */
	public record BorrowAndCorrect(
			int[] target,
			int[] shapeA, int[] shapeB, int[] shapeCaux,
			int discountK, int auxMultiplicity,
			long predictedRank,
			String attribution,
			String description) implements Identity {}

	// ───── Hand-extracted identities ─────

	/**
	 * FMM-Lille's published recipe for {@code ⟨11,11,11⟩=873}:
	 * <pre>
	 *   ⟨11,11,11⟩ = ⟨6,6,6⟩ + 3·⟨5,5,6⟩ + 3·⟨5,6,6⟩
	 *              = 153 + 3·110 + 3·130
	 *              = 873
	 * </pre>
	 * Equivalent under axis-permutation orbits:
	 * <pre>
	 *   ⟨6,6,6⟩=153, ⟨5,6,5⟩=⟨5,5,6⟩=110, ⟨6,5,6⟩=⟨5,6,6⟩=130, etc.
	 * </pre>
	 *
	 * <p>Note the dimension-sum-per-axis is 39, not 11 — this is NOT a
	 * block decomposition, but a Schönhage τ-aggregation where the 7
	 * sub-tensors share input/output variables in an overlap pattern.
	 * Compared to the 8-block Cartesian product of {@code [5,6]³} this
	 * recipe is "missing" the {@code ⟨5,5,5⟩} block — that block is
	 * absorbed via the overlap.</p>
	 */
	public static final DisjointSum FMM_11_11_11_873 = new DisjointSum(
			new int[] { 11, 11, 11 },
			List.of(
					new SubTerm(6, 6, 6, 1),
					new SubTerm(5, 5, 6, 3),
					new SubTerm(5, 6, 6, 3)),
			873L,
			"FMM-Lille",
			"⟨11,11,11⟩=873 = ⟨6,6,6⟩=153 + 3·⟨5,5,6⟩=110 + 3·⟨5,6,6⟩=130. "
					+ "Schönhage τ-aggregation; the 7 sub-tensors share variables (Σ axis-dim = 39 ≠ 11). "
					+ "The 'missing' ⟨5,5,5⟩ Cartesian block is absorbed via overlap.");

	/**
	 * <strong>Sedoglavic 2017 Proposition 1</strong> at (u,v)=(6,5) → {@code ⟨11,11,11⟩=873}:
	 * <pre>
	 *   ⟨u+v, u+v, u+v⟩ ≤ ⟨u,u,u⟩ + 3·⟨u,u,v⟩ + 3·⟨v,v,u⟩    when u > v
	 * </pre>
	 * This is the same numerical recipe as {@link #FMM_11_11_11_873} above —
	 * registered here as a SEPARATE identity to preserve the constructive
	 * attribution (Sedoglavic 2017 hal-01572046v2 is the methodology paper
	 * behind the FMM-Lille catalog for cubic decompositions). The two
	 * entries should always verify to the same predicted rank; if they
	 * diverge we have a catalog inconsistency.
	 *
	 * <p>Source paper: A. Sedoglavic, "A non-commutative algorithm for
	 * multiplying (7×7) matrices using 250 multiplications" (hal-01572046v2,
	 * Dec 2017). Proposition 1 with (u,v)=(4,3) gives the ⟨7,7,7⟩=250 result
	 * of the title; the same proposition for (u,v)=(6,5) lands on
	 * ⟨11,11,11⟩=873 — Lille's published bound.</p>
	 */
	public static final DisjointSum SEDOGLAVIC_PROP1_11_11_11 = new DisjointSum(
			new int[] { 11, 11, 11 },
			List.of(
					new SubTerm(6, 6, 6, 1),
					new SubTerm(6, 6, 5, 3),
					new SubTerm(5, 5, 6, 3)),
			873L,
			"Sedoglavic 2017 hal-01572046v2",
			"⟨u+v,u+v,u+v⟩ ≤ ⟨u,u,u⟩ + 3⟨u,u,v⟩ + 3⟨v,v,u⟩ at (u,v)=(6,5). "
					+ "This is the methodology paper behind the FMM-Lille cubic decomposition catalog.");

	/**
	 * Sedoglavic 2017 Proposition 1 at (u,v)=(4,3) → {@code ⟨7,7,7⟩=248}:
	 * <pre>
	 *   ⟨7,7,7⟩ ≤ ⟨4,4,4⟩=47 + 3·⟨4,4,3⟩=38 + 3·⟨3,3,4⟩=29 = 248
	 * </pre>
	 *
	 * <p><strong>NEW BOUND.</strong> Sedoglavic's paper title cites
	 * ⟨7,7,7⟩=250 because in 2017 the best ⟨4,4,4⟩ was Strassen²=49. With
	 * AlphaTensor 2022's ⟨4,4,4⟩=47 over F₂ plugged in, the same proposition
	 * yields 248 — a 1-multiplication improvement over our current catalog
	 * entry ⟨7,7,7⟩=249. Worth materialising and importing as the new SOTA.</p>
	 *
	 * <p>Note: validity hinges on ⟨4,4,4⟩=47 being usable in the recursion.
	 * AT-F2 47 is over F₂; the Sedoglavic recursion is non-commutative, so
	 * F₂ atoms compose only over F₂ targets. For R/Q/Z the right base is
	 * ⟨4,4,4⟩=48 (AlphaEvolve 2025 over C) — giving 248 + 1·3 = wait, let
	 * me recompute: 48 + 3·38 + 3·29 = 48 + 114 + 87 = 249 (matches current).
	 * So 248 is only achievable over F₂. Tagged commutative=false, field=F2
	 * — distinct from the standard R/Q/Z catalog entry.</p>
	 */
	public static final DisjointSum SEDOGLAVIC_PROP1_7_7_7 = new DisjointSum(
			new int[] { 7, 7, 7 },
			List.of(
					new SubTerm(4, 4, 4, 1),
					new SubTerm(4, 4, 3, 3),
					new SubTerm(3, 3, 4, 3)),
			// Over Q (AE ⟨4,4,4⟩=48): 48 + 3·38 + 3·29 = 249. We register the
			// Q value so verify() against the Q catalog passes. The F₂
			// improvement to 248 (AT-F2 ⟨4,4,4⟩=47) is a separate analysis
			// that needs a field-aware lookup. See #161 for the proper
			// per-field method evaluation.
			249L,
			"Sedoglavic 2017 hal-01572046v2",
			"⟨7,7,7⟩=249 over Q via Sedoglavic Prop 1 at (u,v)=(4,3). "
					+ "Originally 250 with Strassen² ⟨4,4,4⟩=49 (paper title bound); "
					+ "now 249 with AlphaEvolve ⟨4,4,4⟩=48. Over F₂ with AT-F2 ⟨4,4,4⟩=47 "
					+ "the same proposition gives 248 — a new bound queued for #161.");

	/**
	 * Sedoglavic 2017 Prop 1 <strong>doubling extension</strong> for u=v=k:
	 * <pre>
	 *   ⟨2k, 2k, 2k⟩ ≤ ⟨k,k,k⟩ + 3·pair_cost(k,k,k)
	 * </pre>
	 * where {@code pair_cost(k,k,k) = k³ + 3k²} is Pan 1980's pair-fusion
	 * cost of two cyclically-related cubic products. This extends Sedoglavic
	 * Prop 1 (which requires u > v) to u = v = k by using Pan TA pair-fusion
	 * to absorb the 6 same-shape cubic copies into 3 paired computations.
	 *
	 * <p>For k=7: 249 + 3·(343+147) = 249 + 1470 = <strong>1719</strong>,
	 * exactly the ⟨14,14,14⟩=1719 value FMM-Lille publishes and we hold in
	 * the catalog. Per user 2026-06-03: this rank is in our catalog but
	 * without lineage_compact attribution.</p>
	 */
	public static final DisjointSum SEDOGLAVIC_DOUBLING_14_14_14 = new DisjointSum(
			new int[] { 14, 14, 14 },
			List.of(
					new SubTerm(7, 7, 7, 1, SubTerm.Kind.ATOM_RANK),
					new SubTerm(7, 7, 7, 3, SubTerm.Kind.PAN_PAIR_COST)),
			1719L,
			"Sedoglavic 2017 + Pan 1980 TA-pair extension",
			"⟨14,14,14⟩=1719 via ⟨7,7,7⟩=249 + 3·pair_cost(7,7,7)=490. "
					+ "Closes user-observed gap: rank was in catalog (solven-strassen-2026) "
					+ "but lineage_compact attribution was missing.");

	/**
	 * FMM-Lille's published recipe for {@code ⟨6,8,9⟩=296}:
	 * <pre>
	 *   ⟨6,8,9⟩ = (⟨2,3,4⟩=20 − 8) ⊗ ⟨3,3,2⟩=15 + 4·⟨3,3,4⟩=29
	 *           = 12 · 15 + 4 · 29
	 *           = 180 + 116
	 *           = 296
	 * </pre>
	 * vs naïve Kronecker {@code ⟨2,3,4⟩=20 ⊗ ⟨3,3,2⟩=15 → ⟨6,9,8⟩=300}.
	 * The (k=8, c=4) discount saves 4 multiplications.
	 */
	public static final BorrowAndCorrect FMM_6_8_9_296 = new BorrowAndCorrect(
			new int[] { 6, 8, 9 },
			new int[] { 2, 3, 4 },
			new int[] { 3, 3, 2 },
			new int[] { 3, 3, 4 },
			/*discountK*/ 8,
			/*auxMultiplicity*/ 4,
			296L,
			"FMM-Lille",
			"⟨6,8,9⟩=296 = (20−8)·15 + 4·29. Saves 4 mults vs naïve Kron ⟨2,3,4⟩⊗⟨3,3,2⟩=300.");

	/** All identities currently registered. Extend as recipes are extracted. */
	public static final List<Identity> ALL = List.of(
			FMM_11_11_11_873,
			SEDOGLAVIC_PROP1_11_11_11,
			SEDOGLAVIC_PROP1_7_7_7,
			SEDOGLAVIC_DOUBLING_14_14_14,
			FMM_6_8_9_296);

	/**
	 * Identities whose target matches {@code ⟨n,m,p⟩} (under axis-permutation
	 * orbit — the registry holds canonical-sorted targets but the lookup is
	 * permutation-invariant). Returned identities are guaranteed to verify
	 * against the live catalog (their predicted rank arithmetic is valid).
	 *
	 * <p>{@link eu.solven.matmul.recombination.BlockSplitSearch#findBestStrategy}
	 * iterates this list at every target shape and offers each as a
	 * candidate strategy alongside Kron / Concat / recombination.</p>
	 */
	public static List<Identity> applicableTo(int n, int m, int p,
			eu.solven.matmul.recombination.Recombination.SotaResolver sota) {
		int[] sorted = { n, m, p };
		java.util.Arrays.sort(sorted);
		List<Identity> out = new java.util.ArrayList<>();
		for (Identity id : ALL) {
			int[] t = id.target().clone();
			java.util.Arrays.sort(t);
			if (t[0] != sorted[0] || t[1] != sorted[1] || t[2] != sorted[2]) continue;
			if (predictedRankAgainstSota(id, sota) == id.predictedRank()) {
				out.add(id);
			}
		}
		return out;
	}

	/**
	 * Recompute the identity's predicted rank using the supplied SOTA
	 * resolver (catalog ranks at the moment of invocation). Returns
	 * {@code -1} if any sub-shape can't be resolved.
	 */
	public static long predictedRankAgainstSota(Identity id,
			eu.solven.matmul.recombination.Recombination.SotaResolver sota) {
		if (id instanceof DisjointSum ds) {
			long total = 0;
			for (SubTerm t : ds.terms) {
				long termCost;
				if (t.kind == SubTerm.Kind.PAN_PAIR_COST) {
					termCost = (long) t.n * t.m * t.p
							+ (long) t.n * t.m + (long) t.m * t.p + (long) t.p * t.n;
				} else {
					int r = sota.getRank(t.n, t.m, t.p);
					if (r >= Recombination.SotaResolver.UNKNOWN_RANK) return -1;
					termCost = r;
				}
				total += (long) t.multiplicity * termCost;
			}
			return total;
		} else if (id instanceof BorrowAndCorrect bc) {
			int rA = sota.getRank(bc.shapeA[0], bc.shapeA[1], bc.shapeA[2]);
			int rB = sota.getRank(bc.shapeB[0], bc.shapeB[1], bc.shapeB[2]);
			int rC = sota.getRank(bc.shapeCaux[0], bc.shapeCaux[1], bc.shapeCaux[2]);
			if (rA >= Recombination.SotaResolver.UNKNOWN_RANK || rB >= Recombination.SotaResolver.UNKNOWN_RANK
					|| rC >= Recombination.SotaResolver.UNKNOWN_RANK) return -1;
			return (long) (rA - bc.discountK) * rB + (long) bc.auxMultiplicity * rC;
		}
		return -1;
	}

	/**
	 * Verify each identity against the live catalog: assert that the
	 * predicted rank arithmetic uses the catalog's current ranks for the
	 * referenced sub-shapes. Identities that fail verification are
	 * returned so the catalog generator can flag the regression.
	 *
	 * @param lookup catalog rank lookup
	 * @return list of identities whose recipe arithmetic no longer matches
	 *         the catalog (e.g. a sub-shape improved and the recipe became
	 *         a non-tight upper bound, or worsened and the recipe is now
	 *         infeasible)
	 */
	public static List<Identity> verify(FieldAwareLookup lookup) {
		List<Identity> failed = new java.util.ArrayList<>();
		for (Identity id : ALL) {
			long computed;
			if (id instanceof DisjointSum ds) {
				computed = 0;
				for (SubTerm t : ds.terms) {
					long termCost;
					if (t.kind == SubTerm.Kind.PAN_PAIR_COST) {
						// Pan 1980 pair-fusion cost for two cyclically-related
						// products: abc + ab + bc + ca.
						termCost = (long) t.n * t.m * t.p
								+ (long) t.n * t.m + (long) t.m * t.p + (long) t.p * t.n;
					} else {
						int r = lookup.findRank(t.n, t.m, t.p);
						if (r >= Recombination.SotaResolver.UNKNOWN_RANK) {
							failed.add(id);
							computed = -1;
							break;
						}
						termCost = r;
					}
					computed += (long) t.multiplicity * termCost;
				}
				if (computed > 0 && computed != ds.predictedRank) failed.add(id);
			} else if (id instanceof BorrowAndCorrect bc) {
				int rA = lookup.findRank(bc.shapeA[0], bc.shapeA[1], bc.shapeA[2]);
				int rB = lookup.findRank(bc.shapeB[0], bc.shapeB[1], bc.shapeB[2]);
				int rC = lookup.findRank(bc.shapeCaux[0], bc.shapeCaux[1], bc.shapeCaux[2]);
				if (rA >= Recombination.SotaResolver.UNKNOWN_RANK || rB >= Recombination.SotaResolver.UNKNOWN_RANK
						|| rC >= Recombination.SotaResolver.UNKNOWN_RANK) {
					failed.add(id);
					continue;
				}
				computed = (long) (rA - bc.discountK) * rB + (long) bc.auxMultiplicity * rC;
				if (computed != bc.predictedRank) failed.add(id);
			}
		}
		return failed;
	}
}
