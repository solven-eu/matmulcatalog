package eu.solven.matmul.catalog;

import eu.solven.matmul.papers.hopcroftkerr1971.HopcroftKerr2bcAsymmetric;

import java.util.ArrayList;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Map;

/**
 * Structured construction lineage of a materialised scheme. Captures
 * "how was this composite built from leaves" as a {@link Node} DAG with
 * named operations matching {@link Compose} / {@link Recombination}
 * primitives.
 *
 * <p>Two complementary serialisations are supported:</p>
 * <ul>
 *   <li>{@link #prettyString(Node)} — human-readable function-call form,
 *       e.g. {@code ConcatCols(KronProduct(perminov_Z-3x6x3_m40_a194,
 *       fmm_lille-6x3x3_m40), SAME0)}. Used in {@code lineage_str}
 *       field of the scheme JSON. Repeated subtrees are tagged
 *       {@code SAME<i>} on second/later occurrence.</li>
 *   <li>{@link #toJson(Node)} — JSON DAG with explicit
 *       {@code {"op": ..., "args": ...}} nodes and {@code @ref}
 *       back-references for shared subtrees. Used in {@code lineage}
 *       field. Machine-readable; can be fed back to a future
 *       {@code MaterializeFromLineage} to re-derive the scheme.</li>
 * </ul>
 *
 * <p>Operation vocabulary — each entry mirrors a primitive in the
 * composition code base, so lineage and code stay in lockstep:</p>
 *
 * <table>
 *   <caption>Lineage operations</caption>
 *   <tr><th>op</th><th>shape</th><th>source primitive</th></tr>
 *   <tr><td>{@code Atom}</td><td>{@code Atom(ref)}</td>
 *       <td>filename stem of a JSON on disk (e.g. "fmm_lille-2x3x4_m11_a45")</td></tr>
 *   <tr><td>{@code KronProduct}</td><td>{@code KronProduct(outer, inner)}</td>
 *       <td>{@link Compose#kroneckerGeneral}</td></tr>
 *   <tr><td>{@code KronChain}</td><td>{@code KronChain(L1, L2, ...)}</td>
 *       <td>folded k-way Kronecker</td></tr>
 *   <tr><td>{@code ConcatCols}</td><td>{@code ConcatCols(left, right)}</td>
 *       <td>{@link Compose#concatRight} (p-axis tile; was {@code ConcatRight})</td></tr>
 *   <tr><td>{@code ConcatRows}</td><td>{@code ConcatRows(top, bottom)}</td>
 *       <td>{@link Compose#concatBelow} (n-axis tile; was {@code ConcatBelow})</td></tr>
 *   <tr><td>{@code SumInner}</td><td>{@code SumInner(left, right)}</td>
 *       <td>{@link Compose#concatInner} (m-axis contraction sum, C=C1+C2)</td></tr>
 *   <tr><td>{@code Recombination}</td>
 *       <td>{@code Recombination(base, allocA, allocB, allocC, leaves)}</td>
 *       <td>{@link Recombination#constructWithAllocation}</td></tr>
 *   <tr><td>{@code RecombinationWithPair}</td>
 *       <td>{@code RecombinationWithPair(base, pairing, leaves, pairs)}</td>
 *       <td>{@code RecombinationWithPair.constructWithPairing}</td></tr>
 *   <tr><td>{@code AugmentSquareDiscard}</td>
 *       <td>{@code AugmentSquareDiscard(p, n, square)}</td>
 *       <td>{@code HopcroftKerr2bcAsymmetric.buildNaive}</td></tr>
 *   <tr><td>{@code DCE}</td><td>{@code DCE(child)}</td>
 *       <td>dead-code-elimination post-pass</td></tr>
 *   <tr><td>{@code Transpose}</td>
 *       <td>{@code Transpose(child, perm)} (perm = "ABC->CAB" etc.)</td>
 *       <td>tensor-symmetry rewrite</td></tr>
 * </table>
 *
 * <p>Adding a new op: extend the {@link Node} sealed hierarchy and the
 * pretty/JSON renderers below. Do <em>not</em> repurpose an existing op
 * — leave the vocabulary aligned with the primitives so lineage stays
 * reversible.</p>
 */
public final class Lineage {

	private Lineage() {}

	/** Sealed node hierarchy. */
	public sealed interface Node permits Atom, KronProduct, KronChain, ConcatCols,
			ConcatRows, SumInner, RecombinationN, RecombinationTaN, RecombinationWithPairN,
			AugmentSquareDiscard, Dce, Transpose, AxisFlip, AxisPermute, DisjointSum,
			SerendipitousProduct, Project, OrientAs, PeeledViaTa {}

	/**
	 * Peeled matmul {@code ⟨n, n+s, n+s⟩} assembled via rectangular Pan
	 * trilinear aggregation: {@code diag ⟨n,n,n⟩ (=cube) + TA-fused cross-pair
	 * (⟨n,s,n⟩⊕⟨n,n,s⟩) + corner ⟨n,s,s⟩}. Replays via
	 * {@code RectangularTrilinearAggregation.buildPeeledViaTa(n, s, replay(cube),
	 * replay(corner))} — the TA part is deterministic in {@code (n,s)} (integer
	 * ±1 coefficients), so {@code (n,s,cube,corner)} is sufficient for bit-exact
	 * replay. See {@code references/RECTANGULAR_TA.md}.
	 */
	public record PeeledViaTa(int n, int s, Node cube, Node corner) implements Node {}

	/** Atom: reference to a primitive on-disk scheme by its filename stem (no .json). */
	public record Atom(String ref) implements Node {}

	public record KronProduct(Node outer, Node inner) implements Node {}

	public record KronChain(List<Node> factors) implements Node {}

	/**
	 * p-axis (output-column) tile: {@code ⟨n,m,p1⟩ + ⟨n,m,p2⟩ → ⟨n,m,p1+p2⟩},
	 * sharing A, stacking C's columns {@code [C1|C2]}.
	 * Source primitive: {@link Compose#concatRight}. (Formerly {@code ConcatRight}.)
	 */
	public record ConcatCols(Node left, Node right) implements Node {}

	/**
	 * n-axis (output-row) tile: {@code ⟨n1,m,p⟩ + ⟨n2,m,p⟩ → ⟨n1+n2,m,p⟩},
	 * sharing B, stacking C's rows {@code [C1;C2]}.
	 * Source primitive: {@link Compose#concatBelow}. (Formerly {@code ConcatBelow}.)
	 */
	public record ConcatRows(Node top, Node bottom) implements Node {}

	/**
	 * m-axis (inner/contraction) sum: {@code ⟨n,m1,p⟩ + ⟨n,m2,p⟩ → ⟨n,m1+m2,p⟩},
	 * sharing neither operand, <em>accumulating</em> {@code C = C1 + C2}
	 * (not tiling). Source primitive: {@link Compose#concatInner}. This is
	 * the third sibling completing the per-mode direct-sum set; it is a
	 * sum, not a geometric concat, hence the distinct name.
	 */
	public record SumInner(Node left, Node right) implements Node {}

	public record RecombinationN(Node base, int[] allocA, int[] allocB,
			int[] allocC, List<Node> leaves) implements Node {}

	/**
	 * Generic Pan-TA-fused recombination (the generalisation of {@link PeeledViaTa}
	 * to ANY naïve-grid base). Records only {@code base + allocations}; replay re-runs
	 * {@link Recombination#constructWithTaFusion} — which deterministically re-derives
	 * the disjoint cyclic-rotation product pairs to TA-fuse and re-resolves the unpaired
	 * leaves to the catalog NON-COMMUTATIVE best by shape. No leaf list is stored (the
	 * leaves are a pure function of base+allocs+catalog); only {@code base} is pinned.
	 * Used for FMM-style block grids (e.g. ⟨2,3,3⟩) where multiple cross-pairs fuse.
	 *
	 * <p>{@code leaves} are the (NON-COMMUTATIVE-pinned) schemes for the UNPAIRED
	 * products — one per DISTINCT unpaired shape — so replay resolves them exactly
	 * (not via the commutative-blind shape-best) and field inference can read them.
	 * The fused products carry no leaf (the TA block's coefficients are integer ±1,
	 * field-neutral).</p>
	 */
	public record RecombinationTaN(Node base, int[] allocA, int[] allocB,
			int[] allocC, List<Node> leaves) implements Node {}

	public record RecombinationWithPairN(Node base, int[][] pairs, int[] solo,
			List<Node> leaves) implements Node {}

	public record AugmentSquareDiscard(int p, int n, Node square) implements Node {}

	public record Dce(Node child) implements Node {}

	public record Transpose(Node child, String perm) implements Node {}

	/**
	 * Orient the child scheme to the exact shape {@code ⟨n,m,p⟩} (an S₃ axis
	 * permutation of the matmul tensor). Unlike {@link Transpose}'s string perm,
	 * this records the <em>target dims directly</em>, so replay is unambiguous even
	 * when dims repeat: {@code replay(child).orientAs(n,m,p)}. Emitted by the
	 * materialiser whenever it reuses a sub-scheme in a different orientation than
	 * its stored/native one, so {@code replay(lineage)} reproduces the EXACT
	 * orientation of the returned algorithm (the replayable-lineage invariant).
	 */
	public record OrientAs(Node child, int n, int m, int p) implements Node {}

	/**
	 * Axis-flip orbit transform: {@code child} is the canonical scheme;
	 * {@code mask} is the 3-bit axis-flip mask
	 * ({@code bit 0 = swapA, bit 1 = swapB, bit 2 = swapC}). The
	 * factor matrices of the resulting scheme are obtained by reversing
	 * the index order on the indicated axes of {@code child}'s
	 * {@code U, V, W} per
	 * {@link eu.solven.matmul.SymmetryTransforms#axisFlipOrbit}.
	 *
	 * <p>Marks the scheme as "produced by an algorithmic-orbit
	 * transformation, not the canonical leaf". A LineageReplayer that
	 * understands axis-flip can recompute the variant from
	 * {@code child}; a simpler replayer that doesn't can still load the
	 * scheme's on-disk {@code (U, V, W)} directly and treat this lineage
	 * as informational. The whole-scheme materialisation is intentional
	 * — see {@code docs/STRATEGIES.md} §4.1 for the audit-vs-default
	 * discussion.</p>
	 */
	public record AxisFlip(Node child, int mask) implements Node {
		public AxisFlip {
			if (mask < 0 || mask > 7) {
				throw new IllegalArgumentException(
						"axis-flip mask must be in [0,7], got " + mask);
			}
		}
	}

	/**
	 * Per-axis permutation orbit transform: {@code child}'s factor
	 * matrices are obtained by permuting the rows of A/B/C per
	 * {@code permA, permB, permC} respectively (each a permutation of
	 * {@code [0, n)}, {@code [0, m)}, {@code [0, p)}). See
	 * {@link eu.solven.matmul.SymmetryTransforms#permutationOrbit}.
	 *
	 * <p>More expressive than {@link AxisFlip} (which is the
	 * permutation subgroup generated by the axis-reversal {@code J}).</p>
	 */
	public record AxisPermute(Node child, int[] permA, int[] permB, int[] permC)
			implements Node {}

	/**
	 * Serendipitous product (Smith 2002 eq. (69); re-derived Perminov draft
	 * Def 2.12; Sedoglavic):
	 * {@code base ⊗ ⟨n2,m2,p2⟩} where {@code base} is decomposed into elementary
	 * matmul tensors by its <em>buds</em> (rank-one terms sharing a {@code u}/
	 * {@code v}/{@code w} vector up to scaling) and each enlarged block is
	 * realised at its best known rank. Replay recomputes the (deterministic
	 * greedy U→V→W) bud decomposition of {@code base}, so only the base and the
	 * second shape need to be stored. See {@link SerendipitousBudProduct} and
	 * {@code references/SERENDIPITOUS_PARTIAL_PRODUCT.md}.
	 */
	public record SerendipitousProduct(Node base, int n2, int m2, int p2) implements Node {}

	/**
	 * Projection (Perminov draft Def 2.8 / meta-flip-graph {@code Project};
	 * = padding+DCE in reverse). The parent {@code child} is restricted to the
	 * kept indices {@code keepN/keepM/keepP} on its three axes, yielding a smaller
	 * matmul scheme; products whose restricted factor is all-zero are
	 * dead-code-eliminated. Replay = {@link Compose#project}. Which indices are
	 * kept (equivalently, which are dropped) determines the resulting rank.
	 */
	public record Project(Node child, int[] keepN, int[] keepM, int[] keepP) implements Node {}

	/**
	 * τ-theorem-style disjoint-MM-sum decomposition (Pan 1980 / Schönhage
	 * 1981; DIS09 §3; Schwartz-Zwecher 2025). The target matmul is
	 * computed as a sum of smaller matmul sub-tensors, optionally with
	 * "trilinear aggregation" legs that compute several same-cubic-shape
	 * sub-products jointly via Pan's pair-cost formula
	 * {@code abc + ab + bc + ca}.
	 *
	 * <p>{@code children} is the list of sub-MM nodes (each one a
	 * placeholder leaf or a deeper composition); {@code taLegs} groups
	 * indices into {@code children} that are computed jointly via TA
	 * (each TA leg is a list of indices, all sharing the same cubic
	 * shape).</p>
	 *
	 * <p><strong>Materialisation (factor matrices) is a follow-up.</strong>
	 * For now this node is a rank-prediction marker: it encodes the
	 * shape-decomposition + TA structure, and the predicted total rank is
	 * computed by callers using a SOTA lookup.</p>
	 */
	public record DisjointSum(java.util.List<Node> children,
			java.util.List<java.util.List<Integer>> taLegs) implements Node {}

	/**
	 * Render the lineage as a human-readable function-call expression.
	 * Repeated subtrees are <strong>plain-repeated</strong> (no aliasing)
	 * — the string is intended to be read end-to-end without cross-refs,
	 * e.g. shown in SPA tooltips. Use {@link #toJson(Node)} when you
	 * need a DAG that dedups shared subtrees by id.
	 */
	public static String prettyString(Node root) {
		StringBuilder sb = new StringBuilder();
		renderPretty(root, sb);
		return sb.toString();
	}

	private static void renderPretty(Node n, StringBuilder sb) {
		switch (n) {
			case Atom l -> sb.append(shortRefForDisplay(l.ref));
			case KronProduct kp -> {
				sb.append("KronProduct(");
				renderPretty(kp.outer, sb);
				sb.append(", ");
				renderPretty(kp.inner, sb);
				sb.append(")");
			}
			case KronChain kc -> {
				sb.append("KronChain(");
				for (int i = 0; i < kc.factors.size(); i++) {
					if (i > 0) sb.append(", ");
					renderPretty(kc.factors.get(i), sb);
				}
				sb.append(")");
			}
			case ConcatCols cr -> {
				sb.append("ConcatCols(");
				renderPretty(cr.left, sb);
				sb.append(", ");
				renderPretty(cr.right, sb);
				sb.append(")");
			}
			case PeeledViaTa t -> {
				sb.append("PeeledViaTa(n=").append(t.n).append(",s=").append(t.s).append(", cube=");
				renderPretty(t.cube, sb);
				sb.append(", corner=");
				renderPretty(t.corner, sb);
				sb.append(")");
			}
			case ConcatRows cb -> {
				sb.append("ConcatRows(");
				renderPretty(cb.top, sb);
				sb.append(", ");
				renderPretty(cb.bottom, sb);
				sb.append(")");
			}
			case SumInner si -> {
				sb.append("SumInner(");
				renderPretty(si.left, sb);
				sb.append(", ");
				renderPretty(si.right, sb);
				sb.append(")");
			}
			case RecombinationN r -> {
				sb.append("Recombination(base=");
				renderPretty(r.base, sb);
				sb.append(", allocA=").append(java.util.Arrays.toString(r.allocA));
				sb.append(", allocB=").append(java.util.Arrays.toString(r.allocB));
				sb.append(", allocC=").append(java.util.Arrays.toString(r.allocC));
				sb.append(", leaves=[");
				for (int i = 0; i < r.leaves.size(); i++) {
					if (i > 0) sb.append(", ");
					renderPretty(r.leaves.get(i), sb);
				}
				sb.append("])");
			}
			case RecombinationTaN r -> {
				sb.append("RecombinationTa(base=");
				renderPretty(r.base, sb);
				sb.append(", allocA=").append(java.util.Arrays.toString(r.allocA));
				sb.append(", allocB=").append(java.util.Arrays.toString(r.allocB));
				sb.append(", allocC=").append(java.util.Arrays.toString(r.allocC));
				sb.append(", leaves=[");
				for (int i = 0; i < r.leaves.size(); i++) {
					if (i > 0) sb.append(", ");
					renderPretty(r.leaves.get(i), sb);
				}
				sb.append("])");
			}
			case RecombinationWithPairN r -> {
				sb.append("RecombinationWithPair(base=");
				renderPretty(r.base, sb);
				sb.append(", pairs=").append(java.util.Arrays.deepToString(r.pairs));
				sb.append(", solo=").append(java.util.Arrays.toString(r.solo));
				sb.append(", leaves=[");
				for (int i = 0; i < r.leaves.size(); i++) {
					if (i > 0) sb.append(", ");
					renderPretty(r.leaves.get(i), sb);
				}
				sb.append("])");
			}
			case AugmentSquareDiscard a -> {
				sb.append("AugmentSquareDiscard(p=").append(a.p).append(", n=").append(a.n).append(", square=");
				renderPretty(a.square, sb);
				sb.append(")");
			}
			case Dce d -> {
				sb.append("DCE(");
				renderPretty(d.child, sb);
				sb.append(")");
			}
			case Transpose t -> {
				sb.append("Transpose(");
				renderPretty(t.child, sb);
				sb.append(", perm=\"").append(t.perm).append("\")");
			}
			case OrientAs o -> {
				sb.append("OrientAs(");
				renderPretty(o.child, sb);
				sb.append(", ⟨").append(o.n).append(",").append(o.m).append(",").append(o.p).append("⟩)");
			}
			case AxisFlip af -> {
				sb.append("AxisFlip[approximate](");
				renderPretty(af.child, sb);
				sb.append(", mask=").append(af.mask).append(")");
			}
			case AxisPermute ap -> {
				sb.append("AxisPermute[approximate](");
				renderPretty(ap.child, sb);
				sb.append(", permA=").append(java.util.Arrays.toString(ap.permA));
				sb.append(", permB=").append(java.util.Arrays.toString(ap.permB));
				sb.append(", permC=").append(java.util.Arrays.toString(ap.permC));
				sb.append(")");
			}
			case DisjointSum ds -> {
				sb.append("DisjointSum(");
				for (int i = 0; i < ds.children.size(); i++) {
					if (i > 0) sb.append(" + ");
					renderPretty(ds.children.get(i), sb);
				}
				if (!ds.taLegs.isEmpty()) {
					sb.append("; TA-legs=").append(ds.taLegs);
				}
				sb.append(")");
			}
			case SerendipitousProduct sp -> {
				sb.append("Serendipitous(");
				renderPretty(sp.base, sb);
				sb.append(" ⊗ ⟨").append(sp.n2).append(",").append(sp.m2).append(",")
						.append(sp.p2).append("⟩)");
			}
			case Project pr -> {
				sb.append("Project(");
				renderPretty(pr.child, sb);
				sb.append(", keepN=").append(java.util.Arrays.toString(pr.keepN));
				sb.append(", keepM=").append(java.util.Arrays.toString(pr.keepM));
				sb.append(", keepP=").append(java.util.Arrays.toString(pr.keepP)).append(")");
			}
		}
	}

	/**
	 * Compact one-line rendering of the lineage, optimised for size
	 * (folded into {@code catalog.json} alongside per-entry metadata) and
	 * for at-a-glance reading. The grammar:
	 *
	 * <pre>
	 *   Atom(ref)                    → ref                            (drop wrapper)
	 *   Recombination(B, aA, aB, aC, leaves)
	 *                                → R[B; aA | aB | aC]             (DROP leaves —
	 *                                                                  derivable from
	 *                                                                  B + allocs)
	 *   RecombinationWithPair(...)   → R*[B; aA | aB | aC; pairs/solo]
	 *   KronProduct(A, B)            → A ⊗ B
	 *   KronChain(L1,…,Lk)           → L1 ⊗ L2 ⊗ … ⊗ Lk
	 *   ConcatCols(L, R)             → L +p R
	 *   ConcatRows(T, B)             → T +n B
	 *   SumInner(L, R)               → L +m R   (contraction sum, C=C1+C2)
	 *   Transpose(child, perm)       → child^perm
	 *   AxisFlip(child, mask)        → child^J&lt;mask&gt;
	 *   AxisPermute(...)             → child^π(…)
	 *   DisjointSum([c1,…,ck], …)    → c1 ⊕ c2 ⊕ … ⊕ ck
	 *   Dce(child)                   → child                          (elide)
	 * </pre>
	 *
	 * <p>Base aliases are kept as plain strings (no single-letter codes) —
	 * "Makarov" is ambiguous between 3x3x3 and 5x5x5 results, so the full
	 * scheme reference stays in the lineage.
	 *
	 * <p><b>Leaves are dropped from Recombination.</b> They're recoverable
	 * by re-running {@code findBestStrategy} with the recorded
	 * {@code (base, allocations)} tuple and the catalog at replay time;
	 * the trade is ~3× disk savings vs a one-time replay-lookup cost. A
	 * local replay cache can absorb that cost — tracked separately.
	 *
	 * <p>Like {@link #prettyString}, repeated subtrees are plain-repeated
	 * (no aliasing); use {@link #toJson} when DAG dedup matters.
	 */
	public static String prettyCompact(Node root) {
		StringBuilder sb = new StringBuilder();
		renderCompact(root, sb);
		return sb.toString();
	}

	private static void renderCompact(Node n, StringBuilder sb) {
		switch (n) {
			case Atom l -> sb.append(shortRefForDisplay(l.ref));
			case KronProduct kp -> {
				renderCompact(kp.outer, sb);
				sb.append(" ⊗ ");
				renderCompact(kp.inner, sb);
			}
			case KronChain kc -> {
				for (int i = 0; i < kc.factors.size(); i++) {
					if (i > 0) sb.append(" ⊗ ");
					renderCompact(kc.factors.get(i), sb);
				}
			}
			case ConcatCols cr -> {
				renderCompact(cr.left, sb);
				sb.append(" +p ");
				renderCompact(cr.right, sb);
			}
			case PeeledViaTa t -> {
				sb.append("TA[n=").append(t.n).append(",s=").append(t.s).append("; ");
				renderCompact(t.cube, sb);
				sb.append(" ; ");
				renderCompact(t.corner, sb);
				sb.append("]");
			}
			case ConcatRows cb -> {
				renderCompact(cb.top, sb);
				sb.append(" +n ");
				renderCompact(cb.bottom, sb);
			}
			case SumInner si -> {
				renderCompact(si.left, sb);
				sb.append(" +m ");
				renderCompact(si.right, sb);
			}
			case RecombinationN r -> {
				sb.append("R[");
				renderCompact(r.base, sb);
				sb.append("; ").append(joinInts(r.allocA));
				sb.append(" | ").append(joinInts(r.allocB));
				sb.append(" | ").append(joinInts(r.allocC));
				sb.append("]");
			}
			case RecombinationTaN r -> {
				sb.append("Rta[");
				renderCompact(r.base, sb);
				sb.append("; ").append(joinInts(r.allocA));
				sb.append(" | ").append(joinInts(r.allocB));
				sb.append(" | ").append(joinInts(r.allocC));
				sb.append("]");
			}
			case RecombinationWithPairN r -> {
				sb.append("R*[");
				renderCompact(r.base, sb);
				sb.append("; pairs=").append(java.util.Arrays.deepToString(r.pairs));
				sb.append(", solo=").append(java.util.Arrays.toString(r.solo));
				sb.append("]");
			}
			case AugmentSquareDiscard a -> {
				sb.append("AS(").append(a.p).append(",").append(a.n).append(",");
				renderCompact(a.square, sb);
				sb.append(")");
			}
			case Dce d -> renderCompact(d.child, sb);  // elide
			case Transpose t -> {
				renderCompact(t.child, sb);
				sb.append("^").append(t.perm);
			}
			case OrientAs o -> {
				renderCompact(o.child, sb);
				sb.append("→⟨").append(o.n).append(",").append(o.m).append(",").append(o.p).append("⟩");
			}
			case AxisFlip af -> {
				renderCompact(af.child, sb);
				sb.append("^J<").append(af.mask).append(">");
			}
			case AxisPermute ap -> {
				renderCompact(ap.child, sb);
				sb.append("^π(").append(java.util.Arrays.toString(ap.permA));
				sb.append(",").append(java.util.Arrays.toString(ap.permB));
				sb.append(",").append(java.util.Arrays.toString(ap.permC));
				sb.append(")");
			}
			case DisjointSum ds -> {
				for (int i = 0; i < ds.children.size(); i++) {
					if (i > 0) sb.append(" ⊕ ");
					renderCompact(ds.children.get(i), sb);
				}
				if (!ds.taLegs.isEmpty()) {
					sb.append("; TA=").append(ds.taLegs);
				}
			}
			case SerendipitousProduct sp -> {
				renderCompact(sp.base, sb);
				sb.append(" ⊗ˢ⟨").append(sp.n2).append(",").append(sp.m2).append(",")
						.append(sp.p2).append("⟩");
			}
			case Project pr -> {
				renderCompact(pr.child, sb);
				sb.append(" ↓[").append(joinInts(pr.keepN)).append("|")
						.append(joinInts(pr.keepM)).append("|")
						.append(joinInts(pr.keepP)).append("]");
			}
		}
	}

	private static String joinInts(int[] a) {
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < a.length; i++) {
			if (i > 0) sb.append(",");
			sb.append(a[i]);
		}
		return sb.toString();
	}

	/**
	 * Render the lineage as a JSON DAG. Each non-leaf node becomes
	 * {@code {"op": "...", ...}}; repeated nodes after the first
	 * occurrence get {@code {"op": "@ref", "id": "L<i>"}}, with the
	 * first occurrence tagged {@code "id": "L<i>"} in addition to its
	 * op fields.
	 */
	public static String toJson(Node root) {
		IdentityHashMap<Node, Integer> seen = countOccurrences(root);
		Map<Node, Integer> sameIds = new IdentityHashMap<>();
		int[] nextId = { 0 };
		StringBuilder sb = new StringBuilder();
		renderJson(root, sb, seen, sameIds, nextId);
		return sb.toString();
	}

	private static void renderJson(Node n, StringBuilder sb,
			IdentityHashMap<Node, Integer> seen,
			Map<Node, Integer> sameIds, int[] nextId) {
		if (seen.getOrDefault(n, 0) > 1 && sameIds.containsKey(n)) {
			sb.append("{\"op\":\"@ref\",\"id\":\"L").append(sameIds.get(n)).append("\"}");
			return;
		}
		String idBinding = "";
		if (seen.getOrDefault(n, 0) > 1) {
			int id = nextId[0]++;
			sameIds.put(n, id);
			idBinding = ",\"id\":\"L" + id + "\"";
		}
		switch (n) {
			case Atom l ->
				sb.append("{\"op\":\"Atom\",\"ref\":").append(jsonStr(l.ref)).append(idBinding).append("}");
			case KronProduct kp -> {
				sb.append("{\"op\":\"KronProduct\"").append(idBinding).append(",\"outer\":");
				renderJson(kp.outer, sb, seen, sameIds, nextId);
				sb.append(",\"inner\":");
				renderJson(kp.inner, sb, seen, sameIds, nextId);
				sb.append("}");
			}
			case KronChain kc -> {
				sb.append("{\"op\":\"KronChain\"").append(idBinding).append(",\"factors\":[");
				for (int i = 0; i < kc.factors.size(); i++) {
					if (i > 0) sb.append(",");
					renderJson(kc.factors.get(i), sb, seen, sameIds, nextId);
				}
				sb.append("]}");
			}
			case ConcatCols cr -> {
				sb.append("{\"op\":\"ConcatCols\"").append(idBinding).append(",\"left\":");
				renderJson(cr.left, sb, seen, sameIds, nextId);
				sb.append(",\"right\":");
				renderJson(cr.right, sb, seen, sameIds, nextId);
				sb.append("}");
			}
			case ConcatRows cb -> {
				sb.append("{\"op\":\"ConcatRows\"").append(idBinding).append(",\"top\":");
				renderJson(cb.top, sb, seen, sameIds, nextId);
				sb.append(",\"bottom\":");
				renderJson(cb.bottom, sb, seen, sameIds, nextId);
				sb.append("}");
			}
			case SumInner si -> {
				sb.append("{\"op\":\"SumInner\"").append(idBinding).append(",\"left\":");
				renderJson(si.left, sb, seen, sameIds, nextId);
				sb.append(",\"right\":");
				renderJson(si.right, sb, seen, sameIds, nextId);
				sb.append("}");
			}
			case PeeledViaTa t -> {
				sb.append("{\"op\":\"PeeledViaTa\"").append(idBinding)
						.append(",\"n\":").append(t.n).append(",\"s\":").append(t.s).append(",\"cube\":");
				renderJson(t.cube, sb, seen, sameIds, nextId);
				sb.append(",\"corner\":");
				renderJson(t.corner, sb, seen, sameIds, nextId);
				sb.append("}");
			}
			case RecombinationN r -> {
				sb.append("{\"op\":\"Recombination\"").append(idBinding).append(",\"base\":");
				renderJson(r.base, sb, seen, sameIds, nextId);
				sb.append(",\"allocA\":").append(java.util.Arrays.toString(r.allocA));
				sb.append(",\"allocB\":").append(java.util.Arrays.toString(r.allocB));
				sb.append(",\"allocC\":").append(java.util.Arrays.toString(r.allocC));
				sb.append(",\"leaves\":[");
				for (int i = 0; i < r.leaves.size(); i++) {
					if (i > 0) sb.append(",");
					renderJson(r.leaves.get(i), sb, seen, sameIds, nextId);
				}
				sb.append("]}");
			}
			case RecombinationTaN r -> {
				sb.append("{\"op\":\"RecombinationTa\"").append(idBinding).append(",\"base\":");
				renderJson(r.base, sb, seen, sameIds, nextId);
				sb.append(",\"allocA\":").append(java.util.Arrays.toString(r.allocA));
				sb.append(",\"allocB\":").append(java.util.Arrays.toString(r.allocB));
				sb.append(",\"allocC\":").append(java.util.Arrays.toString(r.allocC));
				sb.append(",\"leaves\":[");
				for (int i = 0; i < r.leaves.size(); i++) {
					if (i > 0) sb.append(",");
					renderJson(r.leaves.get(i), sb, seen, sameIds, nextId);
				}
				sb.append("]}");
			}
			case RecombinationWithPairN r -> {
				sb.append("{\"op\":\"RecombinationWithPair\"").append(idBinding).append(",\"base\":");
				renderJson(r.base, sb, seen, sameIds, nextId);
				sb.append(",\"pairs\":").append(java.util.Arrays.deepToString(r.pairs));
				sb.append(",\"solo\":").append(java.util.Arrays.toString(r.solo));
				sb.append(",\"leaves\":[");
				for (int i = 0; i < r.leaves.size(); i++) {
					if (i > 0) sb.append(",");
					renderJson(r.leaves.get(i), sb, seen, sameIds, nextId);
				}
				sb.append("]}");
			}
			case AugmentSquareDiscard a -> {
				sb.append("{\"op\":\"AugmentSquareDiscard\"").append(idBinding)
						.append(",\"p\":").append(a.p).append(",\"n\":").append(a.n).append(",\"square\":");
				renderJson(a.square, sb, seen, sameIds, nextId);
				sb.append("}");
			}
			case Dce d -> {
				sb.append("{\"op\":\"DCE\"").append(idBinding).append(",\"child\":");
				renderJson(d.child, sb, seen, sameIds, nextId);
				sb.append("}");
			}
			case Transpose t -> {
				sb.append("{\"op\":\"Transpose\"").append(idBinding).append(",\"perm\":")
						.append(jsonStr(t.perm)).append(",\"child\":");
				renderJson(t.child, sb, seen, sameIds, nextId);
				sb.append("}");
			}
			case OrientAs o -> {
				sb.append("{\"op\":\"OrientAs\"").append(idBinding)
						.append(",\"n\":").append(o.n).append(",\"m\":").append(o.m)
						.append(",\"p\":").append(o.p).append(",\"child\":");
				renderJson(o.child, sb, seen, sameIds, nextId);
				sb.append("}");
			}
			case AxisFlip af -> {
				sb.append("{\"op\":\"AxisFlip\"").append(idBinding)
						.append(",\"mask\":").append(af.mask)
						.append(",\"approximate\":true,\"child\":");
				renderJson(af.child, sb, seen, sameIds, nextId);
				sb.append("}");
			}
			case AxisPermute ap -> {
				sb.append("{\"op\":\"AxisPermute\"").append(idBinding)
						.append(",\"permA\":").append(java.util.Arrays.toString(ap.permA))
						.append(",\"permB\":").append(java.util.Arrays.toString(ap.permB))
						.append(",\"permC\":").append(java.util.Arrays.toString(ap.permC))
						.append(",\"approximate\":true,\"child\":");
				renderJson(ap.child, sb, seen, sameIds, nextId);
				sb.append("}");
			}
			case DisjointSum ds -> {
				sb.append("{\"op\":\"DisjointSum\"").append(idBinding).append(",\"children\":[");
				for (int i = 0; i < ds.children.size(); i++) {
					if (i > 0) sb.append(",");
					renderJson(ds.children.get(i), sb, seen, sameIds, nextId);
				}
				sb.append("],\"taLegs\":").append(ds.taLegs).append("}");
			}
			case SerendipitousProduct sp -> {
				sb.append("{\"op\":\"SerendipitousProduct\"").append(idBinding)
						.append(",\"n2\":").append(sp.n2).append(",\"m2\":").append(sp.m2)
						.append(",\"p2\":").append(sp.p2).append(",\"base\":");
				renderJson(sp.base, sb, seen, sameIds, nextId);
				sb.append("}");
			}
			case Project pr -> {
				sb.append("{\"op\":\"Project\"").append(idBinding)
						.append(",\"keepN\":[").append(joinInts(pr.keepN)).append("]")
						.append(",\"keepM\":[").append(joinInts(pr.keepM)).append("]")
						.append(",\"keepP\":[").append(joinInts(pr.keepP)).append("]")
						.append(",\"child\":");
				renderJson(pr.child, sb, seen, sameIds, nextId);
				sb.append("}");
			}
		}
	}

	private static IdentityHashMap<Node, Integer> countOccurrences(Node root) {
		IdentityHashMap<Node, Integer> counts = new IdentityHashMap<>();
		walk(root, counts);
		return counts;
	}

	private static void walk(Node n, IdentityHashMap<Node, Integer> counts) {
		counts.merge(n, 1, Integer::sum);
		if (counts.get(n) > 1) return;  // dedup recursion via identity
		switch (n) {
			case Atom l -> { /* no children */ }
			case KronProduct kp -> { walk(kp.outer, counts); walk(kp.inner, counts); }
			case KronChain kc -> { for (Node f : kc.factors) walk(f, counts); }
			case ConcatCols cr -> { walk(cr.left, counts); walk(cr.right, counts); }
			case ConcatRows cb -> { walk(cb.top, counts); walk(cb.bottom, counts); }
			case SumInner si -> { walk(si.left, counts); walk(si.right, counts); }
			case PeeledViaTa t -> { walk(t.cube, counts); walk(t.corner, counts); }
			case RecombinationN r -> {
				walk(r.base, counts);
				for (Node lf : r.leaves) walk(lf, counts);
			}
			case RecombinationTaN r -> {
				walk(r.base, counts);
				for (Node lf : r.leaves) walk(lf, counts);
			}
			case RecombinationWithPairN r -> {
				walk(r.base, counts);
				for (Node lf : r.leaves) walk(lf, counts);
			}
			case AugmentSquareDiscard a -> walk(a.square, counts);
			case Dce d -> walk(d.child, counts);
			case Transpose t -> walk(t.child, counts);
			case OrientAs o -> walk(o.child, counts);
			case AxisFlip af -> walk(af.child, counts);
			case AxisPermute ap -> walk(ap.child, counts);
			case DisjointSum ds -> { for (Node c : ds.children) walk(c, counts); }
			case SerendipitousProduct sp -> walk(sp.base, counts);
			case Project pr -> walk(pr.child, counts);
		}
	}

	/** Convenience: build a leaf node from a filename. The ref is the
	 *  canonical {@code {NxMxP}_m{R}_a{A}} key when the stem matches that
	 *  pattern (i.e. the file encodes the shape, rank, and additions in its
	 *  name) — so lineage references stay stable when the catalog re-tags
	 *  a file's source attribution. Falls back to the full stem when the
	 *  pattern doesn't match (e.g. {@code rosowski_2019_alg1-3x3x3_m21}
	 *  without a-count, or {@code Strassen<2,2,2>=7} short labels). */
	public static Atom atomFromFilename(String filename) {
		String ref = filename;
		if (ref.endsWith(".json")) ref = ref.substring(0, ref.length() - 5);
		int slash = ref.lastIndexOf('/');
		if (slash >= 0) ref = ref.substring(slash + 1);
		return new Atom(canonicalKey(ref));
	}

	private static final java.util.regex.Pattern CANONICAL_KEY =
			java.util.regex.Pattern.compile(".*_(\\d+x\\d+x\\d+_m\\d+_a\\d+)(?:_[A-Za-z0-9.]+)?$");

	/**
	 * Strip the source prefix from a filename stem. {@code fmm_lille-3x3x6_m40_a862}
	 * → {@code 3x3x6_m40_a862}. {@code alphaevolve-4x4x4_m48_a1264_0.5xC}
	 * → {@code 4x4x4_m48_a1264}. The remainder identifies the scheme up to
	 * field-tag equivalence; the source attribution lives separately in the
	 * catalog manifest so lineage references stay stable across re-attributions.
	 *
	 * <p>Returns the input unchanged when no canonical key can be extracted
	 * (short labels, formula names like {@code DIS09Lemma4(n=6)}, refs
	 * lacking the {@code _aN} suffix).</p>
	 */
	public static String canonicalKey(String filenameStem) {
		if (filenameStem == null) return null;
		java.util.regex.Matcher m = CANONICAL_KEY.matcher(filenameStem);
		return m.matches() ? m.group(1) : filenameStem;
	}

	/**
	 * Reduce a pinned ref {@code "{shape}@{hash}"} to the BARE {@code "{shape}"} for the
	 * human-readable {@code lineage_str}/{@code lineage_compact} display — the hash is an
	 * implementation detail for exact resolution and does not belong in the human form
	 * (user 2026-06-13). The machine-readable {@code lineage} JSON DAG ({@link #toJson})
	 * keeps the full ref so replay stays exact.
	 */
	private static String shortRefForDisplay(String ref) {
		int at = ref.indexOf('@');
		return at < 0 ? ref : ref.substring(0, at);
	}

	private static String jsonStr(String s) {
		StringBuilder sb = new StringBuilder("\"");
		for (int i = 0; i < s.length(); i++) {
			char c = s.charAt(i);
			if (c == '"' || c == '\\') sb.append('\\').append(c);
			else if (c == '\n') sb.append("\\n");
			else sb.append(c);
		}
		sb.append("\"");
		return sb.toString();
	}

	/** Convenience: KronChain folded from a list of leaves. */
	public static Node kronChain(Node... factors) {
		if (factors.length == 0) throw new IllegalArgumentException("empty chain");
		if (factors.length == 1) return factors[0];
		return new KronChain(List.of(factors));
	}

	/** Convenience: collect K leaves into a list. */
	public static List<Node> leaves(Node... ns) {
		List<Node> out = new ArrayList<>(ns.length);
		for (Node n : ns) out.add(n);
		return out;
	}

	/** Immediate child nodes of {@code n} (empty for a leaf {@link Atom} or a
	 *  parametric formula). Used by diagnostics to walk any node type uniformly. */
	public static List<Node> childrenOf(Node n) {
		return switch (n) {
			case Atom a -> List.of();
			case KronProduct kp -> List.of(kp.outer(), kp.inner());
			case KronChain kc -> List.copyOf(kc.factors());
			case ConcatCols c -> List.of(c.left(), c.right());
			case ConcatRows c -> List.of(c.top(), c.bottom());
			case SumInner si -> List.of(si.left(), si.right());
			case RecombinationN r -> prepend(r.base(), r.leaves());
			case RecombinationTaN r -> prepend(r.base(), r.leaves());
			case RecombinationWithPairN r -> prepend(r.base(), r.leaves());
			case AugmentSquareDiscard a -> List.of(a.square());
			case Dce d -> List.of(d.child());
			case Transpose t -> List.of(t.child());
			case OrientAs o -> List.of(o.child());
			case AxisFlip af -> List.of(af.child());
			case AxisPermute ap -> List.of(ap.child());
			case DisjointSum ds -> List.copyOf(ds.children());
			case SerendipitousProduct sp -> List.of(sp.base());
			case Project pr -> List.of(pr.child());
			case PeeledViaTa t -> List.of(t.cube(), t.corner());
		};
	}

	private static List<Node> prepend(Node head, List<Node> tail) {
		List<Node> out = new ArrayList<>(tail.size() + 1);
		out.add(head);
		out.addAll(tail);
		return out;
	}
}
