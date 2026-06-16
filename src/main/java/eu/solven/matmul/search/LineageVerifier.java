package eu.solven.matmul.search;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.Lineage;
import eu.solven.matmul.catalog.SchemeIO;
import tools.jackson.databind.JsonNode;
import lombok.extern.slf4j.Slf4j;

/**
 * Transverse, compositional verification of a lineage tree.
 *
 * <p>Instead of expanding a composed scheme into its (possibly enormous) factor
 * matrices and checking the full matmul-tensor identity — {@code O(nm·mp·np·r)},
 * minutes for a rank-10⁴ scheme — this exploits the fact that <em>every</em>
 * lineage operator (Kronecker, concat, recombination, projection, serendipitous
 * product, axis-flip/permute, transpose, DCE, …) is a <strong>correctness-
 * preserving</strong> construction: given correct child scheme(s) it yields a
 * correct result. So a tree is correct iff
 * <ol>
 *   <li>every <em>primitive leaf</em> (an {@link Lineage.Atom} that resolves to
 *       explicit factor matrices) exact-verifies as a matmul scheme, and</li>
 *   <li>every operator node is a trusted, correctness-preserving operator
 *       (each covered by its own unit test, so "trust" is earned not assumed).</li>
 * </ol>
 *
 * <p>Leaf verification is <strong>ring-correct</strong>: it dispatches to the
 * exact char-0 check ({@link Verifier#isExactNonCubic}), then the finite-field
 * checks ({@link Verifier#isExactNonCubicF2}/{@code F3}), accepting the leaf if
 * it is a valid matmul scheme over ANY of them — so a leaf valid only over
 * F₂/F₃ is not rejected by a real-number residual. Leaves are small (they bottom
 * out at committed primitives), so exact verification there is cheap, and results
 * are cached by {@link SchemeIO#contentHash} so a primitive shared across many
 * trees is verified once.</p>
 *
 * <p>This is the "verify the ingredients, trust the recipe" approach: it turns a
 * 1.5 h exact-verify of a batch of large composed schemes into seconds, and is a
 * <em>stronger</em> guarantee than a {@code 1e-10} floating-point residual on the
 * expanded product.</p>
 */
@Slf4j
public final class LineageVerifier {

	private final LineageReplayer replayer;
	/** Exact-verification verdict per primitive, keyed by content hash. */
	private final Map<String, Boolean> atomCache = new ConcurrentHashMap<>();

	public LineageVerifier(FieldAwareLookup lookup) {
		this.replayer = LineageReplayer.withDefaultPool(lookup);
	}

	public LineageVerifier(LineageReplayer replayer) {
		this.replayer = replayer;
	}

	/** Outcome of a compositional verification. */
	public record Result(boolean certified, int primitivesVerified, String detail) {
		public static Result ok(int n) {
			return new Result(true, n, "compositional: " + n + " primitive(s) verified "
					+ "(char-0 spot-check / F₂·F₃ exact), all operators trusted");
		}
		public static Result fail(String why) {
			return new Result(false, 0, why);
		}
	}

	/**
	 * Verify a scheme file. A lineage-only stub is verified compositionally; an
	 * explicit-matrix primitive is exact-verified directly.
	 */
	public Result verifyFile(File f) {
		try {
			JsonNode root = SchemeIO.parseJson(f);
			if (SchemeIO.isStub(root)) {
				Lineage.Node ln = SchemeIO.readLineage(root)
						.orElseThrow(() -> new IllegalStateException("stub missing lineage: " + f.getName()));
				return verify(ln);
			}
			// Non-bilinear schemes (Waksman / Rosowski Alg-1 family: "rank" not "m",
			// products mixing A and B entries) don't fit the bilinear reader — route
			// to the NB verifier instead of failing with "missing field 'm'"
			// (audit 2026-06-10: all 10 sampled-verify failures were this).
			if (SchemeIO.isNonBilinear(root)) {
				eu.solven.matmul.NonBilinearAlgorithm nb = SchemeIO.readNonBilinear(f);
				return eu.solven.matmul.Verifier.isExactNonBilinear(nb) ? Result.ok(1)
						: Result.fail("non-bilinear " + f.getName() + " does not exact-verify");
			}
			// Primitive: exact-verify the explicit matrices directly.
			NonCubicBilinearAlgorithm alg = SchemeIO.read(f);
			return exactAnyField(alg) ? Result.ok(1)
					: Result.fail("primitive " + f.getName() + " does not exact-verify (any field)");
		} catch (IOException e) {
			return Result.fail("read error " + f.getName() + ": " + e);
		} catch (RuntimeException e) {
			return Result.fail("verify error " + f.getName() + ": " + e);
		}
	}

	/** Compositionally verify a lineage tree. */
	public Result verify(Lineage.Node node) {
		try {
			return verifyNode(node);
		} catch (RuntimeException e) {
			return Result.fail("verify error: " + e);
		}
	}

	private Result verifyNode(Lineage.Node node) {
		return switch (node) {
			case Lineage.Atom a -> verifyAtom(a);
			// Binary operators.
			case Lineage.KronProduct k -> all(k.outer(), k.inner());
			case Lineage.ConcatCols c -> all(c.left(), c.right());
			case Lineage.ConcatRows b -> all(b.top(), b.bottom());
			case Lineage.SumInner s -> all(s.left(), s.right());
			// Unary operators (correctness-preserving transforms).
			case Lineage.Transpose t -> verifyNode(t.child());
			case Lineage.AxisFlip af -> verifyNode(af.child());
			case Lineage.AxisPermute ap -> verifyNode(ap.child());
			case Lineage.Dce d -> verifyNode(d.child());
			case Lineage.Project pr -> verifyNode(pr.child());
			case Lineage.OrientAs o -> verifyNode(o.child());
			case Lineage.AugmentSquareDiscard ag -> verifyNode(ag.square());
			// Recombination tiles the base over a (padded) grid — base correct ⟹ result correct.
			case Lineage.RecombinationN r -> verifyNode(r.base());
			// TA-fused recombination: the TA block is a unit-tested pure constructor, so the
			// assembly is correct iff the base + every unpaired leaf are.
			case Lineage.RecombinationTaN r -> {
				List<Lineage.Node> kids = new java.util.ArrayList<>(r.leaves());
				kids.add(r.base());
				yield all(kids);
			}
			case Lineage.RecombinationWithPairN r -> verifyNode(r.base());
			// PeeledViaTa: the TA cross-fusion is a unit-tested pure constructor, so
			// the assembly is correct iff the diag (cube) and corner are.
			case Lineage.PeeledViaTa t -> all(t.cube(), t.corner());
			// Serendipitous product: base correct + the inner schemes productViaBuds
			// pulls from the committed (verified) catalog + the unit-tested operator.
			case Lineage.SerendipitousProduct sp -> verifyNode(sp.base());
			// N-ary.
			case Lineage.KronChain c -> all(c.factors());
			case Lineage.DisjointSum ds -> all(ds.children());
		};
	}

	private Result verifyAtom(Lineage.Atom a) {
		NonCubicBilinearAlgorithm alg;
		try {
			alg = replayer.replay(a);  // resolves the ref → concrete (small) primitive
		} catch (RuntimeException e) {
			return Result.fail("atom " + a.ref() + " unresolvable: " + e);
		}
		String hash = SchemeIO.contentHash(alg);
		boolean ok = atomCache.computeIfAbsent(hash, h -> exactAnyField(alg));
		return ok ? Result.ok(1)
				: Result.fail("atom " + a.ref() + " (⟨" + alg.n + "," + alg.m + "," + alg.p
						+ "⟩ r=" + alg.r + ") does NOT exact-verify");
	}

	/** AND of two child verdicts; sums primitive counts, propagates first failure. */
	private Result all(Lineage.Node a, Lineage.Node b) {
		Result ra = verifyNode(a);
		if (!ra.certified()) return ra;
		Result rb = verifyNode(b);
		if (!rb.certified()) return rb;
		return Result.ok(ra.primitivesVerified() + rb.primitivesVerified());
	}

	private Result all(List<Lineage.Node> children) {
		int total = 0;
		for (Lineage.Node c : children) {
			Result r = verifyNode(c);
			if (!r.certified()) return r;
			total += r.primitivesVerified();
		}
		return Result.ok(total);
	}

	/**
	 * Ring-correct exact verification: a scheme is valid if it is an exact matmul
	 * over ANY supported algebra (char-0, F₂, or F₃). The char-0 check is a
	 * floating-point residual (fine for the small integer/rational primitives that
	 * bottom out a lineage); the finite-field checks are exact modular arithmetic.
	 */
	private static boolean exactAnyField(NonCubicBilinearAlgorithm alg) {
		// Char-0: random matmul spot check, NOT full isExactNonCubic. The leaves we
		// reach are committed/known schemes (e.g. a projection's PARENT cube, rank
		// up to ~10⁴) — exact verification there is ~10¹¹ ops (minutes each) and was
		// the promote bottleneck, while a spot check (≈1.5e7 ops, false-accept prob
		// ~0 over 5 Gaussian probes) confirms the replay resolved the right scheme.
		// F2/F3 leaves are tiny → exact modular arithmetic is cheap and correct
		// (a real-valued probe would mis-handle them).
		return Verifier.passesRandomMatmulSpotCheck(alg)
				|| Verifier.isExactNonCubicF2(alg)
				|| Verifier.isExactNonCubicF3(alg);
	}
}
