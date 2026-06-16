package eu.solven.matmul.docs.generate;

import java.io.File;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.additions.SchemeAdditiveComplexity;
import eu.solven.matmul.catalog.SchemeIO;

/**
 * CLI for the CSE-minimised additive complexity of a scheme (tasks #189/#190).
 * Prints the naive addition count (= {@code Verifier.additionCount}) next to the
 * heuristically-minimal linear-straight-line-program count, per side, and (with
 * {@code --program}) the full reconstructable SLP for each side.
 *
 * <p>Run: {@code mvn -q -o exec:java
 * -Dexec.mainClass=eu.solven.matmul.docs.generate.AdditiveComplexityReport
 * -Dexec.args="src/main/resources/schemes/known/section2/winograd_1971-2x2x2_m7_a24.json [--program]"}</p>
 */
public final class AdditiveComplexityReport {

	private AdditiveComplexityReport() {}

	public static void main(String[] args) throws Exception {
		if (args.length == 0) {
			System.err.println("usage: AdditiveComplexityReport <scheme.json> [<scheme.json> …] [--program]");
			return;
		}
		boolean showProgram = false;
		for (String a : args) {
			if (a.equals("--program")) {
				showProgram = true;
			}
		}
		for (String a : args) {
			if (a.startsWith("--")) {
				continue;
			}
			NonCubicBilinearAlgorithm alg = SchemeIO.read(new File(a));
			SchemeAdditiveComplexity.Result r = SchemeAdditiveComplexity.analyse(alg);
			System.out.printf("%n%s%n  ⟨%d,%d,%d⟩ rank=%d%n  %s%n  (Verifier.additionCount=%d, SLP reconstructs=%b)%n",
					new File(a).getName(), alg.n, alg.m, alg.p, alg.r, r,
					Verifier.additionCount(alg), r.reconstructs(alg));
			if (showProgram) {
				System.out.println("  -- A-side (left factors L_k) --");
				r.aSlp().render("a").forEach(line -> System.out.println("    " + line));
				System.out.println("  -- B-side (right factors R_k) --");
				r.bSlp().render("b").forEach(line -> System.out.println("    " + line));
				System.out.println("  -- output (C_i from products m_k) --");
				r.outSlp().render("m").forEach(line -> System.out.println("    " + line));
			}
		}
	}
}
