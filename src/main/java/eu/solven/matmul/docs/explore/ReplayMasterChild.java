package eu.solven.matmul.docs.explore;

import java.io.File;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.search.LineageReplayer;

/**
 * Throwaway: replay a master scheme's lineage IN THE BRANCH catalog (leaves resolve to our
 * current best) and report the rank it actually reproduces + whether it verifies. If the
 * replayed rank &gt; master's claimed rank, master's number is a phantom relative to our
 * catalog (its leaves are gone / worse here) — same finding as ⟨6,18,26⟩ (claim 1726 →
 * replay 1732). Usage: {@code --file=path.json --claim=RANK}.
 */
public class ReplayMasterChild {

	public static void main(String[] args) {
		String file = null;
		int claim = -1;
		for (String a : args) {
			if (a.startsWith("--file=")) file = a.substring(7);
			else if (a.startsWith("--claim=")) claim = Integer.parseInt(a.substring(8));
		}
		if (file == null) {
			System.err.println("need --file=path.json [--claim=RANK]");
			return;
		}
		FieldAwareLookup lookup = new FieldAwareLookup("Q");
		LineageReplayer rep = LineageReplayer.withDefaultPool(lookup);
		long t0 = System.currentTimeMillis();
		try {
			NonCubicBilinearAlgorithm alg = rep.replayFromFile(new File(file));
			boolean ok = Verifier.isExactNonCubic(alg);
			String verdict = claim < 0 ? "?"
					: (alg.r <= claim ? "REPRODUCES (<= claim)" : "PHANTOM (replay " + alg.r + " > claim " + claim + ")");
			System.out.printf("⟨%d,%d,%d⟩ claim=%d  replayed=%d  verifies=%s  => %s  (%dms)%n",
					alg.n, alg.m, alg.p, claim, alg.r, ok, verdict, System.currentTimeMillis() - t0);
		} catch (Throwable e) {
			System.out.printf("%s replay FAILED: %s%n", file, e);
		}
	}
}
