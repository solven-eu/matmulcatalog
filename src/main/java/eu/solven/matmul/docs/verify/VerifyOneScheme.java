package eu.solven.matmul.docs.verify;

import java.io.File;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.verifiers.Verifier;
import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.search.LineageReplayer;
import lombok.extern.slf4j.Slf4j;

/**
 * Benchmark / run the EXACT {@link Verifier#isExactNonCubic} verification of a
 * single scheme file. Handles both explicit-matrix files and lineage-only stubs
 * (replayed first via {@link LineageReplayer}). Prints the replay time and the
 * exact-verify time separately so you can see where the cost is.
 *
 * <pre>{@code
 * mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.verify.VerifyOneScheme \
 *     -Dexec.args="src/main/resources/schemes/derived/section18/derived_recursive-16x16x18_m2670_a123057.json"
 * }</pre>
 */
@Slf4j
public final class VerifyOneScheme {
	private VerifyOneScheme() {}

	public static void main(String[] args) {
		if (args.length < 1) {
			throw new IllegalArgumentException(
					"usage: VerifyOneScheme <scheme.json> [--field=Q] [--repeat=N]");
		}
		File file = new File(args[0]);
		String fieldTag = "Q";
		int repeat = 1;
		for (int i = 1; i < args.length; i++) {
			if (args[i].startsWith("--field=")) fieldTag = args[i].substring("--field=".length());
			else if (args[i].startsWith("--repeat=")) repeat = Integer.parseInt(args[i].substring("--repeat=".length()));
		}

		FieldAwareLookup committed =
				new FieldAwareLookup(Field.fromTag(fieldTag), java.nio.file.Path.of("src/main/resources/schemes"));
		LineageReplayer replayer = LineageReplayer.withDefaultPool(committed);

		long t0 = System.nanoTime();
		NonCubicBilinearAlgorithm alg = replayer.replayFromFile(file);
		long replayMs = (System.nanoTime() - t0) / 1_000_000L;
		log.info("replayed {} → ⟨{},{},{}⟩ r={}  ({} ms)",
				file.getName(), alg.n, alg.m, alg.p, alg.r, replayMs);

		for (int rep = 0; rep < repeat; rep++) {
			long s = System.nanoTime();
			// Size-aware: exact symbolic proof when the term-map fits, else the O(dim)
			// randomised spot-check — the Verdict says which tier ran (optimality
			// discipline: never report a spot-check as a proof). Unconditional
			// isExactNonCubic OOM'd on dense dim-28+ composites (22 GB term maps;
			// see Verifier.DEFAULT_MAX_EXACT_TERMS and the ⟨30,32,32⟩ precedent).
			Verifier.Verdict v = Verifier.verifyAuto(alg);
			long ms = (System.nanoTime() - s) / 1_000_000L;
			log.info("{}[{}/{}] = {}   (estTerms={}, {} ms){}",
					v.strategy(), rep + 1, repeat, v.ok(), v.estimatedTerms(), ms,
					v.isProof() ? "" : "  [spot-check, NOT an algebraic proof]");
		}
	}
}
