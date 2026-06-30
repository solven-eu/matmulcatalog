package eu.solven.matmul.docs.migrate;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.verifiers.Verifier;
import eu.solven.matmul.catalog.Compose;
import eu.solven.matmul.catalog.MatrixJsonFormatter;
import eu.solven.matmul.catalog.SchemeIO;

/**
 * One-shot: materialise ⟨3,6,6⟩=80 as the trivial inner-dimension doubling of
 * Smirnov's ⟨3,3,6⟩=40 — i.e. ⟨3,6,6⟩ = ⟨3,3,6⟩ ⊕inner ⟨3,3,6⟩ (≡ ⊗⟨1,2,1⟩),
 * rank 40+40 = 80.
 *
 * <p>Background: ⟨3,6,6⟩=80 was previously carried as a raw
 * {@code fmm_lille-3x6x6_m80} import, which mis-credited the rank to the
 * FMM-Lille aggregator. The FMM digest itself cites no reference for it
 * (references: []) precisely because the rank is not a discovery — it is the
 * free split of the shared dimension 6 = 2×3 into two independent ⟨3,3,6⟩
 * products (visible in the original scheme as an odd/even column-parity split).
 * Per "derive what we can, import what we can't", we DERIVE it from the held
 * ⟨3,3,6⟩=40 (Smirnov 2013) and drop the import. The rank content stays
 * Smirnov's; the Source becomes our composition.</p>
 */
public final class MaterialiseConcat366_80 {

	public static void main(String[] args) throws Exception {
		File baseFile =
				eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section6/fmm_lille-3x3x6_m40_a862_b0.json");
		NonCubicBilinearAlgorithm base = SchemeIO.read(baseFile);
		if (base.n != 3 || base.m != 3 || base.p != 6 || base.r != 40) {
			throw new IllegalStateException("expected base ⟨3,3,6⟩=40, got ⟨" + base.n + ","
					+ base.m + "," + base.p + "⟩=" + base.r);
		}

		// Split the shared inner dimension: ⟨3,3,6⟩ ⊕inner ⟨3,3,6⟩ = ⟨3,6,6⟩, r=80.
		NonCubicBilinearAlgorithm concat = Compose.concatInner(base, base);
		if (concat.n != 3 || concat.m != 6 || concat.p != 6 || concat.r != 80) {
			throw new IllegalStateException("expected ⟨3,6,6⟩=80, got ⟨" + concat.n + ","
					+ concat.m + "," + concat.p + "⟩=" + concat.r);
		}
		if (!Verifier.passesRandomMatmulSpotCheck(concat, 20_000, 1e-9)) {
			throw new IllegalStateException("concat ⟨3,6,6⟩=80 failed sampled verification");
		}
		int adds = Verifier.additionCount(concat);

		// Reference the base by its catalog shape token (source-stripped), matching
		// the lineage-ref convention used by other derived schemes.
		String baseRef = "3x3x6-m40";
		String body = SchemeIO.toJsonSparse(concat);
		String meta = "{\n"
				+ "  \"fields\": [\"Q\", \"R\", \"C\"],\n"
				+ "  \"fields_not\": [\"F2\", \"F3\", \"Z\"],\n"
				+ "  \"commutative\": false,\n"
				+ "  \"source\": \"Derived_Concat\",\n"
				+ "  \"lineage_compact\": \"" + baseRef + " +m " + baseRef + " = 3x6x6_m80\",\n"
				+ "  \"lineage\": {\"op\": \"ConcatInner\","
				+ " \"left\": {\"op\": \"Atom\", \"ref\": \"" + baseRef + "\"},"
				+ " \"right\": {\"op\": \"Atom\", \"ref\": \"" + baseRef + "\"}},";
		body = meta + body.substring(body.indexOf('{') + 1);
		String pretty = MatrixJsonFormatter.format(body);

		Path out = Path.of("src/main/resources/schemes/known/section6",
				"derived_concat-3x6x6_m80_a" + adds + ".json");
		Files.writeString(out, pretty);
		System.out.println("wrote " + out + "  (rank=80, adds=" + adds + ", verified)");
	}
}
