package eu.solven.matmul.docs.migrate;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.catalog.Compose;
import eu.solven.matmul.catalog.MatrixJsonFormatter;
import eu.solven.matmul.catalog.SchemeIO;

/**
 * One-shot: materialise ⟨8,8,8⟩=336 = Strassen ⟨2,2,2⟩=7 ⊗ DPS-2025 ⟨4,4,4⟩=48
 * (Kronecker product). This valid Q/R scheme was previously only a
 * derived-bound stub (#158), so {@code FieldAwareLookup} (which reads scheme
 * files) couldn't see it — and the ⟨17,17,17⟩=2930 Winograd[9,8]³ recombination
 * landed on 2937 (⟨8,8,8⟩ falling back to Strassen³=343). Materialising the
 * actual factor matrices restores the 2930 decomposition
 * (336 + 388 + 486 + 4·430 = 2930).
 *
 * <p>Inner uses DPS-2025 (rational, Q/R-valid) rather than AlphaEvolve-48
 * (complex), so the product lifts to R for non-commutative recursion.</p>
 */
public final class MaterialiseKron888_336 {

	public static void main(String[] args) throws Exception {
		NonCubicBilinearAlgorithm strassen = SchemeIO.read(
				eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json"));
		NonCubicBilinearAlgorithm dps48 = SchemeIO.read(
				eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section4/dumas_pernet_sedoglavic_2025-4x4x4_m48_a960.json"));

		NonCubicBilinearAlgorithm kron = Compose.kroneckerGeneral(strassen, dps48);
		if (kron.n != 8 || kron.m != 8 || kron.p != 8 || kron.r != 336) {
			throw new IllegalStateException("expected ⟨8,8,8⟩=336, got ⟨" + kron.n + ","
					+ kron.m + "," + kron.p + "⟩=" + kron.r);
		}
		if (!Verifier.passesRandomMatmulSpotCheck(kron, 20_000, 7L)) {
			throw new IllegalStateException("Kron ⟨8,8,8⟩=336 failed sampled verification");
		}
		int adds = Verifier.additionCount(kron);

		// Emit with the new sparse format + fields[] + lineage.
		String body = SchemeIO.toJsonSparse(kron);
		// Inject metadata fields after the opening brace.
		String meta = "{\n"
				+ "  \"fields\": [\"Q\", \"R\", \"C\"],\n"
				+ "  \"fields_not\": [\"F2\", \"F3\", \"Z\"],\n"
				+ "  \"commutative\": false,\n"
				+ "  \"source\": \"Composed_Strassen_DPS2025\",\n"
				+ "  \"lineage_compact\": \"Kron[strassen-2x2x2_m7; dumas_pernet_sedoglavic_2025-4x4x4_m48] = 8x8x8_m336\",";
		body = meta + body.substring(body.indexOf('{') + 1);
		String pretty = MatrixJsonFormatter.format(body);

		Path out = Path.of("src/main/resources/schemes/known/section8",
				"derived_strassen_dps2025-8x8x8_m336_a" + adds + ".json");
		Files.writeString(out, pretty);
		System.out.println("wrote " + out + "  (rank=336, adds=" + adds + ", verified)");
	}
}
