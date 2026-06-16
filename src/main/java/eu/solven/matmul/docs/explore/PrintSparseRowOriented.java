package eu.solven.matmul.docs.explore;

import java.nio.file.Files;
import java.nio.file.Path;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.catalog.SchemeIO;

/**
 * One-shot: read a legacy sparse-format scheme, emit the row-oriented
 * {@code u_sparse}/{@code v_sparse}/{@code w_sparse} format introduced for
 * task #174, write to {@code target/preview-row-oriented.json} for visual
 * inspection. NOT registered anywhere — manual exec only.
 *
 * <pre>
 * mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.explore.PrintSparseRowOriented \
 *   -Dexec.args=src/main/resources/schemes/known/section4/dumas_pernet_sedoglavic_2025-4x4x4_m48_a960.json
 * </pre>
 */
public final class PrintSparseRowOriented {

	public static void main(String[] args) throws Exception {
		Path in = eu.solven.matmul.catalog.SchemeResolver.byHint(args.length > 0 ? args[0]
				: "src/main/resources/schemes/known/section4/dumas_pernet_sedoglavic_2025-4x4x4_m48_a960.json").toPath();
		NonCubicBilinearAlgorithm alg = SchemeIO.read(Files.readString(in));
		String out = SchemeIO.toJsonSparse(alg);
		Path target = Path.of("target", "preview-row-oriented.json");
		Files.createDirectories(target.getParent());
		Files.writeString(target, out);
		System.out.println("wrote " + target.toAbsolutePath() + " (" + out.length() + " bytes)");
	}
}
