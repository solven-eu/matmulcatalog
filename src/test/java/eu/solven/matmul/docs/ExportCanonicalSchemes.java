package eu.solven.matmul.docs;

import java.io.File;
import java.io.IOException;

import eu.solven.matmul.papers.laderman1976.Laderman23;
import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.papers.strassen1969.Strassen7;
import eu.solven.matmul.catalog.SchemeIO;

/**
 * One-shot exporter: writes our canonical hand-built algorithms to
 * {@code src/main/resources/schemes/} in the dronperminov JSON format so
 * tooling and external consumers can load them.
 *
 * <p>Run manually (it's a {@code main}, not a test) after any change to the
 * hand-built {@code Strassen7} / {@code Laderman23} factor matrices:
 * <pre>java … catalog.ExportCanonicalSchemes</pre>
 */
public class ExportCanonicalSchemes {

	public static void main(String[] args) throws IOException {
		File dir = new File("src/main/resources/schemes");
		dir.mkdirs();
		SchemeIO.write(NonCubicBilinearAlgorithm.fromCubic(Strassen7.get()),
				new File(dir, "strassen_222_r7.json"));
		SchemeIO.write(NonCubicBilinearAlgorithm.fromCubic(Laderman23.get()),
				new File(dir, "laderman_333_r23.json"));
		System.out.println("wrote schemes to " + dir.getAbsolutePath());
	}
}
