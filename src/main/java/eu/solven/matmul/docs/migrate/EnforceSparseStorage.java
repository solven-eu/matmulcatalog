package eu.solven.matmul.docs.migrate;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.catalog.SchemeIO;

/**
 * Enforce {@code SPARSE_DIM_THRESHOLD} on the committed catalog: any dense
 * bilinear scheme with {@code maxDim ≥ threshold} is rewritten into the sparse
 * factor format (and, as a side effect, its coefficients become exact rationals,
 * since the sparse writer goes through {@link SchemeIO#formatCoef}). Already-sparse
 * files, stubs, complex schemes, and small dense files are untouched.
 *
 * <p>Metadata-preserving (delegates to {@link SchemeIO#convertDenseToSparse}).
 * Dry run by default; {@code --apply} writes. On apply, every converted file is
 * re-read and its rank/shape checked; a {@code 1/SAMPLE} fraction also gets a full
 * random spot-check (a full verify of all 1700+ would be slow).</p>
 *
 * <pre>
 *   mvn -q -o exec:java -Dexec.mainClass=eu.solven.matmul.docs.migrate.EnforceSparseStorage -Dexec.args="--apply"
 * </pre>
 */
public final class EnforceSparseStorage {
	private EnforceSparseStorage() {}

	private static final int SPOTCHECK_SAMPLE = 25;  // full-verify 1 in N converted files

	public static void main(String[] args) throws Exception {
		boolean apply = List.of(args).contains("--apply");
		Path root = Path.of("src/main/resources/schemes");
		List<Path> files;
		try (var s = Files.walk(root)) {
			files = s.filter(p -> p.getFileName().toString().endsWith(".json")).sorted().toList();
		}
		long before = 0, after = 0;
		int converted = 0, badRank = 0, verified = 0, verifyFail = 0;
		for (Path p : files) {
			File f = p.toFile();
			long sz = Files.size(p);
			before += sz;
			boolean would;
			try {
				would = SchemeIO.convertDenseToSparse(f, apply);
			} catch (RuntimeException e) {
				System.out.printf("ERROR  %s: %s%n", f.getName(), e.getMessage());
				after += sz;
				continue;
			}
			if (!would) {
				after += sz;
				continue;
			}
			converted++;
			if (apply) {
				after += Files.size(p);
				NonCubicBilinearAlgorithm a = SchemeIO.readBilinear(f);
				int[] decl = declaredShape(f);
				if (decl == null || a.n != decl[0] || a.m != decl[1] || a.p != decl[2]) {
					badRank++;
					System.out.printf("BAD-SHAPE %s → ⟨%d,%d,%d⟩%n", f.getName(), a.n, a.m, a.p);
				} else if (converted % SPOTCHECK_SAMPLE == 0) {
					boolean ok = Verifier.passesRandomMatmulSpotCheck(a);
					if (ok) verified++; else { verifyFail++; System.out.printf("VERIFY-FAIL %s%n", f.getName()); }
				}
			} else {
				after += sz;  // dry: unchanged size estimate
			}
		}
		System.out.printf("%n%s: %d converted dense→sparse, %d bad-shape, %d sampled-verified (%d fail).%n",
				apply ? "APPLIED" : "DRY RUN", converted, badRank, verified, verifyFail);
		if (apply) {
			System.out.printf("disk %.1f MB → %.1f MB%n", before / 1048576.0, after / 1048576.0);
		} else {
			System.out.println("Re-run with --apply to write.");
		}
	}

	private static int[] declaredShape(File f) throws java.io.IOException {
		var n = SchemeIO.parseJson(f).get("n");
		return (n != null && n.isArray() && n.size() == 3)
				? new int[] { n.get(0).asInt(), n.get(1).asInt(), n.get(2).asInt() } : null;
	}
}
