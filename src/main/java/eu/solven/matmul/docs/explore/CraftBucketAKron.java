package eu.solven.matmul.docs.explore;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.verifiers.Verifier;
import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.Compose;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.Lineage;
import eu.solven.matmul.catalog.MatrixJsonFormatter;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.search.LineageReplayer;
import tools.jackson.databind.JsonNode;
import tools.jackson.databind.json.JsonMapper;
import tools.jackson.databind.node.ArrayNode;
import tools.jackson.databind.node.ObjectNode;

/**
 * Craft the 11 "Bucket A" shapes whose published FMM-Lille rank is an exact
 * Kronecker product of two catalog atoms — schemes our {@code findBestStrategy}
 * Kron pass missed (it held both factors but never multiplied them at the right
 * orientation). Each is rebuilt deterministically via {@link Compose#kroneckerGeneral}
 * and written as a lineage-only {@code KronProduct} stub (no FMM dependency).
 *
 * <p>The recipes mirror the FMM-Lille per-format pages, e.g.
 * {@code ⟨12,12,15⟩=1280 = ⟨2,4,5⟩=32 ⊗ ⟨6,3,3⟩=40}. Orientations are explicit:
 * each factor is replayed at the exact orientation so the product lands on the
 * target shape.</p>
 *
 * <pre>mvn -q -o exec:java -Dexec.mainClass=eu.solven.matmul.docs.explore.CraftBucketAKron -Dexec.args="--apply"</pre>
 */
public final class CraftBucketAKron {
	private CraftBucketAKron() {}

	private static final JsonMapper MAPPER = JsonMapper.builder().build();

	/** target shape, expected rank, outer oriented ref, inner oriented ref. */
	private record Recipe(int n, int m, int p, int rank, String outer, String inner) {}

	private static final List<Recipe> RECIPES = List.of(
			new Recipe(4, 9, 9, 225, "2x3x3", "2x3x3"),
			new Recipe(4, 9, 12, 300, "2x3x3", "2x3x4"),
			new Recipe(4, 12, 15, 480, "2x3x3", "2x4x5"),
			new Recipe(6, 6, 16, 385, "2x2x2", "3x3x8"),
			new Recipe(6, 8, 14, 441, "2x2x2", "3x4x7"),
			new Recipe(6, 8, 15, 480, "2x4x5", "3x2x3"),
			new Recipe(6, 12, 12, 560, "2x2x2", "3x6x6"),
			new Recipe(9, 12, 12, 800, "3x2x4", "3x6x3"),
			new Recipe(12, 12, 12, 1040, "2x4x4", "6x3x3"),
			new Recipe(12, 12, 15, 1280, "2x4x5", "6x3x3"),
			new Recipe(14, 15, 16, 2016, "2x5x4", "7x3x4"));

	public static void main(String[] args) throws Exception {
		boolean apply = List.of(args).contains("--apply");
		FieldAwareLookup lookup = new FieldAwareLookup(Field.Q);
		LineageReplayer replayer = LineageReplayer.withDefaultPool(lookup);
		String root = "src/main/resources/schemes";

		int ok = 0, fail = 0;
		for (Recipe r : RECIPES) {
			try {
				NonCubicBilinearAlgorithm outer = replayer.replay(new Lineage.Atom(r.outer()));
				NonCubicBilinearAlgorithm inner = replayer.replay(new Lineage.Atom(r.inner()));
				NonCubicBilinearAlgorithm alg = Compose.kroneckerGeneral(outer, inner);
				boolean shapeOk = alg.n == r.n() && alg.m == r.m() && alg.p == r.p();
				boolean rankOk = alg.r == r.rank();
				boolean verifyOk = Verifier.passesRandomMatmulSpotCheck(alg);
				if (!(shapeOk && rankOk && verifyOk)) {
					System.out.printf("FAIL ⟨%d,%d,%d⟩: got ⟨%d,%d,%d⟩ r=%d shape=%s rank=%s verify=%s%n",
							r.n(), r.m(), r.p(), alg.n, alg.m, alg.p, alg.r, shapeOk, rankOk, verifyOk);
					fail++;
					continue;
				}
				Lineage.Node lineage = new Lineage.KronProduct(
						new Lineage.Atom(r.outer()), new Lineage.Atom(r.inner()));
				List<String> fields = intersectFields(lookup,
						sorted(r.outer()), sorted(r.inner()));
				int maxDim = Math.max(r.n(), Math.max(r.m(), r.p()));
				String fname = String.format("%dx%dx%d-r%d-derived_kron-%s.json",
						r.n(), r.m(), r.p(), r.rank(), SchemeIO.shortHash(alg));
				Path file = Path.of(root, "derived", "section" + maxDim, fname);
				System.out.printf("%-6s ⟨%d,%d,%d⟩=%d  = %s ⊗ %s  fields=%s  → %s%n",
						apply ? "CRAFT" : "(dry)", r.n(), r.m(), r.p(), r.rank(),
						r.outer(), r.inner(), fields, file.getFileName());
				if (apply) {
					Files.createDirectories(file.getParent());
					SchemeIO.writeStub(alg, file.toFile(), lineage);
					stampMetadata(file, fields);
				}
				ok++;
			} catch (RuntimeException e) {
				System.out.printf("FAIL ⟨%d,%d,%d⟩: %s%n", r.n(), r.m(), r.p(), e.toString());
				fail++;
			}
		}
		System.out.printf("%n%s: %d crafted, %d failed%n", apply ? "APPLIED" : "DRY RUN", ok, fail);
		if (!apply) System.out.println("Re-run with --apply to write.");
	}

	/** {@code "2x5x4"} → sorted {@code "2x4x5"} for catalog lookup of the factor. */
	private static String sorted(String ref) {
		String[] p = ref.split("x");
		int[] d = { Integer.parseInt(p[0]), Integer.parseInt(p[1]), Integer.parseInt(p[2]) };
		java.util.Arrays.sort(d);
		return d[0] + "x" + d[1] + "x" + d[2];
	}

	/** A Kron is valid over a field iff BOTH factors are — intersect their fields[]. */
	private static List<String> intersectFields(FieldAwareLookup lookup, String shA, String shB) {
		Set<String> a = readFields(lookup, shA);
		Set<String> b = readFields(lookup, shB);
		Set<String> out = new LinkedHashSet<>(a);
		out.retainAll(b);
		if (out.isEmpty()) out.addAll(List.of("Q", "R", "C"));  // safe char-0 fallback
		return new ArrayList<>(out);
	}

	private static Set<String> readFields(FieldAwareLookup lookup, String shape) {
		String[] p = shape.split("x");
		var path = lookup.findFile(Integer.parseInt(p[0]), Integer.parseInt(p[1]), Integer.parseInt(p[2]));
		Set<String> out = new LinkedHashSet<>();
		if (path.isPresent()) {
			try {
				out.addAll(SchemeIO.fieldTags(SchemeIO.parseJson(path.get().toFile())));
			} catch (Exception ignored) { /* fall back to caller default */ }
		}
		return out;
	}

	/** Add {@code fields[]} (array) + scalar {@code source/commutative/verified} to a stub. */
	private static void stampMetadata(Path file, List<String> fields) throws java.io.IOException {
		ObjectNode root = (ObjectNode) MAPPER.readTree(Files.readString(file));
		ArrayNode arr = root.arrayNode();
		fields.forEach(arr::add);
		root.set("fields", arr);
		root.put("source", "Derived_Kron");
		root.put("commutative", false);
		root.put("verified", true);
		JsonNode reordered = root;
		Files.writeString(file, MatrixJsonFormatter.format(reordered));
		FieldAwareLookup.onSchemeWritten(file.toFile());
	}
}
