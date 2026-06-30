package eu.solven.matmul.docs.explore;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
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
import tools.jackson.databind.json.JsonMapper;
import tools.jackson.databind.node.ArrayNode;
import tools.jackson.databind.node.ObjectNode;

/**
 * Craft the concat residue our {@code findBestStrategy} concat pass misses: a
 * target ⟨n,m,p⟩ that is two equal blocks of a smaller scheme glued along an axis
 * (rank additive). Currently the one shape the materialize sweep left regressed,
 * ⟨3,6,14⟩=188 = ⟨3,6,7⟩=94 {@code +p} ⟨3,6,7⟩=94. Written as a {@code ConcatCols}
 * stub. Same shape as the Kron crafter (CraftBucketAKron) but additive.
 */
public final class CraftConcatResidue {
	private CraftConcatResidue() {}

	private static final JsonMapper MAPPER = JsonMapper.builder().build();

	/** target shape, expected rank, the two p-concat halves (each an oriented ref). */
	private record Recipe(int n, int m, int p, int rank, String left, String right) {}

	private static final List<Recipe> RECIPES = List.of(
			new Recipe(3, 6, 14, 188, "3x6x7", "3x6x7"));

	public static void main(String[] args) throws Exception {
		boolean apply = List.of(args).contains("--apply");
		FieldAwareLookup lookup = new FieldAwareLookup(Field.Q);
		LineageReplayer replayer = LineageReplayer.withDefaultPool(lookup);
		String root = "src/main/resources/schemes";
		int ok = 0, fail = 0;
		for (Recipe r : RECIPES) {
			NonCubicBilinearAlgorithm a = replayer.replay(new Lineage.Atom(r.left()));
			NonCubicBilinearAlgorithm b = replayer.replay(new Lineage.Atom(r.right()));
			NonCubicBilinearAlgorithm alg = Compose.concatRight(a, b);  // glue along p
			boolean okShape = alg.n == r.n() && alg.m == r.m() && alg.p == r.p() && alg.r == r.rank();
			boolean verifyOk = okShape && Verifier.passesRandomMatmulSpotCheck(alg);
			if (!verifyOk) {
				System.out.printf("FAIL ⟨%d,%d,%d⟩: got ⟨%d,%d,%d⟩ r=%d verify=%s%n",
						r.n(), r.m(), r.p(), alg.n, alg.m, alg.p, alg.r, verifyOk);
				fail++;
				continue;
			}
			Lineage.Node lineage = new Lineage.ConcatCols(
					new Lineage.Atom(r.left()), new Lineage.Atom(r.right()));
			List<String> fields = readFields(lookup, r.left());
			int maxDim = Math.max(r.n(), Math.max(r.m(), r.p()));
			String fname = String.format("%dx%dx%d-r%d-derived_concat-%s.json",
					r.n(), r.m(), r.p(), r.rank(), SchemeIO.shortHash(alg));
			Path file = Path.of(root, "derived", "section" + maxDim, fname);
			System.out.printf("%-6s ⟨%d,%d,%d⟩=%d = %s +p %s  fields=%s → %s%n",
					apply ? "CRAFT" : "(dry)", r.n(), r.m(), r.p(), r.rank(),
					r.left(), r.right(), fields, file.getFileName());
			if (apply) {
				Files.createDirectories(file.getParent());
				SchemeIO.writeStub(alg, file.toFile(), lineage);
				ObjectNode node = (ObjectNode) MAPPER.readTree(Files.readString(file));
				ArrayNode arr = node.arrayNode();
				fields.forEach(arr::add);
				node.set("fields", arr);
				node.put("source", "Derived_Concat");
				node.put("commutative", false);
				node.put("verified", true);
				Files.writeString(file, MatrixJsonFormatter.format(node));
				FieldAwareLookup.onSchemeWritten(file.toFile());
			}
			ok++;
		}
		System.out.printf("%n%s: %d crafted, %d failed%n", apply ? "APPLIED" : "DRY RUN", ok, fail);
	}

	private static List<String> readFields(FieldAwareLookup lookup, String shape) {
		String[] p = shape.split("x");
		var path = lookup.findFile(Integer.parseInt(p[0]), Integer.parseInt(p[1]), Integer.parseInt(p[2]));
		Set<String> out = new LinkedHashSet<>();
		if (path.isPresent()) {
			try { out.addAll(SchemeIO.fieldTags(SchemeIO.parseJson(path.get().toFile()))); }
			catch (Exception ignored) { /* default below */ }
		}
		if (out.isEmpty()) out.addAll(List.of("Q", "R", "C"));
		return new java.util.ArrayList<>(out);
	}
}
