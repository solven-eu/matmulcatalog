package eu.solven.matmul.docs.migrate;

import java.io.File;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import tools.jackson.databind.JsonNode;
import tools.jackson.databind.ObjectMapper;
import tools.jackson.databind.node.ObjectNode;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.verifiers.Verifier;
import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.Lineage;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.search.LineageReplayer;

/**
 * Migrate lineages to UNAMBIGUOUS orientation refs:
 * <ol>
 *   <li>every {@code OrientAs} without an {@code axisMap} is augmented with the explicit axis-role
 *       permutation it currently resolves to ({@code child.orientPermFor(n,m,p)});</li>
 *   <li>every {@code Atom("oriented@hash")} whose hash is stored canonically at a DIFFERENT axis
 *       orientation is rewritten to {@code OrientAs(oriented, axisMap, Atom("canonical@hash"))} — so
 *       the Atom ref is a PRECISE {@code shape@hash} (no manifest false-dangling) and the orientation
 *       is explicit (no {@code orientAs} inference that could silently change when dims repeat).</li>
 * </ol>
 * Bit-exact by construction ({@code orientByPerm(perm) == orientAs(n,m,p)} for the captured perm).
 * Re-verifies each file replays to the same rank and passes {@link Verifier#isExactNonCubic} (unless
 * {@code --no-verify}). Args: file paths (or {@code --no-verify} then paths).
 */
public final class AugmentOrientAsAxisMaps {
	private AugmentOrientAsAxisMaps() {}

	private static final Pattern REF = Pattern.compile("^(\\d+)x(\\d+)x(\\d+)@([0-9a-f]{4,})$");
	private static final ObjectMapper MAPPER = new ObjectMapper();

	public static void main(String[] args) throws Exception {
		boolean verify = true;
		java.util.List<String> paths = new java.util.ArrayList<>();
		for (String a : args) { if (a.equals("--no-verify")) verify = false; else paths.add(a); }

		LineageReplayer rep = LineageReplayer.withDefaultPool(new FieldAwareLookup(Field.R));
		int changed = 0, ok = 0, fail = 0;
		for (String p : paths) {
			File f = new File(p);
			try {
				long before = verify ? rep.replayFromFile(f).r : -1;
				ObjectNode root = (ObjectNode) SchemeIO.parseJson(f);
				int[] n = { 0 };
				transform(root.get("lineage"), rep, n);
				if (n[0] == 0) continue;
				changed++;

				// IN-PLACE mutation only: the lineage JSON keeps its original @id/@ref dedup structure
				// (a readLineage→toJson round-trip would FLATTEN shared subtrees into broken
				// Atom("@ref?:L0") placeholders). Every other field (u/v/w, fields, …) is preserved.
				// Display strings are regenerated from the (dedup-intact) in-place JSON so they reflect
				// the new axisMaps; that parse resolves @ref correctly, so no "@ref?:L0" leaks.
				Lineage.Node parsed = SchemeIO.readLineage(root).orElseThrow();
				root.put("lineage_str", Lineage.prettyString(parsed));
				if (root.has("lineage_compact")) root.put("lineage_compact", Lineage.prettyCompact(parsed));
				java.nio.file.Files.writeString(f.toPath(),
						eu.solven.matmul.catalog.MatrixJsonFormatter.format(root));

				if (verify) {
					LineageReplayer rep2 = LineageReplayer.withDefaultPool(new FieldAwareLookup(Field.R));
					NonCubicBilinearAlgorithm after = rep2.replayFromFile(f);
					boolean exact = after.r == before && Verifier.isExactNonCubic(after);
					System.out.printf("%s %s: %d rewrites, r=%d (was %d) verified=%s%n",
							exact ? "✔" : "✗", f.getName(), n[0], after.r, before, exact);
					if (exact) ok++; else fail++;
				} else {
					System.out.printf("• %s: %d rewrites%n", f.getName(), n[0]);
					ok++;
				}
			} catch (Throwable t) {
				System.out.printf("✗ %s: EXC %s%n", f.getName(), t);
				fail++;
			}
		}
		System.out.printf("%n[done] changed=%d ok=%d fail=%d of %d%n", changed, ok, fail, paths.size());
	}

	private static void transform(JsonNode node, LineageReplayer rep, int[] count) {
		if (node == null) return;
		if (node.isArray()) { for (JsonNode c : node) transform(c, rep, count); return; }
		if (!node.isObject()) return;
		ObjectNode obj = (ObjectNode) node;
		String op = obj.path("op").asString("");

		if (op.equals("Atom") && obj.has("ref")) {
			Matcher m = REF.matcher(obj.get("ref").asString());
			if (m.matches()) {
				int rn = Integer.parseInt(m.group(1)), rm = Integer.parseInt(m.group(2)), rp = Integer.parseInt(m.group(3));
				String hash = m.group(4);
				int[] fs = fileShapeForHash(hash);
				if (fs != null && !(fs[0] == rn && fs[1] == rm && fs[2] == rp)) {
					NonCubicBilinearAlgorithm canon = replaySub(atomJson(fs, hash), rep);
					int[] perm = canon.orientPermFor(rn, rm, rp).orElseThrow(
							() -> new IllegalStateException("no S3 orientation " + rn + "x" + rm + "x" + rp + " of " + canon.n + "x" + canon.m + "x" + canon.p));
					JsonNode keepId = obj.get("id"); // preserve dedup @id if this Atom was shared
					obj.removeAll();
					obj.put("op", "OrientAs");
					if (keepId != null) obj.set("id", keepId);
					obj.put("n", rn); obj.put("m", rm); obj.put("p", rp);
					obj.put("axisMap", Lineage.axisMapToStr(perm));
					obj.set("child", atomJson(fs, hash));
					count[0]++;
					return;
				}
			}
		} else if (op.equals("OrientAs") && !obj.hasNonNull("axisMap")) {
			int tn = obj.path("n").asInt(0), tm = obj.path("m").asInt(0), tp = obj.path("p").asInt(0);
			NonCubicBilinearAlgorithm child = replaySub(obj.get("child"), rep);
			int[] perm = child.orientPermFor(tn, tm, tp).orElse(null);
			if (perm != null) { obj.put("axisMap", Lineage.axisMapToStr(perm)); count[0]++; }
			// recurse into child (may contain more)
			transform(obj.get("child"), rep, count);
			return;
		}
		for (JsonNode c : obj.values()) transform(c, rep, count);
	}

	private static ObjectNode atomJson(int[] shape, String hash) {
		ObjectNode a = MAPPER.createObjectNode();
		a.put("op", "Atom");
		a.put("ref", shape[0] + "x" + shape[1] + "x" + shape[2] + "@" + hash);
		return a;
	}

	private static NonCubicBilinearAlgorithm replaySub(JsonNode childJson, LineageReplayer rep) {
		ObjectNode wrap = MAPPER.createObjectNode();
		wrap.set("lineage", childJson);
		return rep.replay(SchemeIO.readLineage(wrap).orElseThrow());
	}

	private static int[] fileShapeForHash(String hash) {
		java.nio.file.Path dir = java.nio.file.Path.of("src/main/resources/schemes");
		String h7 = hash.substring(0, Math.min(7, hash.length()));
		try (var s = java.nio.file.Files.walk(dir)) {
			var hit = s.filter(p -> p.getFileName().toString().endsWith(".json")
					&& p.getFileName().toString().contains(h7)).findFirst();
			if (hit.isPresent()) {
				Matcher m = Pattern.compile("^(\\d+)x(\\d+)x(\\d+)-").matcher(hit.get().getFileName().toString());
				if (m.find()) return new int[] { Integer.parseInt(m.group(1)), Integer.parseInt(m.group(2)), Integer.parseInt(m.group(3)) };
			}
		} catch (Exception e) { /* ignore */ }
		return null;
	}
}
