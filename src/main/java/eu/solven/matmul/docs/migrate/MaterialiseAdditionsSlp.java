package eu.solven.matmul.docs.migrate;

import java.io.File;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;

import tools.jackson.databind.JsonNode;
import tools.jackson.databind.json.JsonMapper;
import tools.jackson.databind.node.ArrayNode;
import tools.jackson.databind.node.ObjectNode;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.additions.SchemeAdditiveComplexity;
import eu.solven.matmul.additions.Slp;
import eu.solven.matmul.catalog.SchemeIO;

/**
 * Dedicated, band-by-band process (tasks #189/#190) that derives the
 * CSE-minimised additive complexity of each dense scheme and stamps the
 * <em>full, replayable SLP</em> plus {@code min_additions} into its JSON — so
 * the additive construction can be reconstructed, not just counted.
 *
 * <p>This is intentionally <b>separate</b> from manifest generation: the
 * {@link SchemeAdditiveComplexity#analyse} CSE is expensive, so it must not run
 * on every catalog rebuild. Run it per dimension band, e.g.
 * {@code --band=14} or {@code --bands=2-16}, after a closure, then regenerate
 * the manifest (which merely reads the precomputed {@code min_additions}).</p>
 *
 * <p>Only schemes that carry explicit factor matrices are processed (the SLP
 * needs them); lineage-only stubs above {@code MATERIALISE_MAX_DIM} are skipped
 * (a future replay-based pass can fill them). Every emitted SLP is
 * self-verified ({@link Slp#reconstructs}) before being written.</p>
 *
 * <p>JSON shape added to each scheme:</p>
 * <pre>
 *   "min_additions": 15,
 *   "slp": {
 *     "additions": 15, "scalar_mults": 2,
 *     "a":   {"inputs": 4, "ops": [["+",4,2,3],["-",5,0,4]], "forms": [0,5,...]},
 *     "b":   {...},
 *     "out": {...}
 *   }
 * </pre>
 * where each op is {@code [kind,target,x,y]} for {@code +/-} or
 * {@code [kind,target,x,coeff]} for {@code *}; ids {@code < inputs} are input
 * variables, the rest are SSA temporaries; {@code forms[f]} is the slot holding
 * form {@code f} ({@code -1} = zero).
 */
public final class MaterialiseAdditionsSlp {

	private MaterialiseAdditionsSlp() {}

	public static void main(String[] args) throws Exception {
		File root = new File("src/main/resources/schemes");
		int bandMin = 2, bandMax = 32;
		boolean skipExisting = false;
		for (String a : args) {
			if (a.startsWith("--schemes-root=")) {
				root = new File(a.substring("--schemes-root=".length()));
			} else if (a.startsWith("--band=")) {
				bandMin = bandMax = Integer.parseInt(a.substring("--band=".length()));
			} else if (a.startsWith("--bands=")) {
				String[] r = a.substring("--bands=".length()).split("-");
				bandMin = Integer.parseInt(r[0]);
				bandMax = Integer.parseInt(r[1]);
			} else if (a.equals("--skip-existing")) {
				skipExisting = true;
			}
		}
		JsonMapper mapper = JsonMapper.builder().build();
		List<File> files = new ArrayList<>();
		collect(root, files);
		// Ascending by max-dim so we proceed band by band (user 2026-06-04).
		files.sort(java.util.Comparator.comparingInt(MaterialiseAdditionsSlp::maxDimOf));

		long t0 = System.nanoTime();
		int processed = 0, stamped = 0, skipped = 0, failed = 0, verifyFail = 0, tooBig = 0;
		int curBand = -1;
		for (File f : files) {
			JsonNode node;
			try {
				node = mapper.readTree(f);
			} catch (Exception e) {
				continue;
			}
			if (node == null || !node.isObject() || !node.has("n")) {
				continue;
			}
			JsonNode nArr = node.get("n");
			if (!nArr.isArray() || nArr.size() != 3) {
				continue;
			}
			int maxDim = Math.max(nArr.get(0).asInt(),
					Math.max(nArr.get(1).asInt(), nArr.get(2).asInt()));
			if (maxDim < bandMin || maxDim > bandMax) {
				continue;
			}
			if (maxDim != curBand) {
				curBand = maxDim;
				System.out.printf("== band max-dim=%d ==%n", maxDim);
			}
			// Skip stubs (no matrices) and complex schemes (additions over C differ).
			if ("stub".equals(node.path("scheme_type").asText(""))
					|| SchemeIO.isComplex(node)) {
				skipped++;
				continue;
			}
			if (skipExisting && node.has("min_additions") && node.has("slp")) {
				skipped++;
				continue;
			}
			NonCubicBilinearAlgorithm alg;
			try {
				alg = SchemeIO.read(f);
			} catch (Exception e) {
				skipped++;
				continue;
			}
			// Cost guard: the greedy CSE is ~ r·nnz²·rounds. Skip pathologically
			// large dense schemes (e.g. ⟨16,16,16⟩=2304 dense) so a single bad
			// scheme can't stall a band; revisit when a faster minimiser lands.
			long nVars = Math.max((long) alg.n * alg.m,
					Math.max((long) alg.m * alg.p, (long) alg.n * alg.p));
			if ((long) alg.r * nVars > 400_000L) {
				tooBig++;
				continue;
			}
			processed++;
			try {
				SchemeAdditiveComplexity.Result r = SchemeAdditiveComplexity.analyse(alg);
				if (!r.reconstructs(alg)) {
					verifyFail++;
					System.err.printf("⟨%d,%d,%d⟩ %s — SLP did not reconstruct, NOT stamped%n",
							alg.n, alg.m, alg.p, f.getName());
					continue;
				}
				ObjectNode obj = (ObjectNode) node;
				obj.put("min_additions", r.minimal());
				// Honesty flag (CLAUDE.md optimality discipline): the value is a
				// greedy upper bound, NOT a proven global minimum → always false for
				// now. Flip to true only when produced by an exact solver.
				obj.put("additions_optimal", false);
				ObjectNode slp = mapper.createObjectNode();
				slp.put("additions", r.minimal());
				slp.put("scalar_mults", r.scalarMults());
				// Honesty tier (CLAUDE.md optimality discipline): this is a greedy,
				// cancellation-free heuristic → a valid UPPER BOUND, NOT a proven
				// minimal LSP (the problem is NP-hard). Labelled so no consumer
				// mistakes it for the optimum.
				slp.put("optimality", "bound");
				slp.put("optimality_method", "greedy cancellation-free CSE (Paar-style)");
				slp.set("a", sideJson(mapper, r.aSlp()));
				slp.set("b", sideJson(mapper, r.bSlp()));
				slp.set("out", sideJson(mapper, r.outSlp()));
				obj.set("slp", slp);
				// The single shared canonical-write entry point — NOT
				// obj.toPrettyString(), whose default spacing/inlining re-introduces
				// the non-canonical style `sanitize` exists to remove.
				eu.solven.matmul.catalog.MatrixJsonFormatter.write(f, obj);
				stamped++;
			} catch (RuntimeException e) {
				failed++;
				System.err.printf("⟨%d,%d,%d⟩ %s — analyse FAILED: %s%n",
						alg.n, alg.m, alg.p, f.getName(), e.getMessage());
			}
			if (processed % 200 == 0) {
				long ms = (System.nanoTime() - t0) / 1_000_000L;
				System.out.printf("[progress] %d processed, %d stamped, %d skipped, %d fail, %dms%n",
						processed, stamped, skipped, failed, ms);
			}
		}
		System.out.printf("%nbands %d-%d: %d processed, %d stamped, %d skipped, %d too-big, %d analyse-fail, %d verify-fail%n",
				bandMin, bandMax, processed, stamped, skipped, tooBig, failed, verifyFail);
	}

	/** Max axis parsed from a {@code …-NxMxP_…} filename (99 if not found, so
	 *  unparseable files sort last). Cheap key for ascending band order. */
	private static int maxDimOf(File f) {
		java.util.regex.Matcher m = java.util.regex.Pattern
				.compile("(\\d+)x(\\d+)x(\\d+)").matcher(f.getName());
		if (!m.find()) {
			return 99;
		}
		return Math.max(Integer.parseInt(m.group(1)),
				Math.max(Integer.parseInt(m.group(2)), Integer.parseInt(m.group(3))));
	}

	private static ObjectNode sideJson(JsonMapper mapper, Slp slp) {
		ObjectNode side = mapper.createObjectNode();
		side.put("inputs", slp.nInputs());
		ArrayNode ops = mapper.createArrayNode();
		for (Slp.Op o : slp.ops()) {
			ArrayNode op = mapper.createArrayNode();
			op.add(String.valueOf(o.kind()));
			op.add(o.target());
			op.add(o.x());
			if (o.kind() == '*') {
				op.add(o.c());
			} else {
				op.add(o.y());
			}
			ops.add(op);
		}
		side.set("ops", ops);
		ArrayNode forms = mapper.createArrayNode();
		for (int fr : slp.formResult()) {
			forms.add(fr);
		}
		side.set("forms", forms);
		return side;
	}

	private static void collect(File dir, List<File> out) {
		File[] kids = dir.listFiles();
		if (kids == null) {
			return;
		}
		for (File k : kids) {
			if (k.isDirectory()) {
				collect(k, out);
			} else if (k.getName().endsWith(".json")) {
				out.add(k);
			}
		}
	}
}
