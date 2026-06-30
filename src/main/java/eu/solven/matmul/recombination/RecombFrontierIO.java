package eu.solven.matmul.recombination;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

import tools.jackson.databind.json.JsonMapper;
import tools.jackson.databind.node.ArrayNode;
import tools.jackson.databind.node.ObjectNode;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.catalog.MatrixJsonFormatter;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.isotropy.RecombinationMultisetOrbit.Result;
import eu.solven.matmul.isotropy.RecombinationMultisetOrbit;

/**
 * Persist the GL-orbit recombination <b>frontier</b> of a base scheme as a sidecar JSON, keyed by
 * the scheme's content hash. One canonical scheme + a transform per frontier multiset (the
 * "1 scheme + transform config" model) — the recombination search reads this instead of
 * re-enumerating orbits or pre-expanding axis-flip pool variants every sweep.
 *
 * <p>Each entry carries the human-readable multiset, its block-index array, and the integer GL
 * transform {@code (X,Y,Z)} that {@link RecombinationMultisetOrbit#materialise} replays to rebuild
 * the concrete orbited scheme on demand (only for a chosen winner). The top-level {@code exhaustive}
 * flag records whether the frontier is certified complete (GL-exact + direction-bound stable) or a
 * partial menu (structural d≥4, or capped) — so a consumer never mistakes a bound for the truth.
 */
public final class RecombFrontierIO {
	private RecombFrontierIO() {}

	/** Default location: {@code src/main/resources/frontiers/{n}x{m}x{p}-{hash7}.json}. */
	public static Path defaultPath(NonCubicBilinearAlgorithm base) {
		String hash7 = SchemeIO.contentHash(base).substring(0, 7);
		return Path.of("src", "main", "resources", "frontiers",
				base.n + "x" + base.m + "x" + base.p + "-" + hash7 + ".json");
	}

	/**
	 * Write the frontier sidecar for {@code base}. {@code exhaustive} should be the certified
	 * completeness (e.g. {@link RecombinationMultisetOrbit#isStable} for the GL path; {@code false}
	 * for structural / capped). {@code method} is a free-text provenance ("GL-exact", "structural").
	 */
	public static Path write(NonCubicBilinearAlgorithm base, Result orbit, boolean exhaustive, String method, int dirBound)
			throws IOException {
		JsonMapper mapper = JsonMapper.builder().build();
		ObjectNode root = mapper.createObjectNode();
		root.put("shape", base.n + "x" + base.m + "x" + base.p);
		root.put("rank", base.r);
		root.put("hash", SchemeIO.contentHash(base));
		root.put("method", method);
		root.put("dirBound", dirBound);
		root.put("exhaustive", exhaustive);
		// CRITICAL honesty: "exhaustive" certifies completeness over THIS base's GL ORBIT only.
		// It is shape-complete iff the shape is known single-orbit — only ⟨2,2,2⟩ (de Groote 1978).
		// For every multi-orbit shape (⟨2,2,3⟩+) another base/orbit may add a non-dominated multiset;
		// the shape frontier is the cross-base-dominated UNION over known bases, never certified total.
		root.put("exhaustiveScope", isKnownSingleOrbit(base.n, base.m, base.p) ? "shape" : "orbit");

		List<String> frontier = orbit.dominanceFrontier();
		root.put("frontierSize", frontier.size());
		root.put("canonicalCount", orbit.canonicalMultisets.size());
		ArrayNode arr = root.putArray("frontier");
		for (String key : frontier) {
			ObjectNode e = mapper.createObjectNode();
			e.put("multiset", RecombinationMultisetOrbit.prettySymbolic(key, "n", "m", "p"));
			int[][] shapes = orbit.representativeShapes.get(key);
			if (shapes != null) {
				ArrayNode bi = e.putArray("blockIndices");
				for (int[] s : shapes) { ArrayNode row = bi.addArray(); for (int v : s) row.add(v); }
			}
			int[][][] xyz = orbit.representativeTransforms.get(key);
			if (xyz != null) {
				ObjectNode t = e.putObject("transform");
				putMatrix(t.putArray("X"), xyz[0]);
				putMatrix(t.putArray("Y"), xyz[1]);
				putMatrix(t.putArray("Z"), xyz[2]);
			} else {
				e.put("transform", (String) null); // structural path: no exact transform captured
			}
			arr.add(e);
		}

		Path out = defaultPath(base);
		Files.createDirectories(out.getParent());
		Files.writeString(out, MatrixJsonFormatter.format(mapper.writeValueAsString(root)));
		return out;
	}

	/** True iff the shape's rank-r decompositions are a SINGLE GL orbit — the only certified case is
	 *  ⟨2,2,2⟩ (de Groote 1978); every other shape is multi-orbit / unknown, so frontier coverage
	 *  there is per-orbit, not per-shape. */
	public static boolean isKnownSingleOrbit(int n, int m, int p) {
		int[] s = { n, m, p };
		java.util.Arrays.sort(s);
		return s[0] == 2 && s[1] == 2 && s[2] == 2;
	}

	/** A loaded sidecar: base dims {@code (bn,bm,bp)} and the frontier block-index arrays {@code [r][3]}. */
	public record Loaded(int[] dims, java.util.List<int[][]> frontier) {}

	/** Parse a frontier sidecar JSON (the {@code blockIndices} per frontier multiset). */
	public static Loaded read(Path file) throws IOException {
		JsonMapper mapper = JsonMapper.builder().build();
		var root = mapper.readTree(Files.readString(file));
		String[] xyz = root.get("shape").asString().split("x");
		int[] dims = { Integer.parseInt(xyz[0]), Integer.parseInt(xyz[1]), Integer.parseInt(xyz[2]) };
		java.util.List<int[][]> frontier = new java.util.ArrayList<>();
		for (var e : root.get("frontier")) {
			var bi = e.get("blockIndices");
			if (bi == null || bi.isNull()) continue;
			int[][] ms = new int[bi.size()][3];
			for (int k = 0; k < bi.size(); k++)
				for (int a = 0; a < 3; a++) ms[k][a] = bi.get(k).get(a).asInt();
			frontier.add(ms);
		}
		return new Loaded(dims, frontier);
	}

	private static void putMatrix(ArrayNode dst, int[][] M) {
		for (int[] row : M) { ArrayNode r = dst.addArray(); for (int v : row) r.add(v); }
	}
}
