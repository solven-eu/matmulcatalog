package eu.solven.matmul.docs.explore;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Optional;
import java.util.Set;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.verifiers.Verifier;
import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.Lineage;
import eu.solven.matmul.catalog.MatrixJsonFormatter;
import eu.solven.matmul.recombination.Recombination;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.search.CitedBound;
import eu.solven.matmul.search.LineageReplayer;
import tools.jackson.databind.JsonNode;
import tools.jackson.databind.json.JsonMapper;
import tools.jackson.databind.node.ArrayNode;
import tools.jackson.databind.node.ObjectNode;

/**
 * Reproduce each genuine {@code fmm_lille} file as a <strong>block
 * recombination</strong> of a catalog outer base. FMM's per-format pages show
 * these are recombinations (e.g. ⟨8,11,11⟩=641 over base ⟨2,3,3⟩ at
 * n=[4,4], m=[4,4,3], p=[4,4,3]). For each file:
 *
 * <ul>
 *   <li>try every outer base in {@link #BASES} × every <em>balanced</em> per-axis
 *       allocation (parts = ⌊d/k⌋/⌈d/k⌉, all distinct orderings — matching FMM's
 *       near-equal block sizes while staying tiny);</li>
 *   <li>predict the rank cheaply via {@link Recombination#recombineWithAllocation};
 *       keep candidates with predicted rank ≤ the FMM rank (<strong>tie-or-better</strong>);</li>
 *   <li>materialise the lowest-rank candidate via {@link Recombination#constructWithAllocation},
 *       verify, write a {@code RecombinationN} stub in {@code derived/}, drop the fmm.</li>
 * </ul>
 *
 * <p>Verify-gated → 0 regressions; converts at the achieved rank (so a strict
 * improvement is taken). {@code --band=N} selects maxDim; {@code --apply} writes.</p>
 */
public final class CraftFmmRecombination {
	private CraftFmmRecombination() {}

	private static final JsonMapper MAPPER = JsonMapper.builder().build();
	private static final String ROOT = "src/main/resources/schemes";

	/** Outer bases observed across the FMM-Lille recombination pages. */
	private static final String[] BASES = {
			"2x2x2", "2x2x3", "2x3x3", "2x4x4", "3x3x3", "2x5x5", "2x3x4", "2x4x5", "3x4x4" };

	public static void main(String[] args) throws Exception {
		boolean apply = List.of(args).contains("--apply");
		int band = List.of(args).stream().filter(a -> a.startsWith("--band="))
				.map(a -> Integer.parseInt(a.substring(7))).findFirst().orElse(8);
		// --shape=NxMxP restricts to a single fmm_lille file (exact orientation),
		// so re-running one shape never reprocesses — or races — its band siblings.
		String onlyShape = List.of(args).stream().filter(a -> a.startsWith("--shape="))
				.map(a -> a.substring(8)).findFirst().orElse(null);
		FieldAwareLookup lookup = new FieldAwareLookup(Field.Q);
		LineageReplayer replayer = LineageReplayer.withDefaultPool(lookup);
		CitedBound sota = new CitedBound(lookup);
		Recombination.AlgorithmLookup recombLookup = (a, b, c) -> {
			try { return Optional.of(replayer.replay(new Lineage.Atom(a + "x" + b + "x" + c))); }
			catch (RuntimeException e) { return Optional.empty(); }
		};
		// resolve each base AND its distinct axis-permutations once (skip any not in
		// the catalog). Rank is S₃-invariant, but the recombination needs the base in
		// the axis order that matches the target's per-axis part-counts — e.g.
		// ⟨9,10,11⟩=651 is a ⟨3,3,2⟩ base (p split into 2), not ⟨2,3,3⟩.
		List<int[]> baseShapes = new ArrayList<>();
		List<NonCubicBilinearAlgorithm> baseAlgs = new ArrayList<>();
		Set<String> seenBase = new LinkedHashSet<>();
		for (String b : BASES) {
			String[] tok = b.split("x");
			int[] dims = { Integer.parseInt(tok[0]), Integer.parseInt(tok[1]), Integer.parseInt(tok[2]) };
			for (int[] pm : distinctPerms(dims)) {
				String ps = pm[0] + "x" + pm[1] + "x" + pm[2];
				if (!seenBase.add(ps)) continue;
				try {
					NonCubicBilinearAlgorithm a = replayer.replay(new Lineage.Atom(ps));
					if (a.n == pm[0] && a.m == pm[1] && a.p == pm[2]) {  // correct orientation
						baseShapes.add(new int[] { a.n, a.m, a.p }); baseAlgs.add(a);
					}
				} catch (RuntimeException ignored) { /* base/orientation absent */ }
			}
		}

		List<Path> files;
		try (var s = Files.walk(Path.of(ROOT))) {
			files = s.filter(p -> p.getFileName().toString().matches(".*fmm_lille.*\\.json")).sorted().toList();
		}
		int converted = 0, kept = 0, improved = 0;
		for (Path p : files) {
			JsonNode root = SchemeIO.parseJson(p.toFile());
			JsonNode nArr = root.get("n");
			if (nArr == null || nArr.size() != 3) continue;
			int n = nArr.get(0).asInt(), m = nArr.get(1).asInt(), pp = nArr.get(2).asInt(), r = root.get("m").asInt();
			if (onlyShape != null) {
				if (!(n + "x" + m + "x" + pp).equals(onlyShape)) continue;
			} else if (Math.max(n, Math.max(m, pp)) != band) {
				continue;
			}

			Best best = findBest(baseShapes, baseAlgs, sota, recombLookup, n, m, pp, r);
			if (best == null) { kept++; continue; }
			Lineage.Node lineage = new Lineage.RecombinationN(new Lineage.Atom(best.base),
					best.aA, best.aB, best.aC, List.of());
			NonCubicBilinearAlgorithm alg = replayer.replay(lineage);
			String tag = best.rank < r ? "IMPROVE" : "tie";
			if (best.rank < r) improved++;
			System.out.printf("%-6s ⟨%d,%d,%d⟩  FMM=%d → %d (%s)  base=%s n=%s m=%s p=%s%n",
					apply ? "CRAFT" : "(dry)", n, m, pp, r, best.rank, tag, best.base,
					java.util.Arrays.toString(best.aA), java.util.Arrays.toString(best.aB),
					java.util.Arrays.toString(best.aC));
			converted++;
			if (apply) {
				int maxDim = Math.max(n, Math.max(m, pp));
				String fname = String.format("%dx%dx%d-r%d-derived_recomb-%s.json", n, m, pp, best.rank,
						SchemeIO.shortHash(alg));
				Path dst = Path.of(ROOT, "derived", "section" + maxDim, fname);
				Files.createDirectories(dst.getParent());
				SchemeIO.writeStub(alg, dst.toFile(), lineage);
				stamp(dst, fieldsOf(root));
				Files.deleteIfExists(p);
			}
		}
		System.out.printf("%nband %d %s: %d converted (%d strict-improvements), %d kept (no recombination tie)%n",
				band, apply ? "APPLIED" : "DRY RUN", converted, improved, kept);
	}

	private record Best(String base, int[] aA, int[] aB, int[] aC, int rank) {}

	/** Lowest-rank verifying recombination ≤ the FMM target across all bases/allocations. */
	private static Best findBest(List<int[]> baseShapes, List<NonCubicBilinearAlgorithm> baseAlgs,
			CitedBound sota, Recombination.AlgorithmLookup lookup, int n, int m, int p, int target) {
		record Cand(NonCubicBilinearAlgorithm base, String ref, int[] aA, int[] aB, int[] aC, int pred) {}
		List<Cand> cands = new ArrayList<>();
		for (int bi = 0; bi < baseShapes.size(); bi++) {
			int[] bs = baseShapes.get(bi); NonCubicBilinearAlgorithm base = baseAlgs.get(bi);
			String ref = bs[0] + "x" + bs[1] + "x" + bs[2];
			for (int[] aA : compositions(n, bs[0]))
				for (int[] aB : compositions(m, bs[1]))
					for (int[] aC : compositions(p, bs[2])) {
						Recombination.Result pr = Recombination.recombineWithAllocation(base, sota, aA, aB, aC);
						if (pr != null && pr.totalRank > 0 && pr.totalRank <= target) {
							cands.add(new Cand(base, ref, aA, aB, aC, (int) pr.totalRank));
						}
					}
		}
		cands.sort(java.util.Comparator.comparingInt(c -> c.pred));
		for (Cand c : cands) {
			try {
				NonCubicBilinearAlgorithm alg = Recombination.constructWithAllocation(c.base, lookup, c.aA, c.aB, c.aC);
				if (alg.n == n && alg.m == m && alg.p == p && alg.r == c.pred
						&& Verifier.passesRandomMatmulSpotCheck(alg)) {
					return new Best(c.ref, c.aA, c.aB, c.aC, c.pred);
				}
			} catch (RuntimeException ignored) { /* not realisable — next candidate */ }
		}
		return null;
	}

	/** Largest sub-block dimension to consider per part. Beyond this the sub-block is
	 *  unlikely to be a good building block and the rank explodes (predict-≤-target
	 *  prunes it anyway, but the cap keeps the enumeration bounded). */
	private static final int MAX_PART = 8;

	/** Every ordered composition of {@code d} into {@code k} positive parts, each in
	 *  [1, {@link #MAX_PART}]. Unbalanced splits included — FMM recipes routinely use
	 *  them (e.g. ⟨5,9,16⟩ splits p=16 as [6,6,4], not the balanced [6,5,5]). Empty
	 *  if {@code d < k}. */
	private static List<int[]> compositions(int d, int k) {
		List<int[]> out = new ArrayList<>();
		if (k <= 0 || d < k) return out;
		comp(d, 0, k, new int[k], out);
		return out;
	}

	private static void comp(int remaining, int idx, int k, int[] cur, List<int[]> out) {
		int partsLeft = k - idx;
		if (partsLeft == 1) {
			if (remaining >= 1 && remaining <= MAX_PART) { cur[idx] = remaining; out.add(cur.clone()); }
			return;
		}
		int max = Math.min(MAX_PART, remaining - (partsLeft - 1));  // leave ≥1 for each remaining part
		for (int v = 1; v <= max; v++) { cur[idx] = v; comp(remaining - v, idx + 1, k, cur, out); }
	}

	/** Distinct axis-permutations of a 3-tuple (1, 3 or 6 of them). */
	private static List<int[]> distinctPerms(int[] d) {
		int[][] all = {
				{ d[0], d[1], d[2] }, { d[0], d[2], d[1] }, { d[1], d[0], d[2] },
				{ d[1], d[2], d[0] }, { d[2], d[0], d[1] }, { d[2], d[1], d[0] } };
		Set<String> seen = new LinkedHashSet<>();
		List<int[]> out = new ArrayList<>();
		for (int[] p : all) if (seen.add(java.util.Arrays.toString(p))) out.add(p);
		return out;
	}

	private static void stamp(Path file, List<String> fields) throws java.io.IOException {
		ObjectNode root = (ObjectNode) MAPPER.readTree(Files.readString(file));
		ArrayNode arr = root.arrayNode();
		fields.forEach(arr::add);
		root.set("fields", arr);
		root.put("source", "Derived_Recombination");
		root.put("commutative", false);
		root.put("verified", true);
		Files.writeString(file, MatrixJsonFormatter.format(root));
		FieldAwareLookup.onSchemeWritten(file.toFile());
	}

	private static List<String> fieldsOf(JsonNode root) {
		Set<String> out = new LinkedHashSet<>(SchemeIO.fieldTags(root));
		if (out.isEmpty()) out.addAll(List.of("Q", "R", "C"));
		return new ArrayList<>(out);
	}
}
