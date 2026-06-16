package eu.solven.matmul.search.als;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.junit.jupiter.api.Tag;
import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.catalog.SchemeIO;

/**
 * One-shot registrar: ensures every ternary {@code ⟨2,2,2⟩} recombination
 * multiset (the 6 found by {@link Ternary2x2x2Orbit}) has a concrete base scheme
 * in the catalog, so the recombination sweep can actually use it. Bases not yet
 * covered are constructed as exact integer schemes (a {@code GL₂(ℤ)³}
 * change-of-basis of Strassen) and written to {@code curated/section2/}.
 */
class RegisterNovel2x2x2Bases {

	private static final String KNOWN2 = "src/main/resources/schemes/known/section2/";
	private static final String CURATED2 = "src/main/resources/schemes/curated/section2/";

	private static int[][][] seedFrom(NonCubicBilinearAlgorithm alg) {
		double[][] U = alg.denseU(), V = alg.denseV(), W = alg.denseW();
		int[][][] seed = new int[7][3][4];
		for (int k = 0; k < 7; k++) {
			for (int a = 0; a < 4; a++) seed[k][0][a] = (int) Math.round(U[a][k]);
			for (int b = 0; b < 4; b++) seed[k][1][b] = (int) Math.round(V[b][k]);
			for (int c = 0; c < 4; c++) seed[k][2][c] = (int) Math.round(W[c][k]);
		}
		return seed;
	}

	private static NonCubicBilinearAlgorithm algFrom(int[][][] s) {
		double[][] U = new double[4][7], V = new double[4][7], W = new double[4][7];
		for (int k = 0; k < 7; k++) {
			for (int a = 0; a < 4; a++) U[a][k] = s[k][0][a];
			for (int b = 0; b < 4; b++) V[b][k] = s[k][1][b];
			for (int c = 0; c < 4; c++) W[c][k] = s[k][2][c];
		}
		return new NonCubicBilinearAlgorithm(2, 2, 2, U, V, W);
	}

	/** Canonical multiset keys already covered by an exact-integer catalog base. */
	private static Map<String, String> coveredKeys() throws Exception {
		Map<String, String> covered = new LinkedHashMap<>();
		for (String dir : new String[] { KNOWN2, CURATED2 }) {
			File d = new File(dir);
			File[] files = d.listFiles((x, n) -> n.endsWith(".json"));
			if (files == null) continue;
			for (File f : files) {
				NonCubicBilinearAlgorithm alg;
				try { alg = SchemeIO.readBilinear(f); } catch (Exception e) { continue; }
				if (alg.n != 2 || alg.m != 2 || alg.p != 2 || alg.r != 7) continue;
				if (!Verifier.isExactNonCubic(alg)) continue; // skip F2-only / commutative
				covered.putIfAbsent(Ternary2x2x2Orbit.canonicalMultisetKey(seedFrom(alg)), f.getName());
			}
		}
		return covered;
	}

	private static String fmtMat(int[] m) {
		return "[[" + m[0] + "," + m[1] + "],[" + m[2] + "," + m[3] + "]]";
	}

	@Test
	@Tag("slow")
	void register() throws Exception {
		NonCubicBilinearAlgorithm strassen = SchemeIO.readBilinear(new File(KNOWN2, "2x2x2-r7-strassen-db11bcc.json"));
		int[][][] strassenSeed = seedFrom(strassen);

		// The 6 ternary multisets (seed-independent).
		List<int[][][]> seeds = new ArrayList<>();
		File[] kf = new File(KNOWN2).listFiles((d, n) -> n.endsWith(".json"));
		java.util.Arrays.sort(kf, Comparator.comparing(File::getName));
		for (File f : kf) {
			NonCubicBilinearAlgorithm a;
			try { a = SchemeIO.readBilinear(f); } catch (Exception e) { continue; }
			if (a.n == 2 && a.m == 2 && a.p == 2 && a.r == 7 && Verifier.isExactNonCubic(a)) seeds.add(seedFrom(a));
		}
		Ternary2x2x2Orbit.Result orbit = Ternary2x2x2Orbit.sweep(seeds, -1, 1);
		Set<String> ternaryKeys = new LinkedHashSet<>(orbit.representatives.keySet());
		Map<String, String> covered = coveredKeys();

		System.out.printf("%n%d ternary multisets; %d already covered by a base.%n", ternaryKeys.size(), covered.size());

		int registered = 0;
		for (String key : ternaryKeys) {
			if (covered.containsKey(key)) {
				System.out.printf("  covered: %-58s  [%s]%n", Ternary2x2x2Orbit.pretty(key), covered.get(key));
				continue;
			}
			// Construct an exact integer base realising this multiset, from Strassen.
			Ternary2x2x2Orbit.Realization real = Ternary2x2x2Orbit.findRealization(strassenSeed, key, 1);
			assertThat(real).as("realization for " + key).isNotNull();
			NonCubicBilinearAlgorithm alg = algFrom(real.scheme);
			assertThat(Verifier.isExactNonCubic(alg)).as("constructed base is exact").isTrue();

			// Tag from the shape distribution: # of full big-cubes (⟨9,9,9⟩) and small-cubes (⟨8,8,8⟩).
			int cube9 = countShape(key, "9,9,9"), cube8 = countShape(key, "8,8,8");
			String tag = "solven_orbit_c9x" + cube9 + "_c8x" + cube8;
			String hash7 = SchemeIO.contentHash(alg).substring(0, 7);
			File out = new File(CURATED2, "2x2x2-r7-" + tag + "-" + hash7 + ".json");

			SchemeIO.write(alg, out);
			String symbolic = symbolic(key);
			StringBuilder note = new StringBuilder();
			note.append("Ternary ⟨2,2,2⟩=7 base realising a recombination multiset not otherwise in the catalog. ");
			note.append("At a 2-part split (n₁≥n₂ per axis) the sub-shape multiset is ").append(symbolic).append(". ");
			note.append("Concretely at (9,8)³ → ⟨17,17,17⟩: ").append(Ternary2x2x2Orbit.pretty(key)).append(". ");
			if (cube9 >= 3) {
				note.append("Of independent interest: the ").append(cube9)
						.append("×⟨n₁,n₁,n₁⟩ cubic group admits 3-way trilinear aggregation (Pan), a cost lever ")
						.append("Strassen's own multiset does not expose. ");
			}
			note.append("Found by exact GL₂(ℚ)³-orbit enumeration (Ternary2x2x2Orbit / RecombinationMultisetOrbit); ")
					.append("de Groote-equivalent to Strassen.");

			Map<String, Object> fields = new LinkedHashMap<>();
			fields.put("source", "solven-orbit-enumeration-2026");
			fields.put("lineage_str", "GL₂(ℤ)³ change-of-basis of Strassen: X=" + fmtMat(real.X)
					+ ", Y=" + fmtMat(real.Y) + ", Z=" + fmtMat(real.Z));
			fields.put("lineage_compact", "strassen-1969-7mult^GL");
			fields.put("discovery_note", note.toString());
			fields.put("commutative", false);
			fields.put("zt", true);
			fields.put("verified", true);
			SchemeIO.updateFields(out, fields, List.of(), true);
			// fields[] / fields_not[] must be JSON arrays — updateFields only handles
			// scalars, so set them with Jackson and re-emit via the shared formatter.
			injectFieldArrays(out, List.of("F2", "F3", "Z", "Q", "R", "C"), List.of());

			// Round-trip: re-read and confirm it is exact and keeps the target multiset.
			NonCubicBilinearAlgorithm back = SchemeIO.readBilinear(out);
			assertThat(Verifier.isExactNonCubic(back)).as("written base re-reads exact").isTrue();
			assertThat(Ternary2x2x2Orbit.canonicalMultisetKey(seedFrom(back))).isEqualTo(key);

			System.out.printf("  REGISTERED: %-50s  → %s%n", Ternary2x2x2Orbit.pretty(key), out.getName());
			registered++;
		}
		System.out.printf("%nRegistered %d new curated base(s).%n", registered);
	}

	/** Set {@code fields[]} / {@code fields_not[]} as real JSON arrays, re-emit canonically. */
	private static void injectFieldArrays(File f, List<String> fields, List<String> fieldsNot) throws Exception {
		tools.jackson.databind.json.JsonMapper mapper = tools.jackson.databind.json.JsonMapper.builder().build();
		tools.jackson.databind.node.ObjectNode root =
				(tools.jackson.databind.node.ObjectNode) mapper.readTree(java.nio.file.Files.readString(f.toPath()));
		tools.jackson.databind.node.ArrayNode fa = root.arrayNode();
		for (String s : fields) fa.add(s);
		root.set("fields", fa);
		tools.jackson.databind.node.ArrayNode fna = root.arrayNode();
		for (String s : fieldsNot) fna.add(s);
		root.set("fields_not", fna);
		java.nio.file.Files.writeString(f.toPath(),
				eu.solven.matmul.catalog.MatrixJsonFormatter.format(root));
	}

	private static int countShape(String key, String shape) {
		int c = 0;
		for (String part : key.split("\\|")) if (part.equals(shape)) c++;
		return c;
	}

	/** Render the (9,8)-key symbolically: 9 ↦ the larger part Nᵢ₁, 8 ↦ smaller Nᵢ₂, per axis n/m/p. */
	private static String symbolic(String key) {
		String[] axis = { "n", "m", "p" };
		Map<String, Integer> counts = new java.util.TreeMap<>();
		for (String part : key.split("\\|")) {
			String[] xyz = part.split(",");
			StringBuilder sb = new StringBuilder();
			for (int a = 0; a < 3; a++) {
				if (a > 0) sb.append(',');
				sb.append(axis[a]).append(xyz[a].equals("9") ? "₁" : "₂");
			}
			counts.merge(sb.toString(), 1, Integer::sum);
		}
		StringBuilder out = new StringBuilder();
		boolean first = true;
		for (var e : counts.entrySet()) {
			if (!first) out.append(" + ");
			out.append(e.getValue()).append("·⟨").append(e.getKey()).append('⟩');
			first = false;
		}
		return out.toString();
	}
}
