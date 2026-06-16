package eu.solven.matmul.catalog;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.FileReader;
import java.io.Reader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.junit.jupiter.api.Tag;
import org.junit.jupiter.api.Test;

import tools.jackson.databind.JsonNode;
import tools.jackson.databind.json.JsonMapper;

/**
 * Catalog closure guard (user 2026-06: "FMM-Lille atom while derived" — reported
 * repeatedly). A catalog entry must NOT be strictly beaten by a one-level
 * composition of other catalog entries: if {@code ⟨n,m,p⟩} can be obtained at
 * an equal-or-better rank by a single concat (direct sum along one axis) or a
 * single Kronecker factorisation of catalogued sub-shapes, then carrying a
 * worse scheme — and especially an opaque <em>imported</em> atom (FMM-Lille,
 * Perminov, …) — is a bug: we should materialise the cheaper DERIVATION (which
 * is replayable and atom-attributed) instead.
 *
 * <p>Rationale per the project rule "derive what we can, import what we can't":
 * an import should be the representative for a shape <em>only</em> when our own
 * derivation does not reach an equal-or-better rank. The canonical incident is
 * ⟨7,7,9⟩: an FMM-Lille atom at 315 while {@code ConcatCols(⟨7,7,2⟩=76,
 * ⟨7,7,7⟩=235)=311} is strictly cheaper.</p>
 *
 * <p><b>Scope.</b> The invariant is enforced up to {@link #GUARD_MAX_DIM} —
 * the band the closure pass keeps closed. Larger shapes (21–32) are an open
 * closure task (their spot-check expansion OOMs; tracked in PENDING_TASKS).
 * Genuinely-unavoidable residuals are listed in {@code
 * onelevel-closure-allowlist.txt} (canonical {@code n x m x p} per line; the
 * list may only shrink). Only char-0 (Z/Q/R/C), non-commutative ranks are
 * compared — concat/Kron of those compose to a valid non-commutative scheme;
 * F₂ and commutative ranks live in separate algebras and are excluded.</p>
 */
@Tag("catalog-iterating")
public class TestCatalogOneLevelClosure {

	/**
	 * Shapes with {@code maxDim ≤ this} must be closed under one-level composition.
	 * The dense / materialised band (≤16) plus the first stub band (17) are fully
	 * closed today; 18–32 is the open stub-closure frontier (the bounded recursive
	 * search doesn't always reach the full-catalog one-level optimum — a direct
	 * one-level-composition closer is the follow-up). Raise this as that frontier
	 * is closed; it should only ever go up.
	 */
	static final int GUARD_MAX_DIM = 17;

	@Test
	public void no_catalog_scheme_strictly_beaten_by_a_one_level_derivation() throws Exception {
		JsonMapper mapper = new JsonMapper();
		JsonNode root;
		try (Reader r = new FileReader(Path.of("docs/catalog.json").toFile())) {
			root = mapper.readTree(r);
		}
		JsonNode schemes = root.get("schemes");
		assertThat(schemes).as("docs/catalog.json must have schemes[]").isNotNull();

		// Best char-0, non-commutative rank per canonical (sorted) shape.
		Map<List<Integer>, Integer> best = new HashMap<>();
		for (JsonNode s : schemes) {
			JsonNode fmt = s.get("format");
			JsonNode rk = s.get("rank");
			if (fmt == null || !fmt.isArray() || fmt.size() != 3 || rk == null || !rk.isInt()) {
				continue;
			}
			if (s.path("commutative").asBoolean(false)) {
				continue;  // commutative ranks don't lift to NC composition
			}
			if (!hasChar0Field(s.get("fields"))) {
				continue;  // F₂-only ranks live in a separate algebra
			}
			List<Integer> key = sortedShape(fmt.get(0).asInt(), fmt.get(1).asInt(), fmt.get(2).asInt());
			int r = rk.asInt();
			best.merge(key, r, Math::min);
		}

		Set<List<Integer>> allow = loadAllowlist();

		List<String> beaten = new ArrayList<>();
		for (Map.Entry<List<Integer>, Integer> e : best.entrySet()) {
			List<Integer> shape = e.getKey();
			int n = shape.get(0), m = shape.get(1), p = shape.get(2);
			if (n == 1 || m == 1 || p == 1) {
				continue;  // width-1 axis: trivial, nothing to compose
			}
			if (Math.max(n, Math.max(m, p)) > GUARD_MAX_DIM) {
				continue;  // out of the guarded band (open closure task)
			}
			if (allow.contains(shape)) {
				continue;  // documented residual
			}
			Integer deriv = bestOneLevel(n, m, p, best);
			if (deriv != null && deriv < e.getValue()) {
				beaten.add(String.format("⟨%d,%d,%d⟩ catalog=%d but one-level derivation=%d",
						n, m, p, e.getValue(), deriv));
			}
		}

		assertThat(beaten)
				.as("catalog entries strictly beaten by a one-level concat/Kron derivation "
						+ "(materialise the cheaper derivation; see TestCatalogOneLevelClosure):\n  "
						+ String.join("\n  ", beaten))
				.isEmpty();
	}

	/** Cheapest one-level concat (sum) or Kronecker (product) bound, or null. */
	private static Integer bestOneLevel(int n, int m, int p, Map<List<Integer>, Integer> best) {
		Integer bound = null;
		// Concat along each axis: split one dim a+b, share the other two.
		bound = min(bound, axisSplit(n, m, p, 'n', best));
		bound = min(bound, axisSplit(n, m, p, 'm', best));
		bound = min(bound, axisSplit(n, m, p, 'p', best));
		// Kronecker: factor each axis into two.
		for (int n1 : divisors(n)) {
			for (int m1 : divisors(m)) {
				for (int p1 : divisors(p)) {
					int n2 = n / n1, m2 = m / m1, p2 = p / p1;
					if ((n1 == 1 && m1 == 1 && p1 == 1) || (n2 == 1 && m2 == 1 && p2 == 1)) {
						continue;
					}
					Integer r1 = best.get(sortedShape(n1, m1, p1));
					Integer r2 = best.get(sortedShape(n2, m2, p2));
					if (r1 != null && r2 != null) {
						bound = min(bound, r1 * r2);
					}
				}
			}
		}
		return bound;
	}

	private static Integer axisSplit(int n, int m, int p, char axis, Map<List<Integer>, Integer> best) {
		int dim = axis == 'n' ? n : axis == 'm' ? m : p;
		Integer bound = null;
		for (int a = 1; a < dim; a++) {
			int b = dim - a;
			Integer r1, r2;
			if (axis == 'n') { r1 = best.get(sortedShape(a, m, p)); r2 = best.get(sortedShape(b, m, p)); }
			else if (axis == 'm') { r1 = best.get(sortedShape(n, a, p)); r2 = best.get(sortedShape(n, b, p)); }
			else { r1 = best.get(sortedShape(n, m, a)); r2 = best.get(sortedShape(n, m, b)); }
			if (r1 != null && r2 != null) {
				bound = min(bound, r1 + r2);
			}
		}
		return bound;
	}

	private static Integer min(Integer a, Integer b) {
		if (a == null) return b;
		if (b == null) return a;
		return Math.min(a, b);
	}

	private static List<Integer> divisors(int x) {
		List<Integer> out = new ArrayList<>();
		for (int i = 1; i <= x; i++) if (x % i == 0) out.add(i);
		return out;
	}

	private static List<Integer> sortedShape(int a, int b, int c) {
		int[] s = { a, b, c };
		java.util.Arrays.sort(s);
		return List.of(s[0], s[1], s[2]);
	}

	private static boolean hasChar0Field(JsonNode fields) {
		if (fields == null || !fields.isArray()) {
			return true;  // unknown → treat as char-0 (don't silently drop)
		}
		for (JsonNode f : fields) {
			String t = f.asText();
			if (t.equals("Z") || t.equals("Q") || t.equals("R") || t.equals("C")) return true;
		}
		return false;
	}

	private static Set<List<Integer>> loadAllowlist() throws Exception {
		Set<List<Integer>> out = new HashSet<>();
		Path f = Path.of("src/test/resources/onelevel-closure-allowlist.txt");
		if (!Files.isReadable(f)) return out;
		for (String line : Files.readAllLines(f)) {
			String s = line.trim();
			if (s.isEmpty() || s.startsWith("#")) continue;
			s = s.split("#")[0].trim();
			var mm = java.util.regex.Pattern.compile("(\\d+)\\s*x\\s*(\\d+)\\s*x\\s*(\\d+)").matcher(s);
			if (mm.find()) {
				out.add(sortedShape(Integer.parseInt(mm.group(1)), Integer.parseInt(mm.group(2)),
						Integer.parseInt(mm.group(3))));
			}
		}
		return out;
	}
}
