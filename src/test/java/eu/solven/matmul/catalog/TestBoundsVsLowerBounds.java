package eu.solven.matmul.catalog;

import static org.assertj.core.api.Assertions.assertThat;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.EnumSet;
import java.util.List;
import java.util.Optional;
import java.util.Set;

import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.KnownAlgorithm.Optimality;
import eu.solven.matmul.catalog.LowerBoundRegistry.Bound;
import org.junit.jupiter.api.Test;
import tools.jackson.databind.JsonNode;
import tools.jackson.databind.json.JsonMapper;

/**
 * Cross-checks every UPPER-bound claim we publish against every published
 * LOWER bound: a rank strictly below a proven floor is a contradiction, not a
 * record.
 *
 * <p><strong>Why this exists</strong> (regression guard, CLAUDE.md "every
 * high-level / silent bug gets a regression test"): {@link KnownAlgorithmCatalog}
 * carried {@code ⟨2,2,2⟩ commutative = 6} — attributed to Hopcroft–Kerr 1971 and
 * Winograd 1971, the two papers that prove 7 is minimal. It survived because
 * {@code docs/lower-bounds.json} (which recorded {@code lb = 7} for that exact
 * shape, from those exact papers) was display-only: the claim and its refutation
 * sat in one repo for months with nothing comparing them. The failure was silent
 * — no scheme file to verify, no crash, just a wrong number flowing into
 * {@code generated/COVERAGE.md} and the theory tables.</p>
 *
 * <p>Scope: the curated Java catalog, {@code docs/cited-bounds.json},
 * {@code docs/derived-from-cited-bounds.json} and the whole scheme manifest
 * {@code docs/catalog.json}.</p>
 */
public class TestBoundsVsLowerBounds {

	private static final JsonMapper MAPPER = JsonMapper.builder().build();

	private final LowerBoundRegistry registry = new LowerBoundRegistry();

	/** One violation, rendered so the failure message names the offending row. */
	private record Violation(String origin, String shape, String algebra, int rank, Bound floor) {
		@Override
		public String toString() {
			return String.format("%s: %s/%s claims rank %d, below %s",
					origin, shape, algebra, rank, floor.label());
		}
	}

	private void check(List<Violation> out, String origin, int[] format, Set<Field> fields,
			boolean commutative, int rank) {
		if (format == null || format.length != 3) return;
		for (Field f : fields) {
			Optional<Bound> floor = registry.binding(format[0], format[1], format[2], f, commutative);
			if (floor.isPresent() && rank < floor.get().lb()) {
				out.add(new Violation(origin,
						String.format("⟨%d,%d,%d⟩", format[0], format[1], format[2]),
						f.tag() + (commutative ? " (cmt)" : " (NC)"), rank, floor.get()));
			}
		}
	}

	@Test
	public void curated_catalog_never_claims_below_a_published_floor() {
		List<Violation> violations = new ArrayList<>();
		for (KnownAlgorithm a : KnownAlgorithmCatalog.all()) {
			check(violations, "KnownAlgorithmCatalog", new int[] { a.n, a.m, a.p },
					EnumSet.of(a.algebra.field()), a.algebra.commutative(), a.rank);
		}
		assertThat(violations).isEmpty();
	}

	/**
	 * The tier must match the evidence in both directions: never claim an
	 * optimum without a matching floor, and never leave a matched rank tagged
	 * as a mere bound.
	 */
	@Test
	public void curated_catalog_optimality_tiers_match_the_lower_bounds() {
		List<String> wrong = new ArrayList<>();
		for (KnownAlgorithm a : KnownAlgorithmCatalog.all()) {
			Optional<Bound> floor =
					registry.binding(a.n, a.m, a.p, a.algebra.field(), a.algebra.commutative());
			boolean matched = floor.isPresent() && floor.get().lb() == a.rank;
			if (matched && a.optimality != Optimality.PROVEN_OPTIMAL) {
				wrong.add(a + " — floor " + floor.get().label() + " matches, tag it PROVEN_OPTIMAL");
			}
			if (!matched && a.optimality == Optimality.PROVEN_OPTIMAL) {
				wrong.add(a + " — claims PROVEN_OPTIMAL but no matching floor ("
						+ floor.map(Bound::label).orElse("none published") + ")");
			}
		}
		assertThat(wrong).isEmpty();
	}

	@Test
	public void cited_bounds_never_claim_below_a_published_floor() {
		assertThat(scanBoundsFile(Path.of("docs/cited-bounds.json"), "cited-bounds")).isEmpty();
	}

	@Test
	public void derived_bounds_never_claim_below_a_published_floor() {
		assertThat(scanBoundsFile(Path.of("docs/derived-from-cited-bounds.json"), "derived-bounds"))
				.isEmpty();
	}

	@Test
	public void catalogued_schemes_never_claim_below_a_published_floor() {
		List<Violation> violations = new ArrayList<>();
		JsonNode root = readJson(Path.of("docs/catalog.json"));
		for (JsonNode s : root.get("schemes")) {
			Set<Field> fields = EnumSet.noneOf(Field.class);
			for (JsonNode f : s.get("fields")) fields.add(Field.fromTag(f.asString()));
			check(violations, "catalog.json", shapeOf(s), fields,
					s.path("commutative").asBoolean(false), s.get("rank").asInt());
		}
		assertThat(violations).isEmpty();
	}

	private List<Violation> scanBoundsFile(Path path, String origin) {
		List<Violation> violations = new ArrayList<>();
		for (JsonNode e : readJson(path).get("entries")) {
			if (!e.hasNonNull("format") || !e.hasNonNull("rank")) continue;
			Set<Field> fields = EnumSet.noneOf(Field.class);
			// "R/Q/Z" is one claim asserted over each of R, Q and Z — it must clear
			// the floor published for every field it names, not just the widest.
			for (String tag : e.get("field").asString().split("/")) fields.add(Field.fromTag(tag));
			check(violations, origin, shapeOf(e), fields,
					e.path("commutative").asBoolean(false), e.get("rank").asInt());
		}
		return violations;
	}

	private static int[] shapeOf(JsonNode node) {
		JsonNode fmt = node.get("format");
		if (fmt == null || !fmt.isArray() || fmt.size() != 3) return null;
		return new int[] { fmt.get(0).asInt(), fmt.get(1).asInt(), fmt.get(2).asInt() };
	}

	private static JsonNode readJson(Path path) {
		try {
			return MAPPER.readTree(Files.readString(path));
		} catch (java.io.IOException e) {
			throw new IllegalStateException("reading " + path, e);
		}
	}
}
