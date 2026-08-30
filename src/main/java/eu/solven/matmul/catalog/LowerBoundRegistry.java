package eu.solven.matmul.catalog;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.EnumSet;
import java.util.List;
import java.util.Optional;
import java.util.Set;

import eu.solven.matmul.algebra.Field;
import tools.jackson.databind.JsonNode;
import tools.jackson.databind.json.JsonMapper;

/**
 * Machine-readable view of {@code docs/lower-bounds.json} — the published
 * <strong>impossibility</strong> results (rank floors), as opposed to the
 * upper bounds carried by schemes / {@code cited-bounds.json} /
 * {@link KnownAlgorithmCatalog}.
 *
 * <p>Until 2026-08 this file was display-only (the SPA rendered it; nothing
 * consulted it). That is how the catalog came to assert a commutative
 * {@code R_c(⟨2,2,2⟩) = 6} while this very file recorded
 * {@code lb = 7, source = "Winograd 1971 + Hopcroft–Kerr 1971"} — the claim
 * and its refutation, in the same repo, never compared. This class exists so
 * that an upper-bound claim below a published floor is a <em>test failure</em>
 * rather than a number nobody cross-reads; see {@code TestBoundsVsLowerBounds}.</p>
 *
 * <h2>Two applicability axes</h2>
 *
 * <p><strong>Field.</strong> A floor proven over {@code K} also binds every
 * <em>sub</em>-field of {@code K} of the same characteristic: a scheme over
 * {@code Z} is valid over {@code R}, so {@code R_R ≤ R_Z}, so
 * {@code R_Z ≥ R_R ≥ lb_R}. It does NOT transfer across characteristics in
 * either direction — the {@code field} token in the JSON names exactly what
 * was proven, and {@link #fieldsBoundBy} only widens it downward inside one
 * characteristic.</p>
 *
 * <p><strong>Model.</strong> {@code model: "bilinear"} (the default) is a floor
 * on the non-commutative rank {@code R} and binds non-commutative claims only —
 * a commutative-only algorithm may legitimately sit below it, since
 * {@code R_c ≤ R}. {@code model: "quadratic"} is a floor on the multiplicative
 * complexity, proven with commutativity of the inputs allowed, and therefore
 * binds <em>both</em> models. This axis is the one the ⟨2,2,2⟩ error hid on:
 * the field was right ({@code "all"}), the model was unrecorded.</p>
 */
public final class LowerBoundRegistry {

	/** Default location of the registry, relative to the repo root. */
	public static final Path DEFAULT_PATH = Path.of("docs/lower-bounds.json");

	/** Algebraic model a floor is proven in — see the class javadoc. */
	public enum Model {
		/** Floor on the non-commutative rank {@code R}. Binds NC claims only. */
		BILINEAR,
		/** Floor on the multiplicative complexity (commutativity allowed). Binds both models. */
		QUADRATIC;

		static Model fromTag(String tag) {
			if (tag == null) return BILINEAR;
			return switch (tag) {
				case "bilinear" -> BILINEAR;
				case "quadratic" -> QUADRATIC;
				default -> throw new IllegalArgumentException("unknown lower-bound model: " + tag);
			};
		}
	}

	/** One published floor. {@code format} is the canonical (sorted) shape. */
	public record Bound(int[] format, Set<Field> fields, Model model, int lb, String source, int year) {

		/** Whether this floor constrains a claim of the given {@code (shape, field, model)}. */
		public boolean binds(int n, int m, int p, Field field, boolean commutative) {
			int[] sorted = { n, m, p };
			java.util.Arrays.sort(sorted);
			if (!java.util.Arrays.equals(sorted, format)) return false;
			if (commutative && model != Model.QUADRATIC) return false;
			return fields.contains(field);
		}

		public String label() {
			return String.format("⟨%d,%d,%d⟩ ≥ %d (%s %d, %s)",
					format[0], format[1], format[2], lb, source, year,
					model == Model.QUADRATIC ? "quadratic/commutativity-allowed" : "bilinear/NC");
		}
	}

	private final List<Bound> bounds;

	public LowerBoundRegistry() {
		this(DEFAULT_PATH);
	}

	public LowerBoundRegistry(Path path) {
		this.bounds = load(path);
	}

	/** All floors, in file order. */
	public List<Bound> bounds() {
		return List.copyOf(bounds);
	}

	/**
	 * The strongest published floor constraining a claim over
	 * {@code (shape, field, commutative)}, or empty when none applies.
	 */
	public Optional<Bound> binding(int n, int m, int p, Field field, boolean commutative) {
		Bound best = null;
		for (Bound b : bounds) {
			if (!b.binds(n, m, p, field, commutative)) continue;
			if (best == null || b.lb() > best.lb()) best = b;
		}
		return Optional.ofNullable(best);
	}

	/** Convenience: the strongest applicable floor, or {@code 0} when none is published. */
	public int lowerBound(int n, int m, int p, Field field, boolean commutative) {
		return binding(n, m, p, field, commutative).map(Bound::lb).orElse(0);
	}

	/**
	 * Fields bound by a floor published over {@code token}. Widens downward
	 * within one characteristic only (see the class javadoc): {@code C} also
	 * binds {@code R, Q, Z}; {@code F2}/{@code F3} bind nothing else.
	 *
	 * <p>Tokens: {@code "all"}, {@code "char0"}, a single field, or a
	 * {@code '/'}-separated set (e.g. {@code "F2/Z/Q/R/C"}).</p>
	 */
	public static Set<Field> fieldsBoundBy(String token) {
		if (token == null) throw new IllegalArgumentException("null field token");
		Set<Field> base = EnumSet.noneOf(Field.class);
		switch (token) {
			case "all" -> base.addAll(EnumSet.allOf(Field.class));
			case "char0" -> base.addAll(EnumSet.of(Field.Z, Field.Q, Field.R, Field.C));
			default -> {
				for (String t : token.split("/")) base.add(Field.fromTag(t));
			}
		}
		Set<Field> out = EnumSet.noneOf(Field.class);
		for (Field f : base) {
			out.add(f);
			// A floor over K binds every sub-field of K: rank can only drop as the
			// field grows, so R_subfield ≥ R_K ≥ lb. fallbackChain() is exactly
			// "fields a scheme may come from to be valid over K" = the sub-fields.
			out.addAll(f.fallbackChain());
		}
		return out;
	}

	private static List<Bound> load(Path path) {
		List<Bound> out = new ArrayList<>();
		JsonNode root;
		try {
			root = JsonMapper.builder().build().readTree(Files.readString(path));
		} catch (IOException e) {
			throw new IllegalStateException("loading " + path, e);
		}
		for (JsonNode e : root.get("entries")) {
			// 'kind' entries record a different quantity (border rank R̃), not a
			// floor on R — skip rather than silently reinterpret them.
			if (e.hasNonNull("kind")) continue;
			JsonNode fmt = e.get("format");
			int[] shape = { fmt.get(0).asInt(), fmt.get(1).asInt(), fmt.get(2).asInt() };
			java.util.Arrays.sort(shape);
			out.add(new Bound(shape,
					fieldsBoundBy(e.get("field").asString()),
					Model.fromTag(e.hasNonNull("model") ? e.get("model").asString() : null),
					e.get("lb").asInt(),
					e.get("source").asString(),
					e.get("year").asInt()));
		}
		return out;
	}
}
