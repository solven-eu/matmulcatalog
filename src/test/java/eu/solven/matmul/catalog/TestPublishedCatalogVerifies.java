package eu.solven.matmul.catalog;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.FileReader;
import java.io.Reader;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.Tag;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.verifiers.Verifier;
import eu.solven.matmul.catalog.SchemeIO;
import tools.jackson.databind.JsonNode;
import tools.jackson.databind.json.JsonMapper;

/**
 * Every scheme referenced from {@code docs/catalog.json} (post-shave)
 * MUST verify. If a scheme is in the manifest, Pages publishes it —
 * publishing a scheme that doesn't actually compute matmul is a
 * correctness bug. Shaved-out files (not referenced by the manifest)
 * are out of scope for this test.
 *
 * <p>The field-appropriate verifier is chosen per row:</p>
 * <ul>
 *   <li>{@code F2}: {@link Verifier#isExactNonCubicF2} — F₂ schemes
 *       use modular arithmetic and need the mod-2 verifier.</li>
 *   <li>Any other field: {@link Verifier#passesRandomMatmulSpotCheck}
 *       — fast O(samples·r·(nm+mp+np)) randomised matmul check.</li>
 * </ul>
 *
 * <p>Complex schemes are skipped (would need a complex-arithmetic
 * verifier; tracked separately in ROADMAP).</p>
 */
@Tag("catalog-iterating")
public class TestPublishedCatalogVerifies {

	@Test
	public void every_manifested_scheme_verifies() throws Exception {
		JsonMapper mapper = new JsonMapper();
		JsonNode root;
		try (Reader r = new FileReader(Path.of("docs/catalog.json").toFile())) {
			root = mapper.readTree(r);
		}
		JsonNode schemes = root.get("schemes");
		assertThat(schemes).as("docs/catalog.json must have schemes[]").isNotNull();

		List<String> failed = new ArrayList<>();
		int total = 0, passed = 0, skipped = 0;
		for (JsonNode s : schemes) {
			JsonNode fileNode = s.get("file");
			if (fileNode == null || fileNode.isNull()) {
				skipped++;
				continue;
			}
			String fileRel = fileNode.asText();
			Path file = Path.of("src/main/resources/schemes", fileRel);
			total++;
			// Skip non-bilinear schemes — they use a different verifier
			// (Verifier.isExactNonBilinear), tracked separately.
			JsonNode rawScheme = mapper.readTree(file.toFile());
			if (rawScheme.has("scheme_type")
					&& "non_bilinear".equals(rawScheme.get("scheme_type").asText())) {
				skipped++;
				continue;
			}

			NonCubicBilinearAlgorithm alg;
			try {
				alg = SchemeIO.readBilinear(file.toFile());
			} catch (Exception e) {
				if (s.has("complex") && s.get("complex").asBoolean(false)) {
					skipped++;
					continue;
				}
				failed.add(fileRel + " (read error: " + e.getMessage() + ")");
				continue;
			}
			// The catalog no longer carries a singular `field` (2026-06-04): the
			// `fields` membership array is authoritative. Pick the verifier from it:
			// any characteristic-0 membership (Z/Q/R) → real-arithmetic spot check
			// (works for integer schemes even when they also list F2); else F2-only
			// → exact F2 verifier; else C-only → skip (no real-arith complex verifier).
			java.util.List<String> fields = new java.util.ArrayList<>();
			if (s.has("fields")) s.get("fields").forEach(f -> fields.add(f.asText()));
			String field = String.join("/", fields);
			boolean hasChar0 = fields.contains("Z") || fields.contains("Q") || fields.contains("R");
			boolean matmulOk;
			if (hasChar0) {
				matmulOk = Verifier.passesRandomMatmulSpotCheck(alg);
			} else if (fields.contains("F2")) {
				matmulOk = Verifier.isExactNonCubicF2(alg);
			} else {
				skipped++;  // C-only (or untagged) — no real-arithmetic complex verifier yet
				continue;
			}
			// Field-coefficient compliance: verify the U/V/W entries actually live
			// in the claimed field. The claim is the NARROWEST field the JSON
			// declares (its `fields` membership) — authoritative, per the
			// 2026-06-04 single-source-of-truth refactor. We deliberately do NOT
			// use the filename heuristic here: it mis-tagged `_Q-…` files as Z and
			// produced ~5100 false compliance failures even though those files
			// correctly declare fields=[Q,R,C]. This still catches genuine
			// over-claims — a scheme declaring Z whose coefficients are rational
			// fails checkAllInField(…, Z).
			eu.solven.matmul.algebra.Field claimedField = null;
			for (String t : new String[] { "Z", "Q", "R", "C", "F3", "F2" }) {
				if (fields.contains(t)) {
					claimedField = eu.solven.matmul.algebra.Field.fromTag(t);
					break;
				}
			}
			if (claimedField == null) {
				// No declared fields (legacy/stub entry) — default to Z (the catalog
				// is integer-coefficient by default). Never parse the filename.
				claimedField = eu.solven.matmul.algebra.Field.Z;
			}
			java.util.List<eu.solven.matmul.algebra.FieldCompliance.Discrepancy> fcDiffs =
					eu.solven.matmul.algebra.FieldCompliance.checkAllInField(alg, claimedField, 1);
			boolean fieldOk = fcDiffs.isEmpty();

			if (matmulOk && fieldOk) {
				passed++;
			} else if (!matmulOk) {
				failed.add(fileRel + " (field=" + field + ", rank=" + alg.r + " — matmul FAIL)");
			} else {
				failed.add(fileRel + " (claimed " + claimedField
						+ " — field-compliance FAIL: " + fcDiffs.get(0) + ")");
			}
		}
		System.out.printf("Published catalog: %d total, %d PASS, %d SKIP, %d FAIL%n",
				total, passed, skipped, failed.size());
		if (!failed.isEmpty()) {
			System.out.println("Failing entries:");
			for (String f : failed.subList(0, Math.min(20, failed.size()))) {
				System.out.println("  • " + f);
			}
		}
		assertThat(failed)
				.as("Every entry in docs/catalog.json (the published catalog) must verify; "
				  + "broken schemes need to be quarantined via FieldAwareLookup.KNOWN_BROKEN_FILES "
				  + "or removed from the manifest")
				.isEmpty();
	}
}
