package eu.solven.matmul.catalog;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.junit.jupiter.api.Test;

import tools.jackson.databind.JsonNode;
import tools.jackson.databind.json.JsonMapper;

/**
 * Sanity check on the year metadata used by Pages' sort-by-year.
 *
 * <p>Mirrors the JS {@code classifySource} (in {@code docs/catalog.js}) to
 * detect entries whose source string yields no year — which makes them
 * collapse together when the user sorts by year. Catches:</p>
 *
 * <ul>
 *   <li>Cited bounds missing the {@code year} JSON field.</li>
 *   <li>Catalog schemes whose source doesn't match a known
 *       publication-year prefix.</li>
 *   <li>Derived bounds whose source string doesn't contain a 4-digit year.</li>
 * </ul>
 *
 * <p>The test prints a worklist for the human (PR author / catalog
 * curator) of every "missing year" entry. It <strong>fails only on
 * unexpected misses</strong> — known community catalogs without formal
 * publication (Perminov / fmm-lille / trivial axis-split bases /
 * "derived-*" entries) are on the allow-list and do not break the
 * build.</p>
 */
public class TestCatalogYearMetadata {

	/** Source prefixes for which a missing year is expected and acceptable. */
	private static final Pattern ALLOWED_NULL_YEAR = Pattern.compile(
			"^(Perminov|Dronperminov|Fmm[-_]lille|Derived|derived|Mul\\d{3}|"
			+ "mul\\d{3}|axis[-_]split|Hopcroft-Kerr ⟨|meta[-_]flip[-_]graph|fmm reduction|"
			// "unknown" = derived-recursive constructions whose provenance was never
			// stamped to a publication — genuinely year-less, like the derived family.
			+ "Solven[-_]|solven |unknown)",
			Pattern.CASE_INSENSITIVE);

	/** Known source → publication year. Sourced from REFERENCES.md. */
	private static final Map<String, Integer> KNOWN_YEARS = new LinkedHashMap<>();
	static {
		KNOWN_YEARS.put("Strassen", 1969);
		KNOWN_YEARS.put("Laderman", 1976);
		KNOWN_YEARS.put("Alphatensor-F2", 2022);
		KNOWN_YEARS.put("Alphatensor-Z", 2022);
		KNOWN_YEARS.put("AlphaTensor", 2022);  // bare label + alphatensor_Q schemes
		// AlphaTensor's published factorisation archive (alphatensor_Q imports).
		KNOWN_YEARS.put("https://github.com/google-deepmind/alphatensor", 2022);
		KNOWN_YEARS.put("Alphaevolve", 2025);
		KNOWN_YEARS.put("AlphaEvolve", 2025);
		KNOWN_YEARS.put("AlphaTensor 2022", 2022);
		KNOWN_YEARS.put("AlphaEvolve 2025", 2025);
		KNOWN_YEARS.put("Moosbauer", 2025);
		KNOWN_YEARS.put("Stapleton 2025", 2025);
		KNOWN_YEARS.put("Smirnov 2013", 2013);
		KNOWN_YEARS.put("Strassen 1969", 1969);
		KNOWN_YEARS.put("DIS09", 2009);  // Dumas-Iliopoulos-Saunders 2009 (cited bound tables)
	}

	private static final Pattern YEAR_IN_STRING = Pattern.compile("\\b(\\d{4})\\b");

	/**
	 * Mirrors {@code docs/catalog.js#classifySource} — returns the year
	 * extractable from the source string, or {@code null}.
	 */
	private static Integer classifyYear(String source) {
		if (source == null) return null;
		for (Map.Entry<String, Integer> e : KNOWN_YEARS.entrySet()) {
			if (source.regionMatches(true, 0, e.getKey(), 0, e.getKey().length())) {
				return e.getValue();
			}
		}
		String stripped = source
				.replaceFirst("(?i)^derived:\\s*", "")
				.replaceFirst("(?i)\\s*\\(via [^)]+\\)\\s*$", "")
				// Source labels keep '_' separators (e.g. Dumas_Pernet_Sedoglavic_2025);
				// underscore is a regex word-char so it hides the \b before a trailing
				// year. Normalise separators to spaces before the year scan.
				.replace('_', ' ')
				.trim();
		Matcher m = YEAR_IN_STRING.matcher(stripped);
		if (m.find()) {
			return Integer.parseInt(m.group(1));
		}
		return null;
	}

	@Test
	public void all_cited_bounds_have_year_field() throws IOException {
		JsonNode root = new JsonMapper().readTree(Path.of("docs/cited-bounds.json").toFile());
		List<String> missing = new ArrayList<>();
		for (JsonNode e : root.get("entries")) {
			if (!e.has("year") || e.get("year").isNull()) {
				missing.add(e.toString());
			}
		}
		assertThat(missing).as("cited-bounds entries missing 'year' field").isEmpty();
	}

	@Test
	public void catalog_sources_classify_to_year_unless_allowed() throws IOException {
		JsonNode root = new JsonMapper().readTree(Path.of("docs/catalog.json").toFile());
		List<String> unexpectedMissing = new ArrayList<>();
		List<String> allowedMissing = new ArrayList<>();
		java.util.Set<String> seenSources = new java.util.HashSet<>();
		for (JsonNode s : root.get("schemes")) {
			String src = s.has("source") ? s.get("source").asText() : null;
			if (src == null || !seenSources.add(src)) continue;
			Integer year = classifyYear(src);
			if (year != null) continue;
			if (ALLOWED_NULL_YEAR.matcher(src).find()) {
				allowedMissing.add(src);
			} else {
				unexpectedMissing.add(src);
			}
		}
		if (!allowedMissing.isEmpty()) {
			System.out.println("Catalog sources WITHOUT year (allow-listed; researchable):");
			allowedMissing.stream().sorted().forEach(s -> System.out.println("  • " + s));
		}
		assertThat(unexpectedMissing)
				.as("Catalog sources missing year that aren't in the allow-list — "
				  + "either add a publication year to KNOWN_YEARS, or whitelist via ALLOWED_NULL_YEAR")
				.isEmpty();
	}

	@Test
	public void derived_bounds_sources_contain_year_or_known_prefix() throws IOException {
		JsonNode root = new JsonMapper().readTree(Path.of("docs/derived-bounds.json").toFile());
		if (root.get("entries") == null) return;
		List<String> unexpectedMissing = new ArrayList<>();
		java.util.Set<String> seenSources = new java.util.HashSet<>();
		for (JsonNode e : root.get("entries")) {
			String src = e.has("source") ? e.get("source").asText() : null;
			if (src == null || !seenSources.add(src)) continue;
			Integer year = classifyYear(src);
			if (year != null) continue;
			if (ALLOWED_NULL_YEAR.matcher(src).find()) continue;
			unexpectedMissing.add(src);
		}
		assertThat(unexpectedMissing)
				.as("Derived-bound sources missing year that aren't in the allow-list")
				.isEmpty();
	}
}
