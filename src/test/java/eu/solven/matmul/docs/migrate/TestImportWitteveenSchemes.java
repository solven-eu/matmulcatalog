package eu.solven.matmul.docs.migrate;

import static org.assertj.core.api.Assertions.assertThat;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.docs.migrate.ImportWitteveenSchemes.Candidate;

/**
 * The scheduled pull from {@code MerlijnW70/fmm-schemes} must pick up exactly the
 * upstream files the catalog does not hold yet — by {@code (shape, rank, class)},
 * not by name — and nothing that is not a scheme.
 */
public class TestImportWitteveenSchemes {

	private static final List<String> LISTING = List.of(
			"README.md",
			"SHA256SUMS",
			".zenodo.json",
			"schemes/11x13x15_m1364_ZT.json",
			"schemes/7x11x15_m772_ZT.json",
			"schemes/2x13x15_m304_ZT.json",
			"schemes/notes.json",
			"scripts/verify.py");

	@Test
	public void plan_lists_only_schemes_not_held_at_their_class() {
		Set<String> have = new HashSet<>();
		KnownSchemeKeys.add(have, 11, 13, 15, 1364, "ZT"); // held (ternary)
		KnownSchemeKeys.add(have, 2, 13, 15, 304, "Z"); // held, but only as a non-ternary integer scheme

		List<Candidate> todo = ImportWitteveenSchemes.plan(LISTING, have, 2, 32);

		assertThat(todo).extracting(Candidate::path)
				.containsExactly("schemes/7x11x15_m772_ZT.json", "schemes/2x13x15_m304_ZT.json");
		assertThat(todo).extracting(Candidate::cls).containsOnly("ZT");
	}

	@Test
	public void plan_is_empty_once_everything_is_held_and_respects_the_dimension_window() {
		Set<String> have = new HashSet<>();
		KnownSchemeKeys.add(have, 15, 11, 13, 1364, "ZT"); // held at another orientation
		KnownSchemeKeys.add(have, 7, 11, 15, 772, "ZT");
		KnownSchemeKeys.add(have, 2, 13, 15, 304, "ZT");
		assertThat(ImportWitteveenSchemes.plan(LISTING, have, 2, 32)).isEmpty();

		assertThat(ImportWitteveenSchemes.plan(LISTING, new HashSet<>(), 2, 14))
				.as("every listed scheme has a dimension of 15 → outside --max-dim=14").isEmpty();
	}
}
