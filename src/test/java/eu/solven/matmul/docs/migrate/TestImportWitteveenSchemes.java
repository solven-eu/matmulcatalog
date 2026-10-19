package eu.solven.matmul.docs.migrate;

import static org.assertj.core.api.Assertions.assertThat;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.docs.migrate.ImportUpstreamSchemeRepos.Candidate;

/**
 * The scheduled pull from the registered upstream repositories must pick up exactly
 * the files the catalog does not hold yet — by {@code (shape, rank, class)}, not by
 * name — and nothing that is not a scheme. One test per naming convention in the
 * registry.
 */
public class TestImportWitteveenSchemes {

	private static final List<String> WITTEVEEN_LISTING = List.of(
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

		List<Candidate> todo = ImportUpstreamSchemeRepos.plan(UpstreamSchemeRepos.WITTEVEEN, WITTEVEEN_LISTING, have, 2, 32);

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
		assertThat(ImportUpstreamSchemeRepos.plan(UpstreamSchemeRepos.WITTEVEEN, WITTEVEEN_LISTING, have, 2, 32)).isEmpty();

		assertThat(ImportUpstreamSchemeRepos.plan(UpstreamSchemeRepos.WITTEVEEN, WITTEVEEN_LISTING, new HashSet<>(), 2, 14))
				.as("every listed scheme has a dimension of 15 → outside --max-dim=14").isEmpty();
	}

	/**
	 * Cloudroom's names carry no coefficient-class tag ({@code 4x6x15_rank261.json}), so
	 * the "already held" test uses the coarsest class: any scheme at that
	 * {@code (shape, rank)} counts, and {@code schemes/ties/…} (a sub-directory) is not a
	 * release file.
	 */
	@Test
	public void cloudroom_names_plan_at_the_coarsest_class() {
		List<String> listing = List.of("README.md", "CONTENTS.md", "schemes/4x6x15_rank261.json",
				"schemes/9x11x11_rank710.json", "schemes/4x4x4_rank48_seeds.tar.xz", "schemes/ties/6x6x6_rank153.json",
				"third_party/perminov/9x11x11_m711_Q.json");
		Set<String> have = new HashSet<>();
		KnownSchemeKeys.add(have, 9, 11, 11, 710, "Q"); // held as a rational scheme — enough

		List<Candidate> todo = ImportUpstreamSchemeRepos.plan(UpstreamSchemeRepos.CLOUDROOM, listing, have, 2, 32);

		assertThat(todo).extracting(Candidate::path).containsExactly("schemes/4x6x15_rank261.json");
		assertThat(todo.get(0).cls()).isEqualTo("Q");
		assertThat(todo.get(0).rank()).isEqualTo(261);
	}

	@Test
	public void registry_keys_are_unique_and_resolvable() {
		assertThat(UpstreamSchemeRepos.ALL).extracting(UpstreamSchemeRepos.UpstreamSchemeRepo::key).doesNotHaveDuplicates();
		assertThat(UpstreamSchemeRepos.byKey("cloudroom").repo()).isEqualTo("matyasstoch/cloudroom-math");
		assertThat(UpstreamSchemeRepos.WITTEVEEN.treeUrl())
				.isEqualTo("https://api.github.com/repos/MerlijnW70/fmm-schemes/git/trees/main?recursive=1");
	}
}
