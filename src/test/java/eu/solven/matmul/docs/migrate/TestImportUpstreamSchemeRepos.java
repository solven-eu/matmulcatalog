package eu.solven.matmul.docs.migrate;

import static org.assertj.core.api.Assertions.assertThat;

import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.docs.migrate.ImportUpstreamSchemeRepos.Candidate;

/**
 * The mizopl / we*star channel (issue #20): two directories, two layouts, one pattern.
 * The addition-reduced files must be planned even when their (shape, rank, class) is
 * held — what they carry is the program — and must stop being planned once that
 * program (by source) is attached.
 */
public class TestImportUpstreamSchemeRepos {
	private static final String DIR = "solutions/fast-matrix-multiplication/";
	private static final List<String> MIZOPL_LISTING = List.of(
			DIR + "README.md",
			DIR + "status.json",
			DIR + "zt_rank_schemes/CERTIFICATE.json",
			DIR + "zt_rank_schemes/9x12x14_m954_ZT.json",
			DIR + "zt_rank_schemes/12x13x16_m1541_Z.json",
			DIR + "slp_addition_reduced/clss_cl1_tellegen_reducer.cpp",
			DIR + "slp_addition_reduced/3x3x7_m49_cr167_cn404_ZT_reduced.json",
			DIR + "slp_addition_reduced/3x5x6_m68_cr245_cn494_ZT_reduced.json");

	@Test
	public void mizopl_pattern_reads_both_directories_and_flags_the_programs() {
		Candidate rank = UpstreamSchemeRepos.MIZOPL.parse(DIR + "zt_rank_schemes/9x12x14_m954_ZT.json").orElseThrow();
		assertThat(rank.n()).isEqualTo(9);
		assertThat(rank.m()).isEqualTo(12);
		assertThat(rank.p()).isEqualTo(14);
		assertThat(rank.rank()).isEqualTo(954);
		assertThat(rank.cls()).isEqualTo(KnownSchemeKeys.ZT);
		assertThat(rank.reduced()).isFalse();

		Candidate z = UpstreamSchemeRepos.MIZOPL.parse(DIR + "zt_rank_schemes/12x13x16_m1541_Z.json").orElseThrow();
		assertThat(z.cls()).isEqualTo(KnownSchemeKeys.Z);
		assertThat(z.reduced()).isFalse();

		Candidate slp = UpstreamSchemeRepos.MIZOPL
				.parse(DIR + "slp_addition_reduced/3x3x7_m49_cr167_cn404_ZT_reduced.json").orElseThrow();
		assertThat(slp.rank()).isEqualTo(49);
		assertThat(slp.cls()).isEqualTo(KnownSchemeKeys.ZT);
		assertThat(slp.reduced()).as("the _reduced suffix marks an addition-reduced program").isTrue();

		assertThat(UpstreamSchemeRepos.MIZOPL.parse(DIR + "zt_rank_schemes/CERTIFICATE.json")).isEmpty();
		assertThat(UpstreamSchemeRepos.MIZOPL.parse(DIR + "slp_addition_reduced/clss_cl1_tellegen_reducer.cpp")).isEmpty();
		assertThat(UpstreamSchemeRepos.MIZOPL.parse(DIR + "status.json")).isEmpty();
	}

	@Test
	public void programs_are_planned_over_held_matrices_until_their_source_is_attached() {
		Set<String> have = new HashSet<>();
		// We hold ⟨3,3,7⟩=49 ternary (Perminov) and ⟨9,12,14⟩=954 ternary (say, via the mirror).
		KnownSchemeKeys.add(have, 3, 3, 7, 49, KnownSchemeKeys.ZT);
		KnownSchemeKeys.add(have, 9, 12, 14, 954, KnownSchemeKeys.ZT);

		List<Candidate> todo = ImportUpstreamSchemeRepos.plan(UpstreamSchemeRepos.MIZOPL, MIZOPL_LISTING, have, Map.of(),
				2, 32);
		assertThat(todo).extracting(Candidate::path).containsExactlyInAnyOrder(
				DIR + "zt_rank_schemes/12x13x16_m1541_Z.json", // rank scheme not held
				DIR + "slp_addition_reduced/3x3x7_m49_cr167_cn404_ZT_reduced.json", // held matrices, new program
				DIR + "slp_addition_reduced/3x5x6_m68_cr245_cn494_ZT_reduced.json"); // not held at all

		// Once the ⟨3,3,7⟩ program from this source is attached, it is no longer planned;
		// a program from ANOTHER source at the same key does not count.
		Map<String, Set<String>> slpHave = Map.of(
				KnownSchemeKeys.key(3, 3, 7, 49, KnownSchemeKeys.ZT), Set.of(UpstreamSchemeRepos.MIZOPL.source()),
				KnownSchemeKeys.key(3, 5, 6, 68, KnownSchemeKeys.ZT), Set.of("Perminov 2023"));
		todo = ImportUpstreamSchemeRepos.plan(UpstreamSchemeRepos.MIZOPL, MIZOPL_LISTING, have, slpHave, 2, 32);
		assertThat(todo).extracting(Candidate::path).containsExactlyInAnyOrder(
				DIR + "zt_rank_schemes/12x13x16_m1541_Z.json",
				DIR + "slp_addition_reduced/3x5x6_m68_cr245_cn494_ZT_reduced.json");

		// The rank-only overload (the two older repositories) never plans a program by mistake:
		// the Witteveen / Cloudroom patterns have no SLP group.
		assertThat(UpstreamSchemeRepos.WITTEVEEN.slpGroup()).isZero();
		assertThat(UpstreamSchemeRepos.CLOUDROOM.slpGroup()).isZero();
	}

	@Test
	public void registry_lists_mizopl_with_its_attribution_key() {
		assertThat(UpstreamSchemeRepos.ALL).extracting(UpstreamSchemeRepos.UpstreamSchemeRepo::key)
				.contains("mizopl").doesNotHaveDuplicates();
		assertThat(UpstreamSchemeRepos.byKey("mizopl").source()).isEqualTo("mizopl 2026");
		assertThat(UpstreamSchemeRepos.byKey("mizopl").citation()).contains("MIT");
		// The same key PerminovKnownAttribution answers for the mirror, so the two never diverge.
		assertThat(eu.solven.matmul.catalog.PerminovKnownAttribution.forPath("schemes/known/mizopl/x.json")
				.orElseThrow().source()).isEqualTo(UpstreamSchemeRepos.MIZOPL.source());
	}
}
