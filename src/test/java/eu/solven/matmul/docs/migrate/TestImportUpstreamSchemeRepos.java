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

	/**
	 * The four sources of 2026-11 (Perminov's issue #18 and his three unmapped mirrors). Two
	 * publish gzipped rank-row JSON and far more schemes above dimension 16 than the catalog
	 * can hold as explicit files, so they carry an explicit-import cap: above it a scheme is
	 * listed in the rank digest, never planned — and the digest keeps the lowest rank per
	 * sorted shape so the scan still sees it.
	 */
	@Test
	public void sidebottom_and_he_are_capped_at_16_and_listed_above() {
		List<String> sidebottom = List.of("README.md", "schemes/9x11x11_705_dyadic/9x11x11_m705_Q.json.gz",
				"schemes/9x11x11_705_dyadic/9x11x11_705_dyadic_L.slp",
				"schemes/13x14x16_1805_ternary/13x14x16_m1805_ZT.json.gz",
				"large/15x16x16_2193_rational/15x16x16_m2193_Q.json.gz",
				"large/14x16x17_2278_rational/14x16x17_m2278_Q.json.gz",
				"large/14x16x17_2278_rational/path_1.json",
				"products/10x10x21_1299_integer/10x10x21_m1299.json.gz");
		UpstreamSchemeRepos.UpstreamSchemeRepo repo = UpstreamSchemeRepos.SIDEBOTTOM;
		assertThat(repo.rankRows()).isTrue();
		assertThat(repo.explicitMaxDim()).isEqualTo(16);
		assertThat(repo.parse("schemes/13x14x16_1805_ternary/13x14x16_m1805_ZT.json.gz").orElseThrow().cls())
				.isEqualTo(KnownSchemeKeys.ZT);
		assertThat(repo.parse("products/10x10x21_1299_integer/10x10x21_m1299.json.gz").orElseThrow().cls())
				.as("no tag: the coarsest class is assumed").isEqualTo(KnownSchemeKeys.Q);
		assertThat(repo.parse("large/14x16x17_2278_rational/path_1.json")).isEmpty();

		List<Candidate> todo = ImportUpstreamSchemeRepos.plan(repo, sidebottom, new HashSet<>(), Map.of(), 2, 32);
		assertThat(todo).extracting(Candidate::path).containsExactlyInAnyOrder(
				"schemes/9x11x11_705_dyadic/9x11x11_m705_Q.json.gz",
				"schemes/13x14x16_1805_ternary/13x14x16_m1805_ZT.json.gz",
				"large/15x16x16_2193_rational/15x16x16_m2193_Q.json.gz");

		UpstreamRepoRanks.RepoRanks digest = UpstreamRepoRanks.of(repo, sidebottom);
		assertThat(digest.ranks()).as("the capped schemes are still listed")
				.containsEntry(eu.solven.matmul.Shape.of(14, 16, 17).canonical(), 2278)
				.containsEntry(eu.solven.matmul.Shape.of(10, 10, 21).canonical(), 1299).hasSize(5);

		UpstreamSchemeRepos.UpstreamSchemeRepo he = UpstreamSchemeRepos.HE_OBLIQUE;
		Candidate dyadic = he.parse("schemes_classic_le16/9x11x11_m702_dyadic.json.gz").orElseThrow();
		assertThat(dyadic.rank()).isEqualTo(702);
		assertThat(dyadic.cls()).as("'dyadic' names the denominators, it is not an integer class")
				.isEqualTo(KnownSchemeKeys.Q);
		assertThat(ImportUpstreamSchemeRepos.plan(he, List.of("schemes_classic_le16/9x11x11_m702_dyadic.json.gz",
				"schemes_large/20x21x22_m4897_Q.json.gz", "verify_all.py"), new HashSet<>(), Map.of(), 2, 32))
				.extracting(Candidate::path).containsExactly("schemes_classic_le16/9x11x11_m702_dyadic.json.gz");
	}

	@Test
	public void rank_digest_keeps_the_lowest_rank_per_sorted_shape_and_round_trips(
			@org.junit.jupiter.api.io.TempDir java.nio.file.Path tmp) throws Exception {
		List<String> listing = List.of("schemes/15x16x16_2200_rational/15x16x16_m2200_Q.json.gz",
				"large/15x16x16_2193_rational/15x16x16_m2193_Q.json.gz",
				"large/14x16x17_2278_rational/14x16x17_m2278_Q.json.gz");
		java.nio.file.Path file = tmp.resolve("ranks.json");
		UpstreamRepoRanks.merge(file, Map.of(UpstreamSchemeRepos.SIDEBOTTOM, listing));
		// a later run that lists another repository only must keep this one
		UpstreamRepoRanks.merge(file, Map.of(UpstreamSchemeRepos.HE_OBLIQUE,
				List.of("schemes_classic_le16/15x16x16_m2173_Q.json.gz")));
		Map<String, UpstreamRepoRanks.RepoRanks> back = UpstreamRepoRanks.read(file);
		assertThat(back).containsOnlyKeys("sidebottom", "he-oblique");
		var s151616 = eu.solven.matmul.Shape.of(15, 16, 16).canonical();
		assertThat(back.get("sidebottom").ranks()).containsEntry(s151616, 2193).hasSize(2);
		assertThat(back.get("sidebottom").paths().get(s151616)).startsWith("large/");
		assertThat(back.get("sidebottom").explicitMaxDim()).isEqualTo(16);
		assertThat(back.get("he-oblique").ranks()).containsEntry(s151616, 2173);
		String once = java.nio.file.Files.readString(file);
		UpstreamRepoRanks.merge(file, Map.of());
		assertThat(java.nio.file.Files.readString(file)).as("an unchanged upstream rewrites the file byte for byte")
				.isEqualTo(once);
	}

	@Test
	public void perminov_mirrors_of_2026_11_are_mapped_and_pulled_by_a_registry_entry() {
		// known/arhancanli and known/stapleton have no JSON origin: the registry reads the mirror itself.
		Candidate canli = UpstreamSchemeRepos.CANLI.parse("schemes/known/arhancanli/9x10x11_m640.json").orElseThrow();
		assertThat(canli.rank()).isEqualTo(640);
		Candidate slp = UpstreamSchemeRepos.STAPLETON
				.parse("schemes/known/stapleton/additions_reduced_ZT/2x4x10_m64_cr207_cn489_ZT_reduced.json").orElseThrow();
		assertThat(slp.reduced()).isTrue();
		assertThat(slp.cls()).isEqualTo(KnownSchemeKeys.ZT);
		assertThat(UpstreamSchemeRepos.CANLI.parse("schemes/known/gsidebottom/13x14x14_m1564_Q.json")).isEmpty();
		for (var e : Map.of("gsidebottom", UpstreamSchemeRepos.SIDEBOTTOM, "arhancanli", UpstreamSchemeRepos.CANLI,
				"stapleton", UpstreamSchemeRepos.STAPLETON).entrySet()) {
			assertThat(eu.solven.matmul.catalog.PerminovKnownAttribution.isMappedSubtree(e.getKey())).isTrue();
			assertThat(eu.solven.matmul.catalog.PerminovKnownAttribution.forPath("schemes/known/" + e.getKey() + "/x.json")
					.orElseThrow().source()).isEqualTo(e.getValue().source());
		}
		assertThat(UpstreamSchemeRepos.ALL).extracting(UpstreamSchemeRepos.UpstreamSchemeRepo::key).doesNotHaveDuplicates();
		assertThat(UpstreamSchemeRepos.ALL).extracting(UpstreamSchemeRepos.UpstreamSchemeRepo::note).doesNotHaveDuplicates();
	}
}
