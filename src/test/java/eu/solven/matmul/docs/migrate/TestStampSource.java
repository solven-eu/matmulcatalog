package eu.solven.matmul.docs.migrate;

import static org.assertj.core.api.Assertions.assertThat;

import org.junit.jupiter.api.Test;

/**
 * Guards the 2026-07-07 silent bug: {@link StampSource} on a NEW-convention
 * filename ({@code {n}x{m}x{p}-r{rank}-{note}-{hash7}}) found no old-convention
 * source prefix, fell back to the whole stem, and stamped 7.6k filename echoes
 * like {@code "source": "10x3x3-R69-Derived-79abb2c"} — wrong attribution, no
 * crash. New-convention stems must never be stamped.
 */
public class TestStampSource {

	@Test
	public void new_convention_stems_are_never_stamped() {
		assertThat(StampSource.sourceForStem("10x3x3-r69-derived-79abb2c")).isEmpty();
		assertThat(StampSource.sourceForStem("5x5x12-r204-perminov_ZT-61a6cb7")).isEmpty();
		assertThat(StampSource.sourceForStem("20x20x25-r5611-derived-0000000")).isEmpty();
	}

	@Test
	public void old_convention_prefixes_still_map_to_attributions() {
		assertThat(StampSource.sourceForStem("alphatensor-Z_2x3x3_r15_a58")).contains("AlphaTensor 2022");
		assertThat(StampSource.sourceForStem("perminov-ZT_5x5x5_r93_a843")).contains("Perminov 2023");
		assertThat(StampSource.sourceForStem("fmm-lille_5x5x12_r204_a2326")).contains("FMM-Lille");
	}

	@Test
	public void echo_detector_matches_exactly_the_bogus_values() {
		assertThat(FixFilenameEchoSources.isFilenameEcho("10x3x3-R69-Derived-79abb2c")).isTrue();
		assertThat(FixFilenameEchoSources.isFilenameEcho("AlphaTensor 2022")).isFalse();
		assertThat(FixFilenameEchoSources.isFilenameEcho("Derived_recursive")).isFalse();
		assertThat(FixFilenameEchoSources.isFilenameEcho("Strassen 1969")).isFalse();
	}
}
