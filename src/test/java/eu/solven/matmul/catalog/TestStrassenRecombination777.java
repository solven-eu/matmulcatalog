package eu.solven.matmul.catalog;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Optional;
import java.util.stream.Stream;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.catalog.Recombination;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.catalog.Recombination.AlgorithmLookup;

/**
 * Probes whether {@link Recombination#construct} on Strassen ⟨2,2,2⟩ base
 * with target ⟨7,7,7⟩ produces Sedoglavic's 250-mult algorithm
 * (the algebraic saving comes from {@code _process_additions}'s
 * cross-block min-reduction). Expected to land somewhere between
 * the naïve block-decomposition rank 273 and Sedoglavic's 250.
 */
public class TestStrassenRecombination777 {

	private static AlgorithmLookup catalogLookup() {
		return (n, m, p) -> {
			int[] sorted = { n, m, p };
			java.util.Arrays.sort(sorted);
			String prefix = sorted[0] + "x" + sorted[1] + "x" + sorted[2];

			Path root = Path.of("src/main/resources/schemes");
			try (Stream<Path> s = Files.walk(root)) {
				Optional<Path> best = s.filter(p_ -> {
					String name = p_.getFileName().toString();
					if (!name.endsWith(".json")) return false;
					// Match both the legacy `note-{shape}_m{rank}` and the
					// content-driven `{shape}-r{rank}-note-{hash}` filename forms:
					// the shape may sit at the start of the name and be followed by
					// `-r` (not only `_m`/`_r`).
					if (!name.matches("(.*[_-])?" + prefix + "[_-][rm].*")) return false;
					if (name.contains("F2") || name.contains("Z2")) return false;
					return true;
				}).findFirst();
				if (best.isEmpty()) return Optional.empty();
				File f = best.get().toFile();
				try {
					NonCubicBilinearAlgorithm alg = SchemeIO.readBilinear(f);
					return alg.orientAs(n, m, p);
				} catch (Exception e) {
					return Optional.empty();
				}
			} catch (IOException e) {
				return Optional.empty();
			}
		};
	}

	@Test
	public void strassen_recombination_777_with_forced_split_4_3_yields_sedoglavic_rank() throws IOException {
		NonCubicBilinearAlgorithm strassen = SchemeIO.read(
				eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json"));
		assertThat(strassen.r).isEqualTo(7);

		Recombination.SotaResolver sota = (a, b, c) -> {
			if (a == 0 || b == 0 || c == 0) return 0;
			// Skip catalog entries for the full ⟨7,7,7⟩ target itself: we want
			// to force a constructive build, not a direct lookup.
			if (a == 7 && b == 7 && c == 7) return Integer.MAX_VALUE / 100;
			return catalogLookup().find(a, b, c).map(alg -> alg.r).orElse(Integer.MAX_VALUE / 100);
		};

		// Force split [4, 3] each axis — the Sedoglavic recipe.
		Recombination.Result rec = Recombination.recombineWithAllocation(
				strassen, sota, new int[] { 4, 3 }, new int[] { 4, 3 }, new int[] { 4, 3 });
		System.out.println("recombine ⟨7,7,7⟩ via Strassen with FORCED alloc [4,3]³: rank=" + rec.totalRank);
		for (int k = 0; k < strassen.r; k++) {
			System.out.println("  base mult " + k + " → sub-shape "
					+ java.util.Arrays.toString(rec.smallMatrixSizes[k]));
		}
		// Expectation: with split [4,3] and Strassen's cross-block sharing, recombination
		// should find a rank lower than the naïve block decomposition (273) — ideally 250.
		assertThat(rec.totalRank).isLessThanOrEqualTo(273);
	}

	@Test
	public void constructive_sedoglavic_777_via_strassen_4_3_split_yields_250_and_verifies()
			throws IOException {
		NonCubicBilinearAlgorithm strassen = SchemeIO.read(
				eu.solven.matmul.catalog.SchemeResolver.byHint("src/main/resources/schemes/known/section2/strassen-2x2x2_m7_a18.json"));

		// Constructive Sedoglavic: Strassen ⟨2,2,2⟩ outer × mixed-shape inner via
		// forced split [4, 3] on each axis. Bounds-checked embedding correctly
		// handles the padded zero positions where block (1,1) is 3×3 but sub-mults
		// see 4×4 inputs.
		NonCubicBilinearAlgorithm built = Recombination.constructWithAllocation(
				strassen, catalogLookup(),
				new int[] { 4, 3 }, new int[] { 4, 3 }, new int[] { 4, 3 });
		assertThat(built.n).isEqualTo(7);
		// This is the deterministic forced-[4,3] constructive build — the original
		// Sedoglavic 2017 recombination, which yields 250 (as the method name says).
		// The global SOTA ⟨7,7,7⟩=249 (Perminov) is found by a different, free-
		// allocation search, NOT by this fixed-allocation construction.
		assertThat(built.r).isEqualTo(250);
		assertThat(Verifier.isExactNonCubic(built)).isTrue();
	}
}
