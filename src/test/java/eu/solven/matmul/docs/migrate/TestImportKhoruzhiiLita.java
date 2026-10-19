package eu.solven.matmul.docs.migrate;

import static org.assertj.core.api.Assertions.assertThat;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

/**
 * {@code ImportKhoruzhiiLita --download}: a new rank for a cube is a NEW file upstream
 * ({@code {N}x{N}x{N}_r{R}.npz}), so the archive name is the idempotence key — the
 * October even cubes ({@code 14x14x14_r1593.npz} next to our {@code _r1594}) must be
 * fetched, everything held must not, and non-archive entries are ignored.
 */
public class TestImportKhoruzhiiLita {

	@Test
	public void missing_archives_are_the_upstream_names_not_on_disk(@TempDir Path dir) throws Exception {
		Files.createFile(dir.resolve("14x14x14_r1594.npz"));
		Files.createFile(dir.resolve("15x15x15_r1977.npz"));
		List<String> upstream = List.of("14x14x14_r1594.npz", "14x14x14_r1593.npz", "15x15x15_r1977.npz",
				"16x16x16_r2236.npz", "README.md", "lita.py");

		assertThat(ImportKhoruzhiiLita.missingArchives(upstream, dir))
				.containsExactly("14x14x14_r1593.npz", "16x16x16_r2236.npz");
	}

	@Test
	public void the_vendored_directory_holds_every_cube_of_both_generations() {
		List<String> upstream = new java.util.ArrayList<>();
		for (int n = 13; n <= 32; n++) {
			long num = (n % 2 == 0) ? (2L * n * n * n + 18L * n * n + 37L * n + 30) : (2L * n * n * n + 21L * n * n + 28L * n - 33);
			long september = num / 6;
			upstream.add(n + "x" + n + "x" + n + "_r" + september + ".npz");
			if (n % 2 == 0) upstream.add(n + "x" + n + "x" + n + "_r" + (september - 1) + ".npz");
		}
		assertThat(ImportKhoruzhiiLita.missingArchives(upstream, eu.solven.matmul.papers.khoruzhii2026.LitaNpzCubes.DIR))
				.as("nothing the catalog already vendors is re-downloaded").isEmpty();
	}
}
