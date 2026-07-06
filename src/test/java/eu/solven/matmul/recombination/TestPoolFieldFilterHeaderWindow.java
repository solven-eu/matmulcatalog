package eu.solven.matmul.recombination;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

/**
 * Guards {@link BlockSplitSearch#isFieldValidLeafNC} against the PAYLOAD-FIRST
 * header-scan leak: imports whose JSON key order puts the u/v/w payload before
 * the metadata (e.g. the {@code flips_mod2} family) keep {@code fields[]} and
 * {@code commutative} at the file TAIL, outside the 4 KB header window. The
 * filter used to treat "no fields[] in the header" as a legacy all-fields atom
 * and silently admitted F₂/C-only bases into char-0 recombination pools — the
 * verify-gate rejected the bogus compositions downstream (e.g. a phantom
 * ⟨20,25,30⟩=8200 over an F₂ ⟨4,4,5⟩=60 base), but the pool stayed polluted.
 *
 * <p>Sibling of {@code eu.solven.matmul.search.TestPoolFieldFilter}, which guards
 * the same filter's earlier singular-"field"-key leak at pool level.</p>
 */
public class TestPoolFieldFilterHeaderWindow {

	/** The real committed regression file: F₂-only, fields[] at byte ~11570. */
	private static final File FLIPS_445 = new File(
			"src/main/resources/schemes/known/section5/4x4x5-r60-flips_mod2-e7a8ee8.json");

	@Test
	public void real_payload_first_f2_file_rejected_from_char0_pools() {
		assertThat(FLIPS_445).exists();
		assertThat(BlockSplitSearch.isFieldValidLeafNC(FLIPS_445, 8, "R")).isFalse();
		assertThat(BlockSplitSearch.isFieldValidLeafNC(FLIPS_445, 8, "Q")).isFalse();
		// … but stays a valid base for an F₂ sweep.
		assertThat(BlockSplitSearch.isFieldValidLeafNC(FLIPS_445, 8, "F2")).isTrue();
	}

	@Test
	public void synthetic_payload_first_f2_rejected(@TempDir Path dir) throws Exception {
		File f = payloadFirst(dir, "\"fields\": [\"F2\"]");
		assertThat(BlockSplitSearch.isFieldValidLeafNC(f, 8, "R")).isFalse();
		assertThat(BlockSplitSearch.isFieldValidLeafNC(f, 8, "F2")).isTrue();
	}

	@Test
	public void synthetic_payload_first_commutative_rejected(@TempDir Path dir) throws Exception {
		// Same window bug, other marker: commutative-only schemes do not lift to NC
		// recombination bases, even when the marker sits after the payload.
		File f = payloadFirst(dir, "\"fields\": [\"Z\"],\n  \"commutative\": true");
		assertThat(BlockSplitSearch.isFieldValidLeafNC(f, 8, "R")).isFalse();
	}

	@Test
	public void header_fields_keep_the_fast_path(@TempDir Path dir) throws Exception {
		// Canonical layout (fields[] in the header): a Z base is valid over R via the
		// inclusion chain Z ⊂ Q ⊂ R, and an F₂-only base is not.
		Path z = dir.resolve("z.json");
		Files.writeString(z, "{\n  \"n\": [2, 2, 2],\n  \"m\": 7,\n  \"fields\": [\"Z\"],\n  \"u\": []\n}");
		assertThat(BlockSplitSearch.isFieldValidLeafNC(z.toFile(), 8, "R")).isTrue();
		Path f2 = dir.resolve("f2.json");
		Files.writeString(f2, "{\n  \"n\": [2, 2, 2],\n  \"m\": 7,\n  \"fields\": [\"F2\"],\n  \"u\": []\n}");
		assertThat(BlockSplitSearch.isFieldValidLeafNC(f2.toFile(), 8, "R")).isFalse();
	}

	/** A >4 KB file whose metadata trails the payload, mimicking the flips imports. */
	private static File payloadFirst(Path dir, String tailMetadata) throws Exception {
		StringBuilder sb = new StringBuilder("{\n  \"n\": [2, 2, 2],\n  \"m\": 7,\n  \"u\": [\n");
		for (int i = 0; i < 600; i++) {
			sb.append("    [1, 0, 0, 1, 0, 0, 1, 0],\n");
		}
		sb.append("    [1, 0, 0, 1, 0, 0, 1, 0]\n  ],\n  ").append(tailMetadata).append("\n}");
		Path f = dir.resolve("payload-first.json");
		Files.writeString(f, sb.toString());
		assertThat(f.toFile().length()).isGreaterThan(4096L);
		return f.toFile();
	}
}
