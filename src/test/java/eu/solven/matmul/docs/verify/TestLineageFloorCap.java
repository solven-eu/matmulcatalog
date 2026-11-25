package eu.solven.matmul.docs.verify;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.Lineage;
import eu.solven.matmul.catalog.SchemeIO;

/**
 * A stub claims what its atoms prove ({@link VerifyLineageFieldCompat#cappedAtLineageFloor}).
 *
 * <p>2026-11-25: two stub writers stamped {@code Z} on stubs whose lineage ran through a base
 * declaring only {@code Q} — {@code MaterialiseSubstitutionProjections} from the built (integer)
 * coefficients, {@code RestampStubRanks} by keeping the old stub's fields after its bare-shape
 * leaves had moved to newly imported Q-only bases. Nothing failed locally; the PR's
 * {@code VerifySchemes --check=lineage} refused four stubs. Both writers now go through the cap
 * tested here, against a two-file catalog (no real catalog needed).</p>
 */
public class TestLineageFloorCap {

	/** Naive ⟨1,1,2⟩ (rank 2); {@code half}: the two factors of product 0 are scaled by ½ and 2. */
	private static NonCubicBilinearAlgorithm naive112(boolean half) {
		double[][] u = { { half ? 0.5 : 1, 1 } };
		double[][] v = { { half ? 2 : 1, 0 }, { 0, 1 } };
		double[][] w = { { 1, 0 }, { 0, 1 } };
		return new NonCubicBilinearAlgorithm(1, 1, 2, u, v, w);
	}

	private static FieldAwareLookup catalogWith(Path root, boolean rational) throws Exception {
		Path dir = Files.createDirectories(root.resolve("known").resolve("section2"));
		NonCubicBilinearAlgorithm alg = naive112(rational);
		File f = dir.resolve(SchemeIO.canonicalName(alg, "test")).toFile();
		SchemeIO.write(alg, f);
		SchemeIO.updateFields(f, Map.of("fields", rational ? List.of("Q", "R", "C") : List.of("F2", "F3", "Z", "Q", "R", "C"),
				"source", "test"), List.of(), true);
		return new FieldAwareLookup("Q", root);
	}

	@Test
	public void z_is_dropped_when_an_atom_declares_only_q(@TempDir Path root) throws Exception {
		FieldAwareLookup lookup = catalogWith(root, true);
		Lineage.Node leaf = new Lineage.Atom("1x1x2");
		assertThat(VerifyLineageFieldCompat.cappedAtLineageFloor(List.of("F2", "F3", "Z", "Q", "R", "C"), leaf, lookup))
				.as("F2 / F3 are not a lineage fact and stay; Z is below the atoms' floor").containsExactly("F2", "F3", "Q",
						"R", "C");
		assertThat(VerifyLineageFieldCompat.cappedAtLineageFloor(List.of("Q", "R", "C"), leaf, lookup))
				.containsExactly("Q", "R", "C");
	}

	@Test
	public void an_integer_atom_leaves_the_claim_untouched(@TempDir Path root) throws Exception {
		FieldAwareLookup lookup = catalogWith(root, false);
		List<String> all = List.of("F2", "F3", "Z", "Q", "R", "C");
		assertThat(VerifyLineageFieldCompat.cappedAtLineageFloor(all, new Lineage.Atom("1x1x2"), lookup))
				.containsExactlyElementsOf(all);
	}
}
