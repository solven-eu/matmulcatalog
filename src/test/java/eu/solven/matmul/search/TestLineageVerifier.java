package eu.solven.matmul.search;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;
import java.nio.file.Path;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.verifiers.Verifier;
import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.Lineage;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.verifiers.LineageVerifier;

/**
 * Compositional verification: correct primitives + correct lineages certify;
 * a tampered primitive does not; an unresolvable atom does not.
 */
public class TestLineageVerifier {

	private final FieldAwareLookup committed = new FieldAwareLookup(Field.Q);
	private final LineageVerifier verifier = new LineageVerifier(committed);

	@Test
	public void correct_primitive_certifies(@TempDir Path tmp) throws Exception {
		NonCubicBilinearAlgorithm s222 = committed.find(2, 2, 2).orElseThrow();
		assertThat(Verifier.isExactNonCubic(s222)).isTrue();  // sanity
		File f = tmp.resolve("good-2x2x2_m7.json").toFile();
		SchemeIO.write(s222, f);

		LineageVerifier.Result r = verifier.verifyFile(f);
		assertThat(r.certified()).isTrue();
		assertThat(r.primitivesVerified()).isGreaterThanOrEqualTo(1);
	}

	@Test
	public void tampered_primitive_fails(@TempDir Path tmp) throws Exception {
		NonCubicBilinearAlgorithm s = committed.find(2, 2, 2).orElseThrow();
		// Flip one W entry → no longer a correct matmul scheme.
		double[][] srcU = s.denseU();
		double[][] srcV = s.denseV();
		double[][] srcW = s.denseW();
		double[][] W = new double[srcW.length][srcW[0].length];
		for (int i = 0; i < W.length; i++) W[i] = srcW[i].clone();
		W[0][0] += 1.0;
		NonCubicBilinearAlgorithm broken = new NonCubicBilinearAlgorithm(s.n, s.m, s.p, srcU, srcV, W);
		assertThat(Verifier.isExactNonCubic(broken)).isFalse();  // sanity
		File f = tmp.resolve("broken-2x2x2_m7.json").toFile();
		SchemeIO.write(broken, f);

		LineageVerifier.Result r = verifier.verifyFile(f);
		assertThat(r.certified()).isFalse();
	}

	@Test
	public void serendipitous_lineage_certifies_compositionally() {
		// The smallest serendipitous identity: base ⟨1,1,2⟩ (bud) ⊗ˢ ⟨2,2,1⟩ → ⟨2,2,2⟩.
		// Compositional verification recurses to the base primitive only.
		Lineage.Node base = new Lineage.Atom("1x1x2");
		Lineage.Node tree = new Lineage.SerendipitousProduct(base, 2, 2, 1);
		LineageVerifier.Result r = verifier.verify(tree);
		assertThat(r.certified()).isTrue();
	}

	@Test
	public void unresolvable_atom_does_not_certify() {
		LineageVerifier.Result r = verifier.verify(new Lineage.Atom("999x998x997-nonexistent"));
		assertThat(r.certified()).isFalse();
	}
}
