package eu.solven.matmul.papers.khoruzhii2026;

import static org.assertj.core.api.Assertions.assertThat;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.stream.Stream;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.io.NpzSchemeReader;
import eu.solven.matmul.search.LineageReplayer;
import eu.solven.matmul.verifiers.Verifier;

/**
 * The {@code TA_lita_npz(n=N,r=R)} data-backed atom: the KGP 2026 LITA cubes as published
 * ({@code src/main/resources/external/khoruzhii-lita/*.npz}), read without NumPy.
 * A whole band of projection stubs pins these cubes by {@code NxNxN@contentHash}, so
 * the reader must be deterministic and convention-exact — a transposed axis or a
 * mis-scaled rational still yields a scheme of the right shape and rank, it just no
 * longer computes matmul (or hashes differently and orphans every dependent stub).
 */
public class TestLitaNpzCubes {

	@Test
	public void reader_yields_an_exact_13x13x13_at_rank_1379() throws Exception {
		File npz = LitaNpzCubes.file(13).orElseThrow().toFile();
		NpzSchemeReader.Loaded loaded = NpzSchemeReader.read(npz);
		NonCubicBilinearAlgorithm alg = loaded.alg();
		assertThat(new int[] { alg.n, alg.m, alg.p }).containsExactly(13, 13, 13);
		assertThat(alg.r).isEqualTo(1379);
		assertThat(loaded.coefficientField()).isEqualTo("Q");
		// Exact BigInteger proof (≈14M sparse terms, a few seconds): the convention check.
		assertThat(Verifier.isExactNonCubic(alg)).as("13³ must compute matmul exactly").isTrue();
		// Rational with an even denominator → not reducible mod 2 as written.
		assertThat(loaded.lcmDenominator().testBit(0)).as("lcm of denominators is even").isFalse();
	}

	@Test
	public void rank_formulas_match_every_archive() throws Exception {
		List<Path> files;
		try (Stream<Path> ls = Files.list(LitaNpzCubes.DIR)) {
			files = ls.filter(p -> p.toString().endsWith(".npz")).toList();
		}
		assertThat(files).as("the 20 September cubes 13 ≤ N ≤ 32 plus the 10 October even ones").hasSize(30);
		for (int n = 13; n <= 32; n++) {
			// September formulas: R_even = N³/3 + 3N² + 37N/6 + 5 ; R_odd = N³/3 + 7N²/2 + 14N/3 − 11/2
			// (README). Since 2026-10-07 ("Improve even LITA construction") every even cube is one
			// below R_even; both archives are kept and the lowest-rank one is the unqualified atom.
			long num = (n % 2 == 0) ? (2L * n * n * n + 18L * n * n + 37L * n + 30) : (2L * n * n * n + 21L * n * n + 28L * n - 33);
			assertThat(num % 6).as("formula integrality at N=%d", n).isZero();
			long september = num / 6;
			assertThat(LitaNpzCubes.file(n, (int) september)).as("September archive for N=%d", n).isPresent();
			long best = n % 2 == 0 ? september - 1 : september;
			String name = LitaNpzCubes.file(n).orElseThrow().getFileName().toString();
			assertThat(name).as("lowest-rank archive for N=%d", n).isEqualTo(n + "x" + n + "x" + n + "_r" + best + ".npz");
		}
	}

	/** Two archives for one cube: each rank-qualified atom replays to ITS archive, bit-exact. */
	@Test
	public void rank_qualified_atoms_keep_both_even_constructions() throws Exception {
		assertThat(LineageReplayer.isParametricRef("TA_lita_npz(n=14,r=1594)")).isTrue();
		assertThat(LitaNpzCubes.file(14, 1594)).isPresent();
		assertThat(LitaNpzCubes.file(14, 1593)).isPresent();
		assertThat(LitaNpzCubes.file(14, 1500)).isEmpty();
		var replayer = LineageReplayer.withDefaultPool(new FieldAwareLookup("Q"));
		NonCubicBilinearAlgorithm september = replayer.replay(new eu.solven.matmul.catalog.Lineage.Atom("TA_lita_npz(n=14,r=1594)"));
		NonCubicBilinearAlgorithm october = replayer.replay(new eu.solven.matmul.catalog.Lineage.Atom("TA_lita_npz(n=14,r=1593)"));
		NonCubicBilinearAlgorithm legacy = replayer.replay(new eu.solven.matmul.catalog.Lineage.Atom("TA_lita_npz(n=14)"));
		assertThat(september.r).isEqualTo(1594);
		assertThat(october.r).isEqualTo(1593);
		assertThat(legacy.r).as("the unqualified atom is the lowest-rank archive").isEqualTo(1593);
		assertThat(SchemeIO.contentHash(october)).isNotEqualTo(SchemeIO.contentHash(september));
	}

	@Test
	public void atom_is_parametric_and_replays_to_the_stubs_stamped_hash() throws Exception {
		assertThat(LineageReplayer.isParametricRef("TA_lita_npz(n=13)")).isTrue();
		assertThat(LineageReplayer.isParametricRef("TA_lita(n=19)")).isTrue();
		Path dir = Path.of("src/main/resources/schemes/known/section13");
		Path stub;
		try (Stream<Path> ls = Files.list(dir)) {
			stub = ls.filter(p -> p.getFileName().toString().startsWith("13x13x13-r1379-khoruzhii_2026-"))
					.findFirst().orElseThrow();
		}
		var root = SchemeIO.parseJson(stub.toFile());
		assertThat(SchemeIO.isStub(root)).as("the cube is a stub, not 117 MB of explicit JSON").isTrue();
		assertThat(root.get("lineage").get("ref").asString()).isEqualTo("TA_lita_npz(n=13,r=1379)");
		NonCubicBilinearAlgorithm replayed =
				LineageReplayer.withDefaultPool(new FieldAwareLookup("Q")).replayFromFile(stub.toFile());
		assertThat(SchemeIO.contentHash(replayed))
				.as("replay must rebuild exactly the hash projection stubs pin as their parent")
				.isEqualTo(SchemeIO.readHash(root));
		assertThat(stub.getFileName().toString()).contains(SchemeIO.readHash(root).substring(0, 7));
	}
}
