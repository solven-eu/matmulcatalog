package eu.solven.matmul.docs.migrate;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.verifiers.Verifier;
import eu.solven.matmul.catalog.Lineage;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.papers.khoruzhii2026.LitaTaConstruction;
import lombok.extern.slf4j.Slf4j;

/**
 * Writes the odd-N LITA cubic lineage-only stubs to {@code schemes/derived/}.
 *
 * <p>Each stub is {@code {NxNxN, m=rank, lineage=Atom("TA_lita(n=N)"),
 * scheme_type="stub", fields=[Q,R,C]}} — the explicit factors are reproduced on
 * demand by {@link eu.solven.matmul.search.LineageReplayer} (which maps
 * {@code TA_lita(n=N)} → {@link LitaTaConstruction#build}). Derived (regenerable),
 * NOT known — {@code known/} is reserved for dense/non-reconstructible imports.</p>
 *
 * <p>Direct + lightweight on purpose: it bypasses the full materialise pipeline
 * (projection-parent PanTA cubes ~300&nbsp;MB each → OOM on the large cubics).
 * It builds one cube at a time, verifies via {@link Verifier#verifyAuto} — an
 * EXACT algebraic proof for the sparse odd cubes, and a random spot-check for the
 * dense even cubes (whose exact term-map would OOM) — and writes. The catalog's
 * own verifiers ({@code VerifyScheme}, {@code VerifyAllSchemes}) spot-check char-0
 * schemes too, so the dense even stubs verify there without OOM.</p>
 *
 * <p>Covers every cubic where LITA beats our catalog: odd {@code 19..31} and even
 * {@code 26,28,30,32} (even {@code 20,22,24} lose to TA_dis, so they're skipped).</p>
 *
 * <pre>mvn -q exec:java -Dexec.mainClass=eu.solven.matmul.docs.migrate.MaterialiseLitaCubes</pre>
 */
@Slf4j
public final class MaterialiseLitaCubes {

	private MaterialiseLitaCubes() {}

	/** Cubics where LITA beats the catalog (odd 19..31 + even 26/28/30/32). */
	private static final int[] N_LIST = { 19, 21, 23, 25, 26, 27, 28, 29, 30, 31, 32 };
	private static final Path ROOT = Path.of("src/main/resources/schemes/derived");

	public static void main(String[] args) throws Exception {
		int written = 0;
		for (int n : N_LIST) {
			NonCubicBilinearAlgorithm alg = LitaTaConstruction.build(n);
			if (alg.r != (int) eu.solven.matmul.papers.khoruzhii2026.LitaTrilinearAggregation.cubicRank(n)) {
				throw new IllegalStateException("⟨" + n + "³⟩ built rank " + alg.r + " != formula");
			}
			// verifyAuto: exact proof for sparse (odd), random spot-check for dense (even).
			Verifier.Verdict verdict = Verifier.verifyAuto(alg);
			if (!verdict.ok()) {
				log.warn("LITA ⟨{}³⟩=r{} did NOT verify ({}) — skipping", n, alg.r, verdict.strategy());
				continue;
			}
			Lineage.Node lineage = new Lineage.Atom("TA_lita(n=" + n + ")");
			// Q-native (÷12 rational): narrow [Q,R,C] to what the coefficients support
			// (drops nothing here — rational is valid over Q⊂R⊂C; not Z/F2/F3).
			List<String> fields = SchemeIO.narrowFieldsToCoefficients(alg, List.of("Q", "R", "C"));
			String hash7 = SchemeIO.contentHash(alg).substring(0, 7);
			String fname = n + "x" + n + "x" + n + "-r" + alg.r + "-derived-" + hash7 + ".json";
			Path dir = ROOT.resolve("section" + n);
			Files.createDirectories(dir);
			File f = dir.resolve(fname).toFile();
			SchemeIO.writeStub(alg, f, lineage, fields);
			written++;
			log.info("wrote STUB {} (r={}, fields={}, verify={}, lineage=TA_lita(n={}))",
					fname, alg.r, fields, verdict.strategy(), n);
		}
		log.info("MaterialiseLitaCubes: wrote {} LITA cubic stubs (odd exact + even spot-check)", written);
	}
}
