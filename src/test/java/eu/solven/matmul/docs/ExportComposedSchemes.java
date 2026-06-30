package eu.solven.matmul.docs;

import java.io.File;
import java.io.IOException;

import eu.solven.matmul.ComplexNonCubicBilinearAlgorithm;
import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.verifiers.Verifier;
import eu.solven.matmul.catalog.Compose;
import eu.solven.matmul.catalog.Compositions;
import eu.solven.matmul.catalog.SchemeIO;

/**
 * One-shot exporter for composed catalog schemes — converts the in-memory
 * compositions in {@link Compositions} into persistent JSON files under
 * {@code src/main/resources/schemes/}, so the next coverage-matrix run picks
 * them up as ✓-verified.
 *
 * <p>Filenames follow the canonical convention
 * {@code derived-<recipe>_<n>x<m>x<p>_r<rank>.json}, with {@code recipe}
 * encoding the constituent schemes (e.g. {@code strassen3}, {@code AT-F2-strassen}).
 * Compositions targeting {@code max-dim ≥ 16} are NOT exported because their
 * dense factor matrices grow to tens of MB per file.</p>
 */
public class ExportComposedSchemes {

	private static final File DIR = new File("src/main/resources/schemes");

	public static void main(String[] args) throws IOException {
		DIR.mkdirs();
		exportReal("derived_strassen3-8x8x8_r343.json", Compositions.strassen3_888());
		exportReal("derived_laderman_strassen-6x6x6_r161.json", Compositions.laderman_strassen_666());
		exportReal("derived_laderman2-9x9x9_r529.json", Compositions.laderman2_999());
		exportF2  ("derived_ATf2_strassen-8x8x8_r329.json", Compositions.alphatensorF2_strassen_888());
		// ⟨16,16,16⟩ via AlphaTensor² over F₂ = 47·47 = 2209. Matches fmm-lille's
		// best after composition; better than Strassen⁴=2401 over Z.
		exportF2  ("derived_ATf2_squared-16x16x16_r2209.json", Compositions.alphatensorF2_squared_16());
		// ⟨16,16,16⟩ over C via AlphaEvolve² = 48² = 2304. Matches fmm-lille's
		// listed Kronecker product. Dense complex; ~28 MB on disk.
		ComplexNonCubicBilinearAlgorithm ae = SchemeIO.readComplex(findAESection4());
		exportComplex("derived_AE2-16x16x16_r2304_0.5xC.json",
				Compose.kroneckerComplex(ae, ae), "0.5*C");
		System.out.printf("%nDone. Exported 6 composed schemes to %s%n", DIR.getAbsolutePath());
	}

	private static void exportComplex(String name, ComplexNonCubicBilinearAlgorithm alg,
			String fieldLabel) throws IOException {
		int maxDim = Math.max(alg.n, Math.max(alg.m, alg.p));
		File out = new File(DIR, "section" + maxDim + "/" + name);
		out.getParentFile().mkdirs();
		SchemeIO.write(alg, out, fieldLabel);
		int adds = Verifier.additionCount(alg);
		// Sampled verify (the full ⟨16,16,16⟩ complex residual is too slow here).
		int wrong = Verifier.residualSampledComplex(alg, 10_000, 0xCAFEBABE);
		System.out.printf("  %s  r=%d  +%d adds  complex-wrong=%d/10000%n",
				name, alg.r, adds, wrong);
	}

	private static void exportReal(String name, NonCubicBilinearAlgorithm alg) throws IOException {
		File out = sectionPathFor(alg, name);
		out.getParentFile().mkdirs();
		SchemeIO.write(alg, out);
		double residual = Verifier.residualNonCubic(alg);
		int adds = Verifier.additionCount(alg);
		System.out.printf("  %s  r=%d  +%d adds  residual=%.2e%n",
				name, alg.r, adds, residual);
	}

	private static void exportF2(String name, NonCubicBilinearAlgorithm alg) throws IOException {
		File out = sectionPathFor(alg, name);
		out.getParentFile().mkdirs();
		SchemeIO.write(alg, out);
		int wrong = Verifier.residualNonCubicF2(alg);
		int adds = Verifier.additionCount(alg);
		System.out.printf("  %s  r=%d  +%d adds  F2-wrong=%d%n",
				name, alg.r, adds, wrong);
	}

	private static File sectionPathFor(NonCubicBilinearAlgorithm alg, String name) {
		int maxDim = Math.max(alg.n, Math.max(alg.m, alg.p));
		return new File(DIR, "section" + maxDim + "/" + name);
	}

	private static File findAESection4() throws IOException {
		try (var s = java.nio.file.Files.walk(new File(DIR, "section4").toPath())) {
			return s.filter(p -> p.getFileName().toString().startsWith("alphaevolve-4x4x4_r48_"))
					.findFirst()
					.orElseThrow(() -> new IOException("AlphaEvolve 4x4x4=48 not in section4/"))
					.toFile();
		}
	}
}
