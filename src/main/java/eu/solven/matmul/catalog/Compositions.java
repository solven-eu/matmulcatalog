package eu.solven.matmul.catalog;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Optional;

import eu.solven.matmul.NonCubicBilinearAlgorithm;

/**
 * Convenience helpers for building large matmul algorithms by composing
 * smaller verified schemes from {@code src/main/resources/schemes/}.
 *
 * <p>For target format {@code ⟨a·b, c·d, e·f⟩} valid compositions are
 * {@code ⟨a, c, e⟩ ⊗ ⟨b, d, f⟩} (Kronecker product, see
 * {@link Compose#kroneckerGeneral}). The composed rank is the product
 * {@code r_outer · r_inner} and the resulting algorithm is correct over
 * any field where both factors are valid.</p>
 *
 * <p>For F₂ composition (e.g. AlphaTensor's `⟨4,4,4⟩=47` Kroneckered with
 * itself), the composed factors stay in {0, 1} and the trilinear identity
 * holds <b>modulo 2</b>. Verify via
 * {@link eu.solven.matmul.verifiers.Verifier#residualNonCubicF2}.</p>
 *
 * <p>For real-arithmetic composition (e.g. Strassen³), the composed factors
 * stay in {-1, 0, +1} and the identity holds exactly over Z.</p>
 */
public final class Compositions {

	private static final File SCHEMES = new File("src/main/resources/schemes");

	private Compositions() {}

	public static NonCubicBilinearAlgorithm loadScheme(String fileName) throws IOException {
		// Schemes live under section{N}/ subdirectories — walk recursively.
		//
		// Resolution is filename-tolerant by design: the 2026-06 catalog migration
		// made filenames pure labels in the content-driven form
		// `{n}x{m}x{p}-r{rank}-{note}-{hash7}.json`, so a legacy ref like
		// `strassen-2x2x2_m7` must NOT be pinned to a literal name (it would
		// re-break on the next rehash/rename). We match by the scheme's IDENTITY —
		// shape + rank + every author/field note token — extracted from the ref and
		// checked against the candidate stem, regardless of token order or hash.
		//   `strassen-2x2x2_m7`        → shape 2x2x2, rank 7,  notes [strassen]
		//   `alphatensor_F2-4x4x4_m47` → shape 4x4x4, rank 47, notes [alphatensor, F2]
		//   `laderman_1976-3x3x3_m23`  → shape 3x3x3, rank 23, notes [laderman, 1976]
		String base = fileName.endsWith(".json")
				? fileName.substring(0, fileName.length() - 5)
				: fileName;
		SchemeRef ref = SchemeRef.parse(base);
		File f;
		try (var s = Files.walk(SCHEMES.toPath())) {
			Optional<Path> hit = s.filter(p -> {
				String name = p.getFileName().toString();
				if (!name.endsWith(".json")) return false;
				String stem = name.substring(0, name.length() - 5);
				// Legacy exact / token-boundary prefix match (still works for any
				// file not yet renamed to the content-driven form), then the
				// content-identity match for migrated files.
				return stem.equals(base) || stem.startsWith(base + "_")
						|| (ref != null && ref.matches(stem));
			}).findFirst();
			if (hit.isEmpty()) {
				throw new IOException("scheme not found: " + fileName);
			}
			f = hit.get().toFile();
		}
		return SchemeIO.readBilinear(f);
	}

	/**
	 * The identity of a scheme as carried by both the legacy and content-driven
	 * filename conventions: its shape {@code n×m×p}, its rank, and the set of
	 * author/field note tokens. A {@code matches} is shape + rank + all notes
	 * present in the candidate stem — order- and hash-independent.
	 */
	private record SchemeRef(String shape, int rank, java.util.List<String> notes) {
		private static final java.util.regex.Pattern SHAPE =
				java.util.regex.Pattern.compile("(\\d+)x(\\d+)x(\\d+)");
		private static final java.util.regex.Pattern RANK =
				java.util.regex.Pattern.compile("[_-][mr](\\d+)\\b");

		static SchemeRef parse(String base) {
			java.util.regex.Matcher ms = SHAPE.matcher(base);
			java.util.regex.Matcher mr = RANK.matcher(base);
			if (!ms.find() || !mr.find()) return null;
			String shape = ms.group(1) + "x" + ms.group(2) + "x" + ms.group(3);
			int rank = Integer.parseInt(mr.group(1));
			// Note tokens = everything that is not the shape or the rank token.
			java.util.List<String> notes = new java.util.ArrayList<>();
			for (String tok : base.replace(ms.group(), " ").replaceAll("[_-][mr]\\d+\\b", " ").split("[_\\-]")) {
				if (!tok.isBlank()) notes.add(tok.toLowerCase(java.util.Locale.ROOT));
			}
			return new SchemeRef(shape, rank, notes);
		}

		boolean matches(String stem) {
			String low = stem.toLowerCase(java.util.Locale.ROOT);
			// Rank token bounded by delimiters so `r7` does not match `r70`.
			boolean rankHit = low.matches(".*[_-][mr]" + rank + "([_-].*|)$");
			if (!low.contains(shape) || !rankHit) return false;
			for (String note : notes) {
				if (!low.contains(note)) return false;
			}
			return true;
		}
	}

	/**
	 * Recursively compose {@code base} with itself {@code levels} times via
	 * Kronecker product. {@code levels=1} returns base; {@code levels=k}
	 * yields an algorithm for the {@code k}-times product format at rank
	 * {@code base.r^k}.
	 */
	public static NonCubicBilinearAlgorithm power(NonCubicBilinearAlgorithm base, int levels) {
		if (levels < 1) throw new IllegalArgumentException("levels >= 1");
		NonCubicBilinearAlgorithm acc = base;
		for (int i = 1; i < levels; i++) {
			acc = Compose.kroneckerGeneral(acc, base);
		}
		return acc;
	}

	// ───────────────────────────────────────────────────────────────────────────
	// Named recipes
	// ───────────────────────────────────────────────────────────────────────────

	/** {@code Strassen³ → ⟨8,8,8⟩ rank 343 over Z}. Three levels of Strassen. */
	public static NonCubicBilinearAlgorithm strassen3_888() throws IOException {
		return power(loadScheme("strassen-2x2x2_m7.json"), 3);
	}

	/** {@code Strassen⁴ → ⟨16,16,16⟩ rank 2,401 over Z}. */
	public static NonCubicBilinearAlgorithm strassen4_16() throws IOException {
		return power(loadScheme("strassen-2x2x2_m7.json"), 4);
	}

	/** {@code Strassen⁵ → ⟨32,32,32⟩ rank 16,807 over Z}. */
	public static NonCubicBilinearAlgorithm strassen5_32() throws IOException {
		return power(loadScheme("strassen-2x2x2_m7.json"), 5);
	}

	/**
	 * {@code AlphaTensor-F₂ ⟨4,4,4⟩=47 ⊗ Strassen ⟨2,2,2⟩=7 → ⟨8,8,8⟩
	 * rank 329 over F₂} — better than Strassen³=343 thanks to AlphaTensor's
	 * `R_{F₂}(⟨4,4,4⟩)` improvement.
	 */
	public static NonCubicBilinearAlgorithm alphatensorF2_strassen_888() throws IOException {
		NonCubicBilinearAlgorithm at = loadScheme("alphatensor_F2-4x4x4_m47.json");
		NonCubicBilinearAlgorithm st = loadScheme("strassen-2x2x2_m7.json");
		return Compose.kroneckerGeneral(at, st);
	}

	/**
	 * {@code AlphaTensor-F₂ ⟨4,4,4⟩=47 ⊗ AlphaTensor-F₂ ⟨4,4,4⟩=47 →
	 * ⟨16,16,16⟩ rank 2,209 over F₂} — strictly better than Strassen⁴=2,401.
	 */
	public static NonCubicBilinearAlgorithm alphatensorF2_squared_16() throws IOException {
		NonCubicBilinearAlgorithm at = loadScheme("alphatensor_F2-4x4x4_m47.json");
		return Compose.kroneckerGeneral(at, at);
	}

	/**
	 * {@code AlphaTensor-F₂ ⟨4,4,4⟩=47 ⊗ AlphaTensor-F₂ ⟨4,4,4⟩=47 ⊗ Strassen
	 * → ⟨32,32,32⟩ rank 15,463 over F₂}. The AlphaEvolve-style headline
	 * recursive composition.
	 */
	public static NonCubicBilinearAlgorithm alphatensorF2_squared_strassen_32() throws IOException {
		NonCubicBilinearAlgorithm at = loadScheme("alphatensor_F2-4x4x4_m47.json");
		NonCubicBilinearAlgorithm st = loadScheme("strassen-2x2x2_m7.json");
		return Compose.kroneckerGeneral(Compose.kroneckerGeneral(at, at), st);
	}

	/** {@code Laderman² → ⟨9,9,9⟩ rank 529 over Z}. */
	public static NonCubicBilinearAlgorithm laderman2_999() throws IOException {
		return power(loadScheme("laderman_1976-3x3x3_m23.json"), 2);
	}

	/**
	 * {@code Laderman ⊗ Strassen → ⟨6,6,6⟩ rank 161 over Z}. Sanity case —
	 * also verifiable.
	 */
	public static NonCubicBilinearAlgorithm laderman_strassen_666() throws IOException {
		return Compose.kroneckerGeneral(
				loadScheme("laderman_1976-3x3x3_m23.json"),
				loadScheme("strassen-2x2x2_m7.json"));
	}
}
