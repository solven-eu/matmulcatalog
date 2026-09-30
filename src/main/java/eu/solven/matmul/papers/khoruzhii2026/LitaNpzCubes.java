package eu.solven.matmul.papers.khoruzhii2026;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Comparator;
import java.util.Optional;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Stream;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.io.NpzSchemeReader;

/**
 * The Khoruzhii–Serafin–Gelß–Pokutta 2026 LITA cubes {@code ⟨N,N,N⟩}, 13 ≤ N ≤ 32, as
 * PUBLISHED (<a href="https://github.com/khoruzhii/lita">khoruzhii/lita</a>
 * {@code schemes/{N}x{N}x{N}_r{R}.npz}, repo state 2026-09-24 "LITA7 odd") — the
 * backing of the {@code TA_lita_npz(n=N)} lineage atom.
 *
 * <p><b>Why a data-backed atom.</b> The explicit factors are far too large for
 * catalog JSON (32³ is 117 MB canonical — over GitHub's 100 MB file limit; all 20
 * total 491 MB), while the upstream sparse-rational archives total 12 MB (each
 * ≤ 2.1 MB). A cube is therefore stored like every other big formula scheme: a
 * tiny hash-stamped stub whose lineage atom rebuilds the matrices on demand —
 * here by loading the archive through {@link NpzSchemeReader} rather than by
 * running a generator. {@link LitaTaConstruction} ({@code TA_lita(n)}) is our Java
 * port of their EARLIER (June 2026) Maple generator and is superseded at every N
 * by these ranks; porting {@code scripts/lita.py} / {@code lita_odd.py} so the atom
 * becomes a true generator ("derive, don't import") is the declared follow-up — it
 * must reproduce these exact matrices (same content hash), because projection
 * stubs pin their cube parent by {@code shape@hash}.</p>
 */
public final class LitaNpzCubes {
	private LitaNpzCubes() {}

	/** Resource directory holding the upstream archives (see its README for provenance). */
	public static final Path DIR = Path.of("src/main/resources/external/khoruzhii-lita");

	private static final Pattern NAME = Pattern.compile("^(\\d+)x(\\d+)x(\\d+)_r(\\d+)\\.npz$");

	/** The best (lowest-rank) archive for {@code ⟨n,n,n⟩}, if one is present. */
	public static Optional<Path> file(int n) {
		if (!Files.isDirectory(DIR)) return Optional.empty();
		try (Stream<Path> ls = Files.list(DIR)) {
			return ls.filter(p -> {
				Matcher m = NAME.matcher(p.getFileName().toString());
				return m.matches() && Integer.parseInt(m.group(1)) == n && Integer.parseInt(m.group(2)) == n
						&& Integer.parseInt(m.group(3)) == n;
			}).min(Comparator.comparingInt(LitaNpzCubes::rankOf));
		} catch (IOException e) {
			return Optional.empty();
		}
	}

	private static int rankOf(Path p) {
		Matcher m = NAME.matcher(p.getFileName().toString());
		return m.matches() ? Integer.parseInt(m.group(4)) : Integer.MAX_VALUE;
	}

	/** Load the published {@code ⟨n,n,n⟩} cube. Deterministic: same archive → same matrices. */
	public static NonCubicBilinearAlgorithm build(int n) {
		Path f = file(n).orElseThrow(() -> new IllegalStateException(
				"TA_lita_npz(n=" + n + "): no " + n + "x" + n + "x" + n + "_r*.npz under " + DIR));
		try {
			NonCubicBilinearAlgorithm alg = NpzSchemeReader.read(f.toFile()).alg();
			if (alg.n != n || alg.m != n || alg.p != n || alg.r != rankOf(f)) {
				throw new IllegalStateException("TA_lita_npz(n=" + n + "): " + f.getFileName() + " holds ⟨" + alg.n
						+ "," + alg.m + "," + alg.p + "⟩ r=" + alg.r + " — name/content mismatch");
			}
			return alg;
		} catch (IOException e) {
			throw new IllegalStateException("TA_lita_npz(n=" + n + "): cannot read " + f, e);
		}
	}
}
