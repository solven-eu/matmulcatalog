package eu.solven.matmul.docs.migrate;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.List;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import java.util.stream.Stream;

import eu.solven.matmul.catalog.SchemeIO;

import lombok.extern.slf4j.Slf4j;
import tools.jackson.databind.JsonNode;

/**
 * The idempotence index of the upstream scheme importers
 * ({@link ImportPerminovSchemes}, {@link ImportWitteveenSchemes}): which
 * {@code (shape, rank, coefficient class)} the catalog already holds under
 * {@code schemes/known/}, read from scheme <b>content</b> ({@code n}, {@code m},
 * {@code fields[]}, {@code commutative}, the stored {@code zt} flag or the
 * coefficients themselves) — never from the filename.
 *
 * <p>The class ladder is {@code ZT ⊂ Z ⊂ Q} (ternary integer ⊂ integer ⊂ rational):
 * a held ternary scheme covers an upstream {@code Z} or {@code Q} file at the same
 * {@code (shape, rank)}, but a held <em>non-ternary</em> scheme does <b>not</b> cover
 * an upstream {@code ZT} file. Until 2026-09-30 the Perminov importer keyed on
 * {@code (shape, rank)} alone (parsed from {@code -perminov_} filenames), so a new
 * ternary scheme at a rank we held only as {@code Z}/{@code Q} was skipped without
 * being downloaded — ⟨2,13,15⟩=304, ⟨2,13,16⟩=324, ⟨2,15,16⟩=374 and four more ZT
 * records never reached the catalog while the sync job reported "0 imported".</p>
 *
 * <p>Shapes are canonical (sorted dims): rank and ternarity are invariant under the
 * transposition symmetries, so a scheme held at any orientation covers the upstream
 * file. Commutative-only schemes and schemes not valid over {@code Q}
 * (F₂-native, C-only) cover nothing — every upstream catalog is NC, char-0.</p>
 */
@Slf4j
final class KnownSchemeKeys {
	private KnownSchemeKeys() {}

	static final String ZT = "ZT";
	static final String Z = "Z";
	static final String Q = "Q";

	/** {@code "{a}x{b}x{c}-r{rank}-{class}"} with {@code a≤b≤c}. */
	static String key(int n, int m, int p, int rank, String cls) {
		int[] d = { n, m, p };
		Arrays.sort(d);
		return d[0] + "x" + d[1] + "x" + d[2] + "-r" + rank + "-" + cls;
	}

	/** Do we hold a scheme at {@code (shape, rank)} whose class is {@code cls} or finer? */
	static boolean covers(Set<String> keys, int n, int m, int p, int rank, String cls) {
		return keys.contains(key(n, m, p, rank, cls));
	}

	/** Record a scheme of class {@code cls}: it also stands for every coarser class. */
	static void add(Set<String> keys, int n, int m, int p, int rank, String cls) {
		keys.add(key(n, m, p, rank, Q));
		if (Z.equals(cls) || ZT.equals(cls)) {
			keys.add(key(n, m, p, rank, Z));
		}
		if (ZT.equals(cls)) {
			keys.add(key(n, m, p, rank, ZT));
		}
	}

	/** Normalise an upstream filename tag / directory segment to a class. */
	static String classOfTag(String tag) {
		String t = tag.toUpperCase(java.util.Locale.ROOT);
		if (t.contains(ZT)) return ZT;
		if (t.contains(Q)) return Q;
		return Z;
	}

	/** Scan {@code known/} (parallel; every {@code *.json} is parsed once). */
	static Set<String> scan(Path known) throws IOException {
		Set<String> keys = ConcurrentHashMap.newKeySet();
		if (!Files.isDirectory(known)) {
			return keys;
		}
		List<Path> files;
		try (Stream<Path> w = Files.walk(known)) {
			files = w.filter(f -> f.getFileName().toString().endsWith(".json")).toList();
		}
		files.parallelStream().forEach(f -> {
			try {
				addFile(keys, f);
			} catch (Exception e) {
				log.debug("KnownSchemeKeys: skipping {}: {}", f.getFileName(), e.toString());
			}
		});
		return keys;
	}

	private static void addFile(Set<String> keys, Path f) throws IOException {
		JsonNode root = SchemeIO.parseJson(f.toFile());
		JsonNode dims = root.get("n");
		JsonNode rank = root.get("m");
		if (dims == null || !dims.isArray() || dims.size() != 3 || rank == null || !rank.isInt()) {
			return;
		}
		if (root.path("commutative").asBoolean(false)) {
			return;
		}
		List<String> fields = SchemeIO.fieldTags(root);
		if (!fields.isEmpty() && !fields.contains(Q)) {
			return; // F₂-native / C-only: no claim over the char-0 classes
		}
		String cls = Q;
		if (fields.contains(Z)) {
			cls = isTernary(root) ? ZT : Z;
		}
		add(keys, dims.get(0).asInt(), dims.get(1).asInt(), dims.get(2).asInt(), rank.asInt(), cls);
	}

	/** Stored {@code zt} flag when stamped, else the coefficients (explicit schemes only). */
	private static boolean isTernary(JsonNode root) throws IOException {
		Boolean stored = SchemeIO.readZT(root);
		if (stored != null) {
			return stored;
		}
		if (SchemeIO.isStub(root) || root.path("complex").asBoolean(false)) {
			return false;
		}
		return SchemeIO.isTernary(SchemeIO.read(root));
	}
}
