package eu.solven.matmul.search;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.EnumMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.concurrent.ConcurrentHashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Stream;

import tools.jackson.core.JsonParser;
import tools.jackson.core.JsonToken;
import tools.jackson.databind.JsonNode;
import tools.jackson.databind.node.ObjectNode;
import tools.jackson.databind.json.JsonMapper;

import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.Lineage;
import eu.solven.matmul.catalog.SchemeIO;
import lombok.extern.slf4j.Slf4j;

/**
 * The outer bases the catalog has already PROVEN useful: every scheme file whose lineage
 * carries a recombination node names the base it was built on, and that base must stay
 * in the default recombination pool whatever the pool's dimension cap says — otherwise
 * the closure sweep cannot re-derive what the catalog holds, and its "no improvement"
 * fixpoint is meaningless for those shapes.
 *
 * <p>Why this exists (2026-10-21): the default pool admits non-cubic bases up to
 * {@code maxBaseDim = 5}. 73 committed lineages stood on 9 dim-6/7 base files (⟨3,4,6⟩ 32,
 * ⟨2,5,6⟩ 17, ⟨3,4,7⟩ 12, ⟨3,5,6⟩, ⟨2,5,7⟩, ⟨2,6,6⟩, ⟨2,4,7⟩, ⟨5,5,7⟩, ⟨5,6,7⟩), every one of
 * them produced by a one-off {@code --base=} run or by an older sweep with a wider cap. Their
 * guards were disk-presence only, precisely because the default sweep could not re-derive
 * them, and the resync playbook's closure step had been declaring fixpoints over a pool
 * that could never reach those bases. The July 2026 follow-up ("a preset admitting a
 * curated dim-6/7 list") was never built; this is the evidence-driven version of it: the
 * list is read from the catalog, not curated by hand, and the guard
 * ({@code TestProvenOuterBases}) fails the moment a committed stub stands on a base the
 * default pool does not admit.</p>
 *
 * <p><b>Resolution.</b> A base reference is read through the replayer's single leaf grammar
 * ({@link LineageReplayer#classifyLeaf}): a pinned {@code shape@hash} names the file by
 * content; a legacy {@code ext[<canonical stem>.json] / <shape>=rank :: ORBIT} label names
 * it through its stem's hash7 (the file that produced the stub) AND through what the
 * replayer resolves it to today (the shape-best); a bare / legacy-stem ref resolves to the
 * best materialised non-commutative file at the shape under the file's own fields. Named
 * pool bases ({@code Strassen<2,2,2>=7}), elementary blocks ({@code naive-…}) and
 * parametric constructors are built in code and need no file.</p>
 *
 * <p><b>Cost.</b> The scan opens every JSON under the root but parses only the prefix up to
 * the {@code lineage} / {@code fields} keys (the canonical formatter writes the payload
 * last), so it costs about what the lookup index costs; the result is cached per root.</p>
 */
@Slf4j
public final class ProvenOuterBases {

	/** One proven base file: the hash prefix the catalog names it by, its (oriented) shape,
	 *  how many lineage files stand on it and one example. */
	public record Proven(String hashPrefix, int[] shape, int uses, Path example) {
		public String shapeKey() {
			return shape[0] + "x" + shape[1] + "x" + shape[2];
		}
		public int maxDim() {
			return Math.max(shape[0], Math.max(shape[1], shape[2]));
		}
	}

	/** A base reference that resolves to no file at all (dangling pin, unknown shape). */
	public record Unresolved(String ref, Path file) {}

	private static final Pattern PINNED = Pattern.compile("^(\\d+)x(\\d+)x(\\d+)@([0-9a-f]{7,})$");
	private static final Pattern EXT_LABEL = Pattern.compile("^ext\\[([^\\]]+?)\\]");

	private static final Map<Path, ProvenOuterBases> CACHE = new ConcurrentHashMap<>();
	private static final JsonMapper MAPPER = JsonMapper.builder().build();

	/** hash7 → every proven entry whose prefix starts with it (a pin may carry a longer prefix). */
	private final Map<String, List<Proven>> byHash7;
	private final List<Unresolved> unresolved;
	private final int builtInCode;
	private final int lineageFiles;

	private ProvenOuterBases(Map<String, List<Proven>> byHash7, List<Unresolved> unresolved, int builtInCode,
			int lineageFiles) {
		this.byHash7 = byHash7;
		this.unresolved = unresolved;
		this.builtInCode = builtInCode;
		this.lineageFiles = lineageFiles;
	}

	/** The proven bases of the schemes under {@code root}, scanned once per root (see {@link #invalidate}). */
	public static ProvenOuterBases of(Path root) {
		return CACHE.computeIfAbsent(root.toAbsolutePath().normalize(), ProvenOuterBases::scan);
	}

	/** Drop the cached scan of {@code root} (after writing stubs under it). */
	public static void invalidate(Path root) {
		CACHE.remove(root.toAbsolutePath().normalize());
	}

	/** Does the catalog name a file with this content hash as a recombination base? */
	public boolean admits(String contentHash) {
		if (contentHash == null || contentHash.length() < 7) return false;
		List<Proven> cands = byHash7.get(contentHash.substring(0, 7));
		if (cands == null) return false;
		for (Proven p : cands) {
			if (contentHash.startsWith(p.hashPrefix()) || p.hashPrefix().startsWith(contentHash)) return true;
		}
		return false;
	}

	/** Every proven base, most-used first. */
	public List<Proven> all() {
		List<Proven> out = new ArrayList<>();
		byHash7.values().forEach(out::addAll);
		out.sort((a, b) -> Integer.compare(b.uses(), a.uses()));
		return out;
	}

	/** The proven bases whose max dimension exceeds {@code maxBaseDim} (cubic ones: {@code maxCubicBaseDim}). */
	public List<Proven> aboveCap(int maxBaseDim, int maxCubicBaseDim) {
		List<Proven> out = new ArrayList<>();
		for (Proven p : all()) {
			boolean cubic = p.shape()[0] == p.shape()[1] && p.shape()[1] == p.shape()[2];
			int cap = cubic ? Math.max(maxBaseDim, maxCubicBaseDim) : maxBaseDim;
			if (p.maxDim() > cap) out.add(p);
		}
		return out;
	}

	/** The largest dimension among proven bases (0 when none) — the header cap a pool scan must read up to. */
	public int maxDim() {
		int m = 0;
		for (List<Proven> l : byHash7.values()) for (Proven p : l) m = Math.max(m, p.maxDim());
		return m;
	}

	public List<Unresolved> unresolved() {
		return unresolved;
	}

	/** Base references satisfied by code (named pool bases, naive blocks, parametric constructors). */
	public int builtInCode() {
		return builtInCode;
	}

	/** How many files under the root carried a lineage. */
	public int lineageFiles() {
		return lineageFiles;
	}

	// ── scan ─────────────────────────────────────────────────────────────────────────

	private record BaseRef(String ref, Path file, List<String> fieldTags) {}

	static ProvenOuterBases scan(Path root) {
		long t0 = System.nanoTime();
		List<Path> paths;
		try (Stream<Path> walk = Files.walk(root)) {
			paths = walk.filter(p -> p.getFileName().toString().endsWith(".json")).sorted().toList();
		} catch (IOException e) {
			throw new IllegalStateException("walking " + root, e);
		}
		List<BaseRef> refs = Collections.synchronizedList(new ArrayList<>());
		java.util.concurrent.atomic.AtomicInteger withLineage = new java.util.concurrent.atomic.AtomicInteger();
		paths.parallelStream().forEach(p -> {
			try {
				JsonNode head = readHead(p.toFile());
				if (head == null) return;
				Optional<Lineage.Node> lineage = SchemeIO.readLineage(head);
				if (lineage.isEmpty()) return;
				withLineage.incrementAndGet();
				List<String> tags = SchemeIO.fieldTags(head);
				List<Lineage.Node> bases = new ArrayList<>();
				collectRecombinationBases(lineage.get(), bases, Collections.newSetFromMap(new java.util.IdentityHashMap<>()));
				for (Lineage.Node b : bases) {
					Lineage.Node leaf = unwrap(b);
					if (leaf instanceof Lineage.Atom a) refs.add(new BaseRef(a.ref(), p, tags));
					// A composed base (an inline recombination / Kron as the outer base) is
					// not a pool file; its own recombination nodes were collected above.
				}
			} catch (Exception e) {
				log.debug("ProvenOuterBases: skipping {} ({})", p, e.toString());
			}
		});

		Resolver resolver = new Resolver(root);
		Map<String, Proven> byKey = new java.util.LinkedHashMap<>();   // hashPrefix+shape → entry
		List<Unresolved> unresolved = new ArrayList<>();
		int builtInCode = 0;
		for (BaseRef r : refs) {
			List<Target> targets = resolver.resolve(r);
			if (targets.isEmpty()) {
				unresolved.add(new Unresolved(r.ref(), r.file()));
				continue;
			}
			for (Target t : targets) {
				if (t.builtInCode()) {
					builtInCode++;
					continue;
				}
				String key = t.hashPrefix() + "|" + t.shape()[0] + "x" + t.shape()[1] + "x" + t.shape()[2];
				Proven cur = byKey.get(key);
				byKey.put(key, new Proven(t.hashPrefix(), t.shape(), cur == null ? 1 : cur.uses() + 1,
						cur == null ? r.file() : cur.example()));
			}
		}
		Map<String, List<Proven>> byHash7 = new java.util.HashMap<>();
		for (Proven p : byKey.values()) {
			byHash7.computeIfAbsent(p.hashPrefix().substring(0, 7), k -> new ArrayList<>()).add(p);
		}
		ProvenOuterBases out = new ProvenOuterBases(byHash7, unresolved, builtInCode, withLineage.get());
		log.info("ProvenOuterBases: {} lineage files under {} → {} base refs, {} proven base files ({} built in code,"
				+ " {} unresolved) in {} ms", withLineage.get(), root, refs.size(), byKey.size(), builtInCode,
				unresolved.size(), (System.nanoTime() - t0) / 1_000_000L);
		return out;
	}

	/** The scheme's base-reference-bearing nodes: every recombination node's {@code base}. */
	private static void collectRecombinationBases(Lineage.Node n, List<Lineage.Node> out,
			java.util.Set<Lineage.Node> seen) {
		if (n == null || !seen.add(n)) return;
		switch (n) {
			// An all-unit allocation over naive-1x1x1 leaves is the materialiser's spelling
			// of "this catalog scheme, as a leaf" (every product multiplies 1×1 blocks): the
			// base is a leaf reference, not a pool base anyone recombined on.
			case Lineage.RecombinationN x -> {
				if (!identity(x.allocA(), x.allocB(), x.allocC())) out.add(x.base());
			}
			case Lineage.RecombinationTaN x -> {
				if (!identity(x.allocA(), x.allocB(), x.allocC())) out.add(x.base());
			}
			case Lineage.RecombinationWithPairN x -> {
				if (!identity(x.allocA(), x.allocB(), x.allocC())) out.add(x.base());
			}
			default -> {
				// transforms, products, concats, projections: no base of their own
			}
		}
		for (Lineage.Node c : Lineage.childrenOf(n)) collectRecombinationBases(c, out, seen);
	}

	private static boolean identity(int[]... allocs) {
		for (int[] alloc : allocs) {
			if (alloc == null) return false;
			for (int a : alloc) if (a != 1) return false;
		}
		return true;
	}

	/** Strip the transform wrappers an outer base may carry (orientation only). */
	static Lineage.Node unwrap(Lineage.Node n) {
		while (true) {
			switch (n) {
				case Lineage.Transpose t -> n = t.child();
				case Lineage.OrientAs o -> n = o.child();
				case Lineage.AxisFlip f -> n = f.child();
				case Lineage.AxisPermute p -> n = p.child();
				case Lineage.Dce d -> n = d.child();
				default -> {
					return n;
				}
			}
		}
	}

	/**
	 * Read a scheme file's metadata prefix: {@code lineage} and {@code fields}, stopping at the
	 * first payload key so an explicit scheme's matrices are never tokenised. Returns
	 * {@code null} when the file carries no lineage before its payload.
	 */
	static JsonNode readHead(File f) throws IOException {
		ObjectNode head = MAPPER.createObjectNode();
		boolean lineage = false, fields = false;
		try (JsonParser p = MAPPER.createParser(f)) {
			if (p.nextToken() != JsonToken.START_OBJECT) return null;
			String name;
			while ((name = p.nextName()) != null) {
				p.nextToken();
				switch (name) {
					case "lineage" -> {
						head.set("lineage", p.readValueAsTree());
						lineage = true;
					}
					case "fields" -> {
						head.set("fields", p.readValueAsTree());
						fields = true;
					}
					case "u", "v", "w", "u_sparse", "v_sparse", "w_sparse", "U", "V", "W" -> {
						return lineage ? head : null;
					}
					default -> p.skipChildren();
				}
				if (lineage && fields) return head;
			}
		}
		return lineage ? head : null;
	}

	// ── resolution ───────────────────────────────────────────────────────────────────

	private record Target(String hashPrefix, int[] shape, boolean builtInCode) {
		static Target code() {
			return new Target(null, null, true);
		}
	}

	/** Resolves base refs under the file's own fields through the replayer's leaf grammar. */
	private static final class Resolver {
		private final Path root;
		private final Map<Field, LineageReplayer> replayers = new EnumMap<>(Field.class);
		private final Map<Field, FieldAwareLookup> lookups = new EnumMap<>(Field.class);
		/** (field, shape) → hash of the best materialised NC file — the shape-best resolution. */
		private final Map<String, Optional<String>> shapeBest = new java.util.HashMap<>();

		Resolver(Path root) {
			this.root = root;
		}

		List<Target> resolve(BaseRef r) {
			String ref = r.ref();
			List<Target> out = new ArrayList<>();
			// Legacy pool label "ext[<stem>.json] / <shape>=rank :: ORBIT": the stem names the
			// file that produced the stub (hash7 for canonical stems); keep that AND what the
			// ref resolves to today, since the replay follows the latter.
			Matcher ext = EXT_LABEL.matcher(ref);
			if (ext.find()) {
				String stemRef = Lineage.atomFromFilename(ext.group(1)).ref();
				Matcher pin = PINNED.matcher(stemRef);
				if (pin.matches()) out.add(pinned(pin));
			}
			Matcher pin = PINNED.matcher(ref);
			if (pin.matches()) {
				out.add(pinned(pin));
				return out;
			}
			for (Field f : fieldsToTry(r.fieldTags())) {
				LineageReplayer.LeafResolution lr = replayer(f).classifyLeaf(ref);
				switch (lr.cls()) {
					case PINNED -> out.add(new Target(SchemeIO.contentHash(lr.pinned().alg()), lr.shape(), false));
					case NAMED, NAIVE, PARAMETRIC -> out.add(Target.code());
					case SHAPE_BEST, LEGACY_SHAPE_BEST -> {
						int[] sh = lr.shape();
						String key = f + ":" + sh[0] + "x" + sh[1] + "x" + sh[2];
						Optional<String> hash = shapeBest.computeIfAbsent(key, k -> lookup(f)
								.findWithSource(sh[0], sh[1], sh[2])
								.map(ws -> SchemeIO.contentHash(ws.alg())));
						hash.ifPresent(h -> out.add(new Target(h, sh, false)));
					}
					case PINNED_STUB, DANGLING_PIN, UNRESOLVABLE -> {
						// not a pool file / nothing to admit
					}
				}
			}
			return out;
		}

		private static Target pinned(Matcher pin) {
			int[] sh = { Integer.parseInt(pin.group(1)), Integer.parseInt(pin.group(2)), Integer.parseInt(pin.group(3)) };
			return new Target(pin.group(4), sh, false);
		}

		private LineageReplayer replayer(Field f) {
			return replayers.computeIfAbsent(f, k -> LineageReplayer.withDefaultPool(lookup(k)));
		}

		private FieldAwareLookup lookup(Field f) {
			return lookups.computeIfAbsent(f, k -> new FieldAwareLookup(k, root));
		}
	}

	/**
	 * The lookups a file's leaves resolve under (mirrors {@code AuditLineageRefs.fieldsToTry}):
	 * C covers every characteristic-0 tag, F₂ and F₃ are their own universes.
	 */
	static List<Field> fieldsToTry(List<String> tags) {
		List<Field> out = new ArrayList<>();
		boolean char0 = tags.stream().anyMatch(t -> t.equals("Z") || t.equals("Q") || t.equals("R")
				|| t.equals("C") || t.equals("ZT"));
		if (char0 || tags.isEmpty()) out.add(Field.C);
		if (tags.contains("F2") || tags.contains("Z2") || tags.isEmpty()) out.add(Field.F2);
		if (tags.contains("F3") || tags.contains("Z3") || tags.isEmpty()) out.add(Field.F3);
		return out;
	}
}
