package eu.solven.matmul.catalog;

import eu.solven.matmul.recombination.Recombination;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.concurrent.ConcurrentHashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Stream;

import com.google.common.cache.Cache;
import com.google.common.cache.CacheBuilder;
import tools.jackson.databind.JsonNode;
import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.algebra.Field;

/**
 * Field-aware catalog lookup that implements the practical fallback chain
 * documented in {@code RANK_KNOWLEDGE.md §1.2bis "Practical fallback rule"}.
 *
 * <p>Downward transfer: any scheme valid over field {@code K} is automatically
 * valid over any super-field {@code K' ⊇ K}. So when a caller asks for a
 * scheme over {@code C} and none exists, an {@code R}-class scheme is a
 * mathematically valid substitute. {@code F_p} schemes have no fallback —
 * they live in their own characteristic and don't lift to / from
 * characteristic 0.</p>
 *
 * <p>Per-field chain (priority order):</p>
 * <ul>
 *   <li>{@code "Z"} → only {@code Z}</li>
 *   <li>{@code "Q"} → {@code Q}, {@code Z}</li>
 *   <li>{@code "R"} → {@code R}, {@code Q}, {@code Z}</li>
 *   <li>{@code "C"} → {@code C}, {@code R}, {@code Q}, {@code Z}</li>
 *   <li>{@code "F2"} → only {@code F_2} / {@code Z_2}</li>
 * </ul>
 *
 * <p>Returns the lowest-rank scheme available in the chain (so a richer
 * field's catalog-entry is preferred when it improves rank).</p>
 */
public final class FieldAwareLookup implements Recombination.AlgorithmLookup {

	private final eu.solven.matmul.algebra.Field field;
	private final Path schemesRoot;

	/** Type-safe constructor. */
	public FieldAwareLookup(eu.solven.matmul.algebra.Field field) {
		this(field, Path.of("src/main/resources/schemes"));
	}

	public FieldAwareLookup(eu.solven.matmul.algebra.Field field, Path schemesRoot) {
		this.field = field;
		this.schemesRoot = schemesRoot;
	}

	/** The field this lookup resolves over. */
	public eu.solven.matmul.algebra.Field field() {
		return field;
	}

	/** @deprecated use {@link #FieldAwareLookup(eu.solven.matmul.algebra.Field)} */
	@Deprecated
	public FieldAwareLookup(String field) {
		this(eu.solven.matmul.algebra.Field.fromTag(field));
	}

	/** @deprecated use {@link #FieldAwareLookup(eu.solven.matmul.algebra.Field, Path)} */
	@Deprecated
	public FieldAwareLookup(String field, Path schemesRoot) {
		this(eu.solven.matmul.algebra.Field.fromTag(field), schemesRoot);
	}

	/**
	 * Classifies a scheme file by its name into a {@link eu.solven.matmul.algebra.Field}.
	 * Filename tokens (case-insensitive for F₂):
	 * {@code F2}/{@code f2}/{@code ATf2}/{@code Z2} → {@code F2};
	 * {@code _C} or {@code xC} → {@code C}; else defaults to {@code Z_Q_R}.
	 */
	/**
	 * Field for a catalog file, preferring the JSON's unified {@code fields[]}
	 * (task #174) over the filename. Since #173 dropped trailing field tokens
	 * from filenames, the JSON body is now the authoritative source; the
	 * filename classifier is a fallback for un-migrated / synthetic names.
	 *
	 * <p>From {@code fields[]} we pick the STRICTEST characteristic-0 tag
	 * (Z ⊂ Q ⊂ R ⊂ C) so an integer scheme indexes as Z, not C. A
	 * prime-field-only scheme indexes as F2 / F3.</p>
	 */
	static eu.solven.matmul.algebra.Field classifyField(Path path, String name) {
		// Fast path: the narrowest field is already computed per scheme in the
		// docs/catalog.json manifest. Reading it once (one 7 MB metadata file)
		// avoids a full SchemeIO.parseJson of EVERY scheme JSON during buildIndex
		// — those parses (incl. the big matrix payloads of maxDim≥16 schemes) were
		// the dominant cost of the first lookup. Files absent from the manifest
		// (e.g. fresh staging writes) fall through to a per-file parse below.
		eu.solven.matmul.algebra.Field fromManifest = manifestFields().get(name);
		if (fromManifest != null) {
			return fromManifest;
		}
		try {
			java.util.List<String> tags = SchemeIO.fieldTags(SchemeIO.parseJson(path.toFile()));
			for (String t : new String[] { "Z", "Q", "R", "C", "F3", "F2" }) {
				if (tags.contains(t)) return eu.solven.matmul.algebra.Field.fromTag(t);
			}
		} catch (Exception e) {
			// fall through to the integer-by-default classification
		}
		// Content unreadable and absent from the manifest — default to Z
		// (the catalog is integer-coefficient by default). Never parse the filename.
		return eu.solven.matmul.algebra.Field.Z;
	}

	/** Per-file narrowest field, loaded ONCE from {@code docs/catalog.json}. Keyed
	 *  by the scheme file's basename. Empty if the manifest is missing/unreadable,
	 *  in which case {@link #classifyField} falls back to per-file parsing. */
	private static volatile Map<String, eu.solven.matmul.algebra.Field> MANIFEST_FIELDS;

	/** Full manifest entry per scheme file (keyed by the entry's relative {@code file}),
	 *  loaded ONCE from {@code docs/catalog.json}. Lets {@link #buildIndex} read
	 *  shape/rank/fields/commutative WITHOUT parsing the ~13k schemes the manifest
	 *  already describes — only disk files ABSENT from the manifest are parsed. */
	private static volatile Map<String, JsonNode> MANIFEST_BY_FILE;

	/** @return relative-{@code file} → manifest entry, or empty if the manifest is
	 *  missing/unreadable (then {@link #buildIndex} parses every file, as before). */
	private static Map<String, JsonNode> manifestByFile() {
		Map<String, JsonNode> m = MANIFEST_BY_FILE;
		if (m != null) return m;
		synchronized (FieldAwareLookup.class) {
			if (MANIFEST_BY_FILE != null) return MANIFEST_BY_FILE;
			Map<String, JsonNode> loaded = new HashMap<>();
			Path manifest = Path.of("docs", "catalog.json");
			if (Files.isReadable(manifest)) {
				try {
					JsonNode schemes = SchemeIO.parseJson(manifest.toFile()).get("schemes");
					if (schemes != null && schemes.isArray()) {
						for (JsonNode e : schemes) {
							JsonNode f = e.get("file");
							if (f != null && f.isTextual()) loaded.put(f.asString(), e);
						}
					}
				} catch (Exception ex) {
					LOG.warn("manifestByFile: docs/catalog.json unreadable ({}) — buildIndex will"
							+ " parse every scheme file", ex.toString());
				}
			}
			MANIFEST_BY_FILE = loaded;
			return loaded;
		}
	}

	private static Map<String, eu.solven.matmul.algebra.Field> manifestFields() {
		Map<String, eu.solven.matmul.algebra.Field> m = MANIFEST_FIELDS;
		if (m != null) return m;
		synchronized (FieldAwareLookup.class) {
			if (MANIFEST_FIELDS != null) return MANIFEST_FIELDS;
			Map<String, eu.solven.matmul.algebra.Field> loaded = new HashMap<>();
			Path manifest = Path.of("docs", "catalog.json");
			if (Files.isReadable(manifest)) {
				try {
					JsonNode root = SchemeIO.parseJson(manifest.toFile());
					JsonNode schemes = root.get("schemes");
					if (schemes != null && schemes.isArray()) {
						for (JsonNode e : schemes) {
							JsonNode fileNode = e.get("file");
							JsonNode fieldsNode = e.get("fields");
							if (fileNode == null || fieldsNode == null || !fieldsNode.isArray()) continue;
							String base = fileNode.asString();
							int slash = base.lastIndexOf('/');
							if (slash >= 0) base = base.substring(slash + 1);
							eu.solven.matmul.algebra.Field narrow = narrowestField(fieldsNode);
							if (narrow != null) loaded.put(base, narrow);
						}
					}
				} catch (Exception ex) {
					// Manifest unreadable/partial — keep whatever loaded; the
					// per-file fallback in classifyField preserves correctness.
				}
			}
			MANIFEST_FIELDS = loaded;
			return loaded;
		}
	}

	/** Per-sorted-shape max projection margin μ, loaded ONCE from {@code docs/catalog.json}.
	 *  Keyed {@code "a x b x c"} (a≤b≤c — μ = max over axes is permutation-invariant, so the
	 *  sorted key is exact). Empty if the manifest is missing/unreadable, in which case the
	 *  projection prune is simply disabled (callers fall back to replaying the parent). */
	private static volatile Map<String, Integer> MANIFEST_MARGIN;

	private static Map<String, Integer> manifestMargin() {
		Map<String, Integer> m = MANIFEST_MARGIN;
		if (m != null) return m;
		synchronized (FieldAwareLookup.class) {
			if (MANIFEST_MARGIN != null) return MANIFEST_MARGIN;
			Map<String, Integer> loaded = new HashMap<>();
			Path manifest = Path.of("docs", "catalog.json");
			if (Files.isReadable(manifest)) {
				try {
					JsonNode root = SchemeIO.parseJson(manifest.toFile());
					JsonNode schemes = root.get("schemes");
					if (schemes != null && schemes.isArray()) {
						for (JsonNode e : schemes) {
							JsonNode fmt = e.get("format");
							JsonNode pm = e.get("projection_margin");
							if (fmt == null || !fmt.isArray() || fmt.size() != 3
									|| pm == null || !pm.isIntegralNumber()) continue;
							int[] d = { fmt.get(0).asInt(), fmt.get(1).asInt(), fmt.get(2).asInt() };
							Arrays.sort(d);
							loaded.merge(d[0] + "x" + d[1] + "x" + d[2], pm.asInt(), Math::max);
						}
					}
				} catch (Exception ex) {
					// Manifest unreadable/partial — leave empty; the projection prune
					// disables itself per-shape (margin -1 → caller replays), so a
					// missing manifest only costs speed, never correctness.
				}
			}
			MANIFEST_MARGIN = loaded;
			return loaded;
		}
	}

	/**
	 * Upper bound on the projection margin μ over all catalog schemes for shape
	 * {@code ⟨n,m,p⟩} (max over that shape's schemes), read once from
	 * {@code docs/catalog.json}. Drives the downward-projection prune: a parent of
	 * rank {@code R}, projected by dropping one index on each of {@code k} increased
	 * axes, has rank {@code ≥ R − k·μ} — every dropped index DCEs at most μ products,
	 * a rigorous bound for the {@code ≤1-row-per-axis} drops the projector performs.
	 * If even {@code R − k·μ ≥} the incumbent, no projection of any same-shape parent
	 * can win, so the (expensive stub) replay is skipped.
	 *
	 * <p>Returns {@code -1} when the shape carries no manifest margin — the caller
	 * MUST then NOT prune (fall back to replay), so an un-enriched parent is never
	 * dropped on missing data.</p>
	 */
	public int projectionMarginUpperBound(int n, int m, int p) {
		int[] d = { n, m, p };
		Arrays.sort(d);
		Integer v = manifestMargin().get(d[0] + "x" + d[1] + "x" + d[2]);
		return v == null ? -1 : v;
	}

	/** Basenames flagged {@code corrupted:true} in the manifest — lineage stubs whose
	 *  recipe can no longer be replayed. Loaded ONCE; empty if manifest missing. */
	private static volatile java.util.Set<String> MANIFEST_CORRUPTED;

	/**
	 * The set of scheme basenames the manifest marks {@code corrupted:true}. Such a
	 * scheme carries a <em>claimed</em> rank that is NOT currently reproducible, so
	 * it must be treated as ABSENT by this gating lookup — it cannot shadow a fresh,
	 * verifiable discovery (the dual of the phantom-win bug). It still appears in the
	 * catalog/SPA, flagged. Falls back to empty (include everything) if the manifest
	 * is unreadable, so a missing manifest never hides real schemes.
	 */
	private static java.util.Set<String> manifestCorrupted() {
		java.util.Set<String> s = MANIFEST_CORRUPTED;
		if (s != null) return s;
		synchronized (FieldAwareLookup.class) {
			if (MANIFEST_CORRUPTED != null) return MANIFEST_CORRUPTED;
			java.util.Set<String> loaded = new java.util.HashSet<>();
			Path manifest = Path.of("docs", "catalog.json");
			if (Files.isReadable(manifest)) {
				try {
					JsonNode root = SchemeIO.parseJson(manifest.toFile());
					// Preferred: the top-level corrupted_files array (the COMPLETE set,
					// including same-shape siblings dropped by the manifest's shape-dedup).
					JsonNode corr = root.get("corrupted_files");
					if (corr != null && corr.isArray()) {
						for (JsonNode b : corr) {
							if (b.isTextual()) loaded.add(b.asString());
						}
					}
					// Fallback for older manifests without that array: scan schemes[].
					JsonNode schemes = root.get("schemes");
					if (loaded.isEmpty() && schemes != null && schemes.isArray()) {
						for (JsonNode e : schemes) {
							JsonNode corrupt = e.get("corrupted");
							JsonNode fileNode = e.get("file");
							if (fileNode == null || corrupt == null || !corrupt.asBoolean(false)) continue;
							String base = fileNode.asString();
							int slash = base.lastIndexOf('/');
							if (slash >= 0) base = base.substring(slash + 1);
							loaded.add(base);
						}
					}
				} catch (Exception ex) {
					// Manifest unreadable — fall back to empty (include everything).
				}
			}
			MANIFEST_CORRUPTED = loaded;
			return loaded;
		}
	}

	/** Narrowest field tag (Z &lt; Q &lt; R &lt; C, then F3, F2) in a manifest
	 *  {@code fields} array, or null if none recognised. */
	private static eu.solven.matmul.algebra.Field narrowestField(JsonNode fieldsArr) {
		java.util.Set<String> tags = new java.util.HashSet<>();
		for (JsonNode t : fieldsArr) {
			if (t.isTextual()) tags.add(t.asString());
		}
		for (String t : new String[] { "Z", "Q", "R", "C", "F3", "F2" }) {
			if (tags.contains(t)) return eu.solven.matmul.algebra.Field.fromTag(t);
		}
		return null;
	}

	/**
	 * Field-inference result from a {@link Lineage.Node} tree. {@link #field}
	 * is the narrowest characteristic-0 field (Z &lt; Q &lt; R &lt; C) all leaves
	 * agree on, OR a finite field (F2, F3) if every leaf is in that
	 * characteristic, OR {@link Optional#empty()} when we couldn't decide
	 * (mixed characteristics, or a leaf with an unrecognised field tag).
	 *
	 * <p>{@link #unknownLeaves} carries leaves whose ref didn't resolve to
	 * a recognised filename token — used by callers to log "audit me".</p>
	 */
	public record InferredField(Optional<Field> field, List<String> unknownLeaves) {}

	/**
	 * Canonical narrowest-field → inclusion-correct {@code fields[]} expansion,
	 * ordered {@code [F2,F3,Z,Q,R,C]}. The SINGLE source of truth shared by
	 * {@code StampFields} (lineage stub-stamper), the materialiser's born-stamped
	 * stub write ({@code RecursiveMaterialiser.writeToDisk}), and any other caller
	 * — so the expansion can never drift between them (these stampers historically
	 * disagreed; this method exists to prevent a recurrence).
	 *
	 * <p>F₂/F₃ are separate prime fields and never lift to characteristic 0; a Z
	 * (integer-exact) scheme reduces mod ANY prime, so Z grants F₂ and F₃ as a
	 * theorem. The expansion is exactly as trustworthy as the input field claim —
	 * a scheme that over-claims Z while holding rational coefficients would emit a
	 * wrong set, so callers must infer the field faithfully (content for matrices,
	 * lineage for stubs).</p>
	 */
	public static List<String> inclusionFieldNames(Field f) {
		return switch (f) {
			case Z -> List.of("F2", "F3", "Z", "Q", "R", "C");
			case Q -> List.of("Q", "R", "C");
			case R -> List.of("R", "C");
			case C -> List.of("C");
			case F2 -> List.of("F2");
			case F3 -> List.of("F3");
		};
	}

	/**
	 * Content-only {@code fields[]} for a composed scheme = the set-INTERSECTION
	 * of its leaf atoms' stamped {@code fields[]} (a composition is valid over a
	 * field iff EVERY leaf is). Each atom is resolved by SHAPE — via the
	 * content-driven index, NEVER the filename (CLAUDE.md: no filename-derived
	 * metadata) — to the lowest-rank non-commutative entry in this lookup's field
	 * chain (the same scheme a build resolves) and its {@code fields[]} read from
	 * CONTENT. Lineage-only (no matrices), so valid for matrix-less stubs.
	 *
	 * <p>Returns an EMPTY list when any leaf is unresolvable / unstamped — the
	 * caller MUST then refuse to stamp rather than over-claim. This replaces the
	 * single-{@link Field} inference, which (a) collapsed sets like {@code
	 * [F3,Q,R,C]} to one field and (b) fell back to {@code [Z]} on failure, both
	 * producing Z over-claims on rational-leaf schemes (e.g. ⟨5,32,32⟩ via the
	 * rational ⟨3,8,8⟩=145 leaf must be {@code [F3,Q,R,C]}, not {@code [Z,…]}).</p>
	 */
	public List<String> fieldNamesFromLineage(Lineage.Node lineage) {
		List<List<String>> leafFields = new ArrayList<>();
		collectLeafFields(lineage, leafFields);
		if (leafFields.isEmpty()) {
			return List.of();
		}
		java.util.LinkedHashSet<String> acc = null;
		for (List<String> leaf : leafFields) {
			if (leaf.isEmpty()) {
				// FLOOR (user 2026-06-13: "always at least shrink into the requested
				// field"): a build over `field` resolved this leaf to a field-valid
				// scheme, so the leaf is ≥ that field — fall back to the field's
				// inclusion chain. Safe: never claims F2/F3/Z from an unresolved leaf
				// (so no over-claim), only possibly under-claims. Recovers the win
				// instead of dropping it; still no over-claim.
				leaf = inclusionFieldNames(field);
			}
			if (acc == null) {
				acc = new java.util.LinkedHashSet<>(leaf);
			} else {
				acc.retainAll(leaf);
			}
		}
		List<String> out = new ArrayList<>();
		for (String t : List.of("F2", "F3", "Z", "Q", "R", "C")) {
			if (acc.contains(t)) {
				out.add(t);
			}
		}
		return out;
	}

	private static final Pattern SHAPE_IN_REF = Pattern.compile("(\\d+)x(\\d+)x(\\d+)");
	private static final Pattern ANGLE_SHAPE = Pattern.compile("<(\\d+),(\\d+),(\\d+)>");

	/** Field-set contributed by one lineage leaf:
	 *  <ul>
	 *    <li>a SHAPED atom (ref carries {@code NxMxP} / {@code <n,m,p>}) → its content
	 *        {@code fields[]} via {@link #bestFieldsAtShape} (or empty → caller floors);</li>
	 *    <li>a trivial dim-1 leaf → integer (full chain);</li>
	 *    <li>a parametric FORMULA atom (ref has no shape, e.g. {@code DIS09Lemma4(n=20)})
	 *        → the construction's INTRINSIC field — Pan/Islam trilinear-aggregation cube
	 *        is Q-strict (divides by {@code n+1}; lifts to R/C; NOT valid over F₂/Z), so
	 *        {@code [Q,R,C]} — NOT whatever ⟨n,n,n⟩ scheme the catalog happens to carry
	 *        (which would over-claim integer for a rational leaf, and is why a shape-based
	 *        lookup is wrong here). Without this, the formula leaf contributed no shape,
	 *        {@code leafFields} was empty, and the write-time field stamp failed loud.</li>
	 *  </ul>
	 *  Empty = unresolved → caller floors to the field's inclusion chain. */
	private List<String> atomFields(String ref) {
		int[] s = shapeOfRef(ref);
		if (s != null) {
			return shapeFields(s[0], s[1], s[2]);
		}
		if (ref.startsWith("DIS09Lemma4")) {
			return List.of("Q", "R", "C"); // Pan/Islam TA cube — Q-strict
		}
		if (ref.contains("naive")) {
			return List.of("F2", "F3", "Z", "Q", "R", "C"); // elementary products — integer
		}
		return List.of(); // unrecognised formula → floor
	}

	/** A trivial leaf (some dim = 1) is a naive outer/inner product — always
	 *  integer — so it carries the full chain without a catalog entry; otherwise
	 *  resolve by shape (content) to the build's scheme. */
	private List<String> shapeFields(int n, int m, int p) {
		return (n == 1 || m == 1 || p == 1)
				? List.of("F2", "F3", "Z", "Q", "R", "C")
				: bestFieldsAtShape(n, m, p);
	}

	/** Collect every leaf's field-set in the lineage (mirrors {@code Lineage.walk}). */
	private void collectLeafFields(Lineage.Node n, List<List<String>> out) {
		switch (n) {
			case Lineage.Atom a -> out.add(atomFields(a.ref()));
			case Lineage.KronProduct kp -> { collectLeafFields(kp.outer(), out); collectLeafFields(kp.inner(), out); }
			case Lineage.KronChain kc -> { for (Lineage.Node f : kc.factors()) collectLeafFields(f, out); }
			case Lineage.ConcatCols c -> { collectLeafFields(c.left(), out); collectLeafFields(c.right(), out); }
			case Lineage.ConcatRows c -> { collectLeafFields(c.top(), out); collectLeafFields(c.bottom(), out); }
			case Lineage.SumInner si -> { collectLeafFields(si.left(), out); collectLeafFields(si.right(), out); }
			case Lineage.RecombinationN r -> { collectLeafFields(r.base(), out); for (Lineage.Node lf : r.leaves()) collectLeafFields(lf, out); }
			case Lineage.RecombinationTaN r -> { collectLeafFields(r.base(), out); for (Lineage.Node lf : r.leaves()) collectLeafFields(lf, out); }
			case Lineage.RecombinationWithPairN r -> { collectLeafFields(r.base(), out); for (Lineage.Node lf : r.leaves()) collectLeafFields(lf, out); }
			case Lineage.AugmentSquareDiscard a -> collectLeafFields(a.square(), out);
			case Lineage.Dce d -> collectLeafFields(d.child(), out);
			case Lineage.Transpose t -> collectLeafFields(t.child(), out);
			case Lineage.OrientAs o -> collectLeafFields(o.child(), out);
			case Lineage.AxisFlip af -> collectLeafFields(af.child(), out);
			case Lineage.AxisPermute ap -> collectLeafFields(ap.child(), out);
			case Lineage.DisjointSum ds -> { for (Lineage.Node c : ds.children()) collectLeafFields(c, out); }
			case Lineage.SerendipitousProduct sp -> { collectLeafFields(sp.base(), out); out.add(shapeFields(sp.n2(), sp.m2(), sp.p2())); }
			case Lineage.Project pr -> collectLeafFields(pr.child(), out);
			// PeeledViaTa: field-relevant leaves are the cube + corner (TA part is integer ±1).
			case Lineage.PeeledViaTa t -> { collectLeafFields(t.cube(), out); collectLeafFields(t.corner(), out); }
		}
	}

	private static int[] shapeOfRef(String ref) {
		Matcher m = SHAPE_IN_REF.matcher(ref);
		if (m.find()) {
			return new int[] { Integer.parseInt(m.group(1)), Integer.parseInt(m.group(2)), Integer.parseInt(m.group(3)) };
		}
		Matcher a = ANGLE_SHAPE.matcher(ref);
		if (a.find()) {
			return new int[] { Integer.parseInt(a.group(1)), Integer.parseInt(a.group(2)), Integer.parseInt(a.group(3)) };
		}
		return null;  // unrecognised ref → unresolvable
	}

	/** Content {@code fields[]} of the lowest-rank non-commutative entry at
	 *  {@code ⟨n,m,p⟩} within this lookup's field chain (the scheme a build
	 *  resolves); empty if none. Reads the chosen file's {@code fields[]}. */
	List<String> bestFieldsAtShape(int n, int m, int p) {
		int[] s = { n, m, p };
		Arrays.sort(s);
		List<FileEntry> bucket = fileIndex().get(s[0] + "x" + s[1] + "x" + s[2]);
		if (bucket == null) {
			return List.of();
		}
		List<eu.solven.matmul.algebra.Field> chain = field.fallbackChain();
		FileEntry best = null;
		for (FileEntry e : bucket) {
			if (!chain.contains(e.field) || e.commutative) {
				continue;
			}
			if (best == null || e.rank < best.rank) {
				best = e;
			}
		}
		if (best == null) {
			return List.of();
		}
		try {
			return SchemeIO.fieldTags(SchemeIO.parseJson(best.path.toFile()));
		} catch (Exception ex) {
			return List.of();
		}
	}

	/**
	 * Walk a {@link Lineage.Node} tree, classify each {@link Lineage.Atom}
	 * by its on-disk filename, and return the narrowest field the
	 * composition is valid over (i.e. the intersection of the leaves'
	 * "valid-over" sets).
	 *
	 * <p>Combination semantics — all composition primitives
	 * (KronProduct, KronChain, ConcatCols, ConcatRows, SumInner, RecombinationN,
	 * RecombinationWithPairN, AugmentSquareDiscard, Dce, Transpose) are
	 * "field-preserving": the composed scheme's coefficients are
	 * polynomials in the leaves' coefficients, so the result lives in the
	 * smallest field containing all leaves. With our containment lattice
	 * Z ⊂ Q ⊂ R ⊂ C, "smallest containing" is {@code MAX} in that order.
	 * F2 / F3 don't share a containment relation with characteristic-0
	 * fields; mixing them yields {@link Optional#empty()} (= "unknown").
	 * Per project convention: if any leaf is F2-only, the composition
	 * is F2-only (we assume the char-0 leaves were Z-coefficient and
	 * reducible mod 2).</p>
	 *
	 * <p>Leaf resolution:</p>
	 * <ul>
	 *   <li>If the leaf ref ends in "-direct" / "-naive" or is a bare
	 *       shape (matched by {@link LineageReplayer}-style pattern), we
	 *       resolve the ref to an on-disk file via {@code lookupForFile}
	 *       (by shape) and recurse into that file's own lineage if it has
	 *       one; otherwise we classify the file's name.</li>
	 *   <li>If the leaf ref matches a known filename stem (we add ".json"
	 *       and search for the file), we classify by filename + recurse
	 *       into that file's lineage if present.</li>
	 *   <li>Named-base refs (e.g. {@code Strassen<2,2,2>=7}) are assumed
	 *       Z-coefficient (these are the canonical hand-coded
	 *       integer-coefficient bases).</li>
	 *   <li>Otherwise the leaf is added to {@link InferredField#unknownLeaves}
	 *       and treated conservatively (returns "unknown").</li>
	 * </ul>
	 */
	public InferredField inferFieldFromLineage(Lineage.Node root) {
		List<String> unknown = new ArrayList<>();
		java.util.Set<Path> visited = new java.util.HashSet<>();
		Optional<Field> f = inferFieldRecursive(root, unknown, visited);
		return new InferredField(f, unknown);
	}

	/** {@code Optional.empty()} means "unknown / can't decide". */
	private Optional<Field> inferFieldRecursive(Lineage.Node n, List<String> unknown,
			java.util.Set<Path> visited) {
		return switch (n) {
			case Lineage.Atom l -> classifyLeaf(l.ref(), unknown, visited);
			case Lineage.KronProduct kp -> combine(
					inferFieldRecursive(kp.outer(), unknown, visited),
					inferFieldRecursive(kp.inner(), unknown, visited));
			case Lineage.KronChain kc -> {
				Optional<Field> acc = Optional.of(Field.Z);
				for (Lineage.Node f : kc.factors()) {
					acc = combine(acc, inferFieldRecursive(f, unknown, visited));
				}
				yield acc;
			}
			case Lineage.ConcatCols cr -> combine(
					inferFieldRecursive(cr.left(), unknown, visited),
					inferFieldRecursive(cr.right(), unknown, visited));
			case Lineage.ConcatRows cb -> combine(
					inferFieldRecursive(cb.top(), unknown, visited),
					inferFieldRecursive(cb.bottom(), unknown, visited));
			case Lineage.SumInner si -> combine(
					inferFieldRecursive(si.left(), unknown, visited),
					inferFieldRecursive(si.right(), unknown, visited));
			case Lineage.RecombinationN r -> {
				Optional<Field> acc = inferFieldRecursive(r.base(), unknown, visited);
				for (Lineage.Node lf : r.leaves()) {
					acc = combine(acc, inferFieldRecursive(lf, unknown, visited));
				}
				yield acc;
			}
			case Lineage.RecombinationTaN r -> {
				// TA block is integer ±1 (field-neutral); field = base ∩ unpaired leaves.
				Optional<Field> acc = inferFieldRecursive(r.base(), unknown, visited);
				for (Lineage.Node lf : r.leaves()) {
					acc = combine(acc, inferFieldRecursive(lf, unknown, visited));
				}
				yield acc;
			}
			case Lineage.RecombinationWithPairN r -> {
				Optional<Field> acc = inferFieldRecursive(r.base(), unknown, visited);
				for (Lineage.Node lf : r.leaves()) {
					acc = combine(acc, inferFieldRecursive(lf, unknown, visited));
				}
				yield acc;
			}
			case Lineage.AugmentSquareDiscard a -> inferFieldRecursive(a.square(), unknown, visited);
			case Lineage.Dce d -> inferFieldRecursive(d.child(), unknown, visited);
			case Lineage.Transpose t -> inferFieldRecursive(t.child(), unknown, visited);
			case Lineage.OrientAs o -> inferFieldRecursive(o.child(), unknown, visited);
			// Axis-flip / axis-permute are symmetry rewrites of an underlying scheme —
			// they only rearrange (U, V, W) entries, so they preserve the field exactly.
			case Lineage.AxisFlip af -> inferFieldRecursive(af.child(), unknown, visited);
			case Lineage.AxisPermute ap -> inferFieldRecursive(ap.child(), unknown, visited);
			// DisjointSum: the composed field is the union (max-field) over all children.
			case Lineage.DisjointSum ds -> {
				Optional<eu.solven.matmul.algebra.Field> acc = Optional.of(eu.solven.matmul.algebra.Field.Z);
				for (Lineage.Node c : ds.children()) {
					acc = combine(acc, inferFieldRecursive(c, unknown, visited));
				}
				yield acc;
			}
			// Serendipitous product: field = base ⊕ second-scheme family.
			case Lineage.SerendipitousProduct sp -> combine(
					inferFieldRecursive(sp.base(), unknown, visited),
					classifyLeaf(sp.n2() + "x" + sp.m2() + "x" + sp.p2(), unknown, visited));
			// Projection: restricting indices + DCE preserves the field of the parent.
			case Lineage.Project pr -> inferFieldRecursive(pr.child(), unknown, visited);
			// PeeledViaTa: the TA cross-fusion is integer (±1) so doesn't narrow; the
			// field is the composition of the diag (cube) and corner.
			case Lineage.PeeledViaTa t -> combine(
					inferFieldRecursive(t.cube(), unknown, visited),
					inferFieldRecursive(t.corner(), unknown, visited));
		};
	}

	private static final Pattern LEAF_SHAPE_REF = Pattern.compile(
			"(\\d+)x(\\d+)x(\\d+)(?:@[0-9a-f]+)?(?:-direct|-naive)?");

	private Optional<Field> classifyLeaf(String ref, List<String> unknown,
			java.util.Set<Path> visited) {
		// Named base: Strassen<n,m,p>=r, Laderman<3,3,3>=23, etc. → Z-coefficient.
		if (ref.contains("<") && ref.contains("=")) {
			return Optional.of(Field.Z);
		}
		// Terminal ground-truth leaf ("naive-NxMxP", synthesised by trivialOneAxis /
		// naive grids): all-ones INTEGER coefficients, valid over every field → Z.
		// Without this it fell through to unknownLeaves and the whole lineage's
		// field inference bailed (mirrors atomFields' naive branch; same blind spot
		// as the SELF-SHAPE naive exemption in RecursiveMaterialiser).
		if (ref.contains("naive")) {
			return Optional.of(Field.Z);
		}
		// Shape ref: resolve to a file via the catalog lookup by shape.
		Matcher m = LEAF_SHAPE_REF.matcher(ref);
		if (m.matches()) {
			int n = Integer.parseInt(m.group(1));
			int mm = Integer.parseInt(m.group(2));
			int p = Integer.parseInt(m.group(3));
			Optional<Path> path = findFile(n, mm, p);
			if (path.isPresent()) {
				return classifyFileWithLineage(path.get(), unknown, visited);
			}
			unknown.add(ref);
			return Optional.empty();
		}
		// Filename-stem ref: look up by stem across the index.
		Optional<Path> stemMatch = findFileByStem(ref);
		if (stemMatch.isPresent()) {
			return classifyFileWithLineage(stemMatch.get(), unknown, visited);
		}
		unknown.add(ref);
		return Optional.empty();
	}

	private Optional<Field> classifyFileWithLineage(Path path, List<String> unknown,
			java.util.Set<Path> visited) {
		// Avoid cycles (a stub referencing itself transitively).
		if (!visited.add(path)) {
			return Optional.of(fieldFromContent(path));
		}
		// Try to read the file's own lineage and recurse — leaves of a
		// composed stub might themselves be composed.
		try {
			tools.jackson.databind.JsonNode jroot = SchemeIO.parseJson(path.toFile());
			Optional<Lineage.Node> nested = SchemeIO.readLineage(jroot);
			if (nested.isPresent()) {
				return inferFieldRecursive(nested.get(), unknown, visited);
			}
		} catch (Exception e) {
			// Fall through to content-based field classification.
		}
		return Optional.of(fieldFromContent(path));
	}

	private Optional<Path> findFileByStem(String stem) {
		String wanted = stem.endsWith(".json") ? stem : stem + ".json";
		// Walk the cached index buckets; cheap because they're already in memory.
		for (List<FileEntry> bucket : fileIndex().values()) {
			for (FileEntry e : bucket) {
				if (e.path.getFileName().toString().equals(wanted)) {
					return Optional.of(e.path);
				}
			}
		}
		return Optional.empty();
	}

	/**
	 * Combine two leaf-field opinions per the lattice rules:
	 * <ul>
	 *   <li>both characteristic-0 → MAX in Z &lt; Q &lt; R &lt; C ordering;</li>
	 *   <li>both same finite field → that field;</li>
	 *   <li>{@code Optional.empty()} on either side → empty (propagate
	 *       "unknown" conservatively);</li>
	 *   <li>F2 mixed with char-0 → F2 (per project convention — assume the
	 *       char-0 leaf is Z-coefficient and reduces mod 2);</li>
	 *   <li>F3 mixed with char-0 → F3 (same convention);</li>
	 *   <li>F2 mixed with F3 → empty (no common ground).</li>
	 * </ul>
	 */
	public static Optional<Field> combine(Optional<Field> a, Optional<Field> b) {
		if (a.isEmpty() || b.isEmpty()) return Optional.empty();
		Field fa = a.get(), fb = b.get();
		if (fa == fb) return Optional.of(fa);
		boolean aFinite = (fa == Field.F2 || fa == Field.F3);
		boolean bFinite = (fb == Field.F2 || fb == Field.F3);
		if (aFinite && bFinite) return Optional.empty();  // F2 vs F3
		if (aFinite) return Optional.of(fa);  // F-side dominates char-0
		if (bFinite) return Optional.of(fb);
		// Both characteristic 0: pick max in the Z < Q < R < C lattice.
		return Optional.of(maxChar0(fa, fb));
	}

	private static Field maxChar0(Field a, Field b) {
		int ra = char0Rank(a), rb = char0Rank(b);
		return (ra >= rb) ? a : b;
	}

	private static int char0Rank(Field f) {
		return switch (f) {
			case Z -> 0;
			case Q -> 1;
			case R -> 2;
			case C -> 3;
			default -> -1;  // unreachable on the char-0 branch
		};
	}

	private static final org.slf4j.Logger LOG = org.slf4j.LoggerFactory.getLogger(FieldAwareLookup.class);

	/**
	 * Known-broken scheme files: spot-check failures from
	 * {@link AuditAlphaTensorSchemes}. These are AT-Z derived-improvement
	 * entries whose factor matrices in our local import don't actually
	 * compute matmul. Until the import is fixed (or the files are
	 * re-fetched from a clean source), the lookup skips them so
	 * downstream constructions don't propagate the bug.
	 *
	 * <p>Tracking: see ROADMAP.md "AT bulk-import provenance audit" and
	 * DISCOVERIES_PENDING_ANALYSIS.md for the audit script and remediation
	 * plan.</p>
	 */
	/** Test if a scheme file is on the known-broken denylist. */
	public static boolean isKnownBroken(String filename) {
		return KNOWN_BROKEN_FILES.contains(filename);
	}

	private static final java.util.Set<String> KNOWN_BROKEN_FILES = java.util.Set.of();

	/**
	 * One-time index of every scheme file on disk, keyed by sorted
	 * "{@code nxmxp}" shape. Per-shape candidates are pre-sorted by
	 * filename-claimed rank ascending so {@link #find} can short-circuit
	 * on the first orientable hit per allowed field.
	 *
	 * <p>Built lazily by {@link #fileIndex()} the first time a lookup is
	 * needed — avoids walk-on-init cost for constructors that never use
	 * {@code find}. Shared across all {@link FieldAwareLookup} instances
	 * pointing at the same {@code schemesRoot}.</p>
	 */
	private static final Map<Path, Map<String, List<FileEntry>>> INDEX_CACHE = new ConcurrentHashMap<>();

	/**
	 * Cache of parsed schemes (stored in {@link CompactScheme} form for
	 * memory efficiency) by canonical file path. Bounded LRU sized for
	 * the working set of a recursive sweep — each entry is small (kBs
	 * for typical Strassen-style integer schemes, more for Q-rational
	 * dense ones). Schemes are {@link CompactScheme#expand expanded} to
	 * full {@link NonCubicBilinearAlgorithm} only at the use site.
	 */
	private static final int PARSE_CACHE_CAPACITY = 16*4096;
	private static final Cache<Path, CompactScheme> PARSE_CACHE = CacheBuilder.newBuilder().maximumSize(PARSE_CACHE_CAPACITY).build();

	private record FileEntry(Path path, eu.solven.matmul.algebra.Field field, int rank,
			boolean commutative) {}

	private Map<String, List<FileEntry>> fileIndex() {
		return INDEX_CACHE.computeIfAbsent(schemesRoot, FieldAwareLookup::buildIndex);
	}

	private static Map<String, List<FileEntry>> buildIndex(Path root) {
		// CONTENT-DRIVEN (2026-06): shape (n), rank (m), field (fields[]) and
		// commutativity are read from each scheme's JSON, NOT its filename. Filenames
		// are now pure labels — no `_m{rank}` token, no `_Z_`/`_ZT_` field token. This
		// is what lets the catalog be renamed freely. Requires `fields` stamped on every
		// scheme (StampFields / the materialiser); imports always carried it.
		// FAST PATH: read shape/rank/fields/commutative from docs/catalog.json for the
		// files it covers (no per-file JSON parse); parse ONLY the disk files the
		// manifest OMITS (~2k worse-rank derived alternates today; → ~0 once the catalog
		// keep-closure prunes/keeps the reference-set). For a non-canonical root (a
		// staging schemes dir the manifest does not describe) the relative-path lookup
		// simply misses and every file is parsed, as before.
		boolean canonical = root.normalize().equals(Path.of("src/main/resources/schemes").normalize());
		Map<String, JsonNode> manifest = canonical ? manifestByFile() : Map.of();

		List<Path> paths;
		try (Stream<Path> walk = Files.walk(root)) {
			paths = walk.filter(p -> p.getFileName().toString().endsWith(".json"))
					.collect(java.util.stream.Collectors.toList());
		} catch (IOException e) {
			throw new RuntimeException("walking " + root, e);
		}
		Map<String, List<FileEntry>> out = new java.util.concurrent.ConcurrentHashMap<>();
		java.util.concurrent.atomic.AtomicInteger fromManifest = new java.util.concurrent.atomic.AtomicInteger();
		java.util.concurrent.atomic.AtomicInteger parsed = new java.util.concurrent.atomic.AtomicInteger();
		paths.parallelStream().forEach(p -> {
			String name = p.getFileName().toString();
			// Corrupted stubs are treated as ABSENT for gating (#68).
			if (KNOWN_BROKEN_FILES.contains(name) || manifestCorrupted().contains(name)) return;
			try {
				// Prefer the manifest entry (no parse); fall back to parsing the file.
				JsonNode meta = manifest.get(root.relativize(p).toString());
				int[] dims;
				int rank;
				java.util.List<String> tags;
				boolean commutative;
				if (meta != null) {
					JsonNode fmt = meta.get("format");
					if (fmt == null || !fmt.isArray() || fmt.size() != 3) return;
					dims = new int[] { fmt.get(0).asInt(), fmt.get(1).asInt(), fmt.get(2).asInt() };
					rank = meta.has("rank") ? meta.get("rank").asInt() : -1;
					tags = SchemeIO.fieldTags(meta);
					commutative = meta.has("commutative") && meta.get("commutative").asBoolean();
					fromManifest.incrementAndGet();
				} else {
					JsonNode r = SchemeIO.parseJson(p.toFile());
					JsonNode n = r.get("n");
					if (n == null || !n.isArray() || n.size() != 3) return;
					dims = new int[] { n.get(0).asInt(), n.get(1).asInt(), n.get(2).asInt() };
					rank = r.has("m") ? r.get("m").asInt() : (r.has("rank") ? r.get("rank").asInt() : -1);
					tags = SchemeIO.fieldTags(r);
					commutative = r.has("commutative") && r.get("commutative").asBoolean();
					parsed.incrementAndGet();
				}
				Arrays.sort(dims);
				if (rank < 0) return;
				if (tags.isEmpty()) {
					// A scheme with no verified fields[] is a DATA BUG (backfill it via
					// BackfillMissingFields). Treat it as ABSENT for field gating rather
					// than silently defaulting to Z — that default (the ⟨3,3,8⟩ bogus-Z
					// bug) let a Q-only derived scheme satisfy Z queries it never belonged
					// to. Fail loud + skip so it can never leak into the wrong field.
					LOG.warn("buildIndex: {} has no fields[]; treating as ABSENT (run BackfillMissingFields)", name);
					return;
				}
				String key = dims[0] + "x" + dims[1] + "x" + dims[2];
				eu.solven.matmul.algebra.Field fld = fieldFromTags(tags);
				out.computeIfAbsent(key, k -> java.util.Collections.synchronizedList(new ArrayList<>()))
						.add(new FileEntry(p, fld, rank, commutative));
			} catch (Exception e) {
				// unreadable / malformed → skip (transparently absent, as before)
			}
		});
		// Sort each shape's candidates by rank ascending, then by path for a STABLE
		// tie-break (the parallel walk above is unordered — without this, "catalog-best
		// at a shape" among same-rank schemes would be non-deterministic across runs).
		for (List<FileEntry> entries : out.values()) {
			entries.sort(Comparator.comparingInt(FileEntry::rank)
					.thenComparing(e -> e.path().getFileName().toString()));
		}
		LOG.info("buildIndex (content-driven): {} shapes / {} files ({} from manifest, {} parsed).",
				out.size(), paths.size(), fromManifest.get(), parsed.get());
		return out;
	}

	/** Narrowest single field from a scheme's verified {@code fields[]} tags, for the
	 *  index field filter ({@code requested.fallbackChain().contains(e.field)}). Char-0
	 *  narrowest first (Z⊂Q⊂R⊂C); F₂/F₃ only when no char-0 tag is present. */
	private static eu.solven.matmul.algebra.Field fieldFromTags(java.util.List<String> tags) {
		if (tags.contains("Z")) return eu.solven.matmul.algebra.Field.Z;
		if (tags.contains("Q")) return eu.solven.matmul.algebra.Field.Q;
		if (tags.contains("R")) return eu.solven.matmul.algebra.Field.R;
		if (tags.contains("C")) return eu.solven.matmul.algebra.Field.C;
		if (tags.contains("F2")) return eu.solven.matmul.algebra.Field.F2;
		if (tags.contains("F3")) return eu.solven.matmul.algebra.Field.F3;
		return eu.solven.matmul.algebra.Field.Z;
	}

	/** Narrowest field of a scheme file, read from its JSON {@code fields[]}
	 *  content (never the filename). Defaults to Z (the catalog is
	 *  integer-coefficient by default) when the file is unreadable. */
	static eu.solven.matmul.algebra.Field fieldFromContent(Path path) {
		try {
			return fieldFromTags(SchemeIO.fieldTags(SchemeIO.parseJson(path.toFile())));
		} catch (Exception e) {
			return eu.solven.matmul.algebra.Field.Z;
		}
	}

	private static NonCubicBilinearAlgorithm parseCached(Path path) {
		CompactScheme hit = PARSE_CACHE.getIfPresent(path);
		if (hit != null) return hit.expand();
		File f = path.toFile();
		NonCubicBilinearAlgorithm alg;
		try {
			// Parse once; stubs (lineage-only) are not yet handled by the
			// runtime materialiser — return null so they're transparently
			// skipped (the catalog still contains the scheme as a derivable
			// artifact; on-demand materialisation lands in a follow-up).
			tools.jackson.databind.JsonNode root = SchemeIO.parseJson(f);
			if (SchemeIO.isStub(root)) return null;
			alg = SchemeIO.isReduced(root) ? SchemeIO.readReduced(root) : SchemeIO.read(root);
		} catch (Exception e) {
			return null;  // sentinel — broken file
		}
		PARSE_CACHE.put(path, CompactScheme.of(alg));
		return alg;  // return the just-loaded expanded form to avoid double-allocation
	}

	@Override
	public Optional<NonCubicBilinearAlgorithm> find(int n, int m, int p) {
		return findWithSource(n, m, p).map(WithSource::alg);
	}

	/**
	 * Rank-only fast path: returns the catalog rank for shape
	 * {@code ⟨n,m,p⟩} (via sorted-shape key + field-chain filter) WITHOUT
	 * parsing or expanding the scheme's factor matrices. Avoids the
	 * {@code CompactScheme.expand()} allocation explosion that dominates
	 * inner-loop search cost when callers only need the rank.
	 *
	 * <p>Rank is permutation-invariant under axis-permutation, so no
	 * orientation step is needed — the {@code FileEntry.rank} value
	 * applies to any S₃ permutation of the canonical sorted shape.
	 *
	 * @return {@code rank} as an int, or {@code Recombination.SotaResolver.UNKNOWN_RANK}
	 *         if no entry is found in the field-chain. The sentinel value
	 *         is what callers expect for "unknown" in the SOTA resolver.
	 */
	public int findRank(int n, int m, int p) {
		int[] sorted = { n, m, p };
		Arrays.sort(sorted);
		String key = sorted[0] + "x" + sorted[1] + "x" + sorted[2];
		List<FileEntry> candidates = fileIndex().get(key);
		if (candidates != null) {
			List<eu.solven.matmul.algebra.Field> allowedFields = field.fallbackChain();
			// candidates are sorted by rank ascending — return the first matching field
			// that is NOT commutative. Commutative ranks (Waksman, Rosowski, Makarov)
			// are invalid for non-commutative recombination — using them silently
			// would under-count the true rank and produce bogus "wins" (see
			// 2026-06-02 incident where mask5 Winograd "discovered"
			// ⟨17,17,17⟩=2868 by importing Waksman/Rosowski sub-shape ranks).
			for (FileEntry e : candidates) {
				if (!allowedFields.contains(e.field)) continue;
				if (e.commutative) continue;
				return e.rank;
			}
		}
		// No non-commutative catalog scheme over this field → the NAÏVE rank a·b·c, which is
		// ALWAYS achievable as a formula (the trivial scheme; no dense representation, so even
		// a ⟨1,1234,5678⟩ that can never be on disk has rank 1·1234·5678). NEVER the old
		// MAX/100 sentinel: a sentinel poisoned sum-based lower bounds (AllocationOptimizer
		// dropped good bases over a ⟨…,1,…⟩ relaxation block → the ⟨12,13,13⟩=1274 loss);
		// naïve is a valid, constructive bound. (user 2026-06-15)
		return naiveRank(sorted[0], sorted[1], sorted[2]);
	}

	/** Distinct in-range non-degenerate shapes already warned about (dedup + cap). */
	private static final java.util.Set<String> WARNED_MISSES = ConcurrentHashMap.newKeySet();
	private static final int WARN_CAP = 100;

	/** Naïve {@code a·b·c} rank (a ≤ b ≤ c) with a throttled "thin catalog" WARN. A miss is
	 *  EXPECTED — and stays silent — for a degenerate shape ({@code a == 1}, the B&B's
	 *  relaxation blocks) or an out-of-range one ({@code c > MAX_DIM}, too big to ever store);
	 *  a non-degenerate IN-range miss is a real gap worth surfacing once. */
	private static int naiveRank(int a, int b, int c) {
		long naive = (long) a * b * c;
		if (a >= 2 && c <= CatalogLimits.MAX_DIM && WARNED_MISSES.size() < WARN_CAP
				&& WARNED_MISSES.add(a + "x" + b + "x" + c)) {
			LOG.warn("findRank: no scheme for ⟨{},{},{}⟩ (in range) — returning naïve {} (thin catalog)",
					a, b, c, naive);
			if (WARNED_MISSES.size() == WARN_CAP) {
				LOG.warn("findRank: further thin-catalog warnings suppressed (cap {})", WARN_CAP);
			}
		}
		return naive > Integer.MAX_VALUE ? Integer.MAX_VALUE : (int) naive;
	}

	/**
	 * Rank lookup that <strong>includes</strong> commutative-only schemes
	 * (Waksman 1970, Rosowski 2019 Thm 2/3, Makarov 1986, Islam 2009). Use
	 * this at the TOP LEVEL of a search when the caller accepts a commutative
	 * algorithm (e.g. computing scalar matmul, not recursive matmul over a
	 * matrix ring).
	 *
	 * <p>NEVER use as the SOTA resolver inside a recombination — see the
	 * 2026-06-02 ⟨17,17,17⟩=2868 incident: a commutative sub-product rank
	 * silently pulled into NC recombination under-counts the true rank.
	 * The standard {@link #findRank} keeps NC purity for that case.</p>
	 *
	 * @return smallest rank including commutative entries, or
	 *         {@code Recombination.SotaResolver.UNKNOWN_RANK} if no entry is found
	 */
	public int findRankAllowCommutative(int n, int m, int p) {
		int[] sorted = { n, m, p };
		Arrays.sort(sorted);
		String key = sorted[0] + "x" + sorted[1] + "x" + sorted[2];
		List<FileEntry> candidates = fileIndex().get(key);
		if (candidates == null) return Recombination.SotaResolver.UNKNOWN_RANK;
		List<eu.solven.matmul.algebra.Field> allowedFields = field.fallbackChain();
		for (FileEntry e : candidates) {
			if (!allowedFields.contains(e.field)) continue;
			return e.rank; // first match (sorted by rank ascending)
		}
		return Recombination.SotaResolver.UNKNOWN_RANK;
	}

	/**
	 * Rank lookup restricted to <strong>commutative-only</strong> entries.
	 * Useful for "what's the best commutative bound at this shape?" queries,
	 * and to assess the commutative-vs-NC gap at a given format.
	 *
	 * @return smallest commutative-only rank, or {@code Recombination.SotaResolver.UNKNOWN_RANK}
	 *         if no commutative entry is found
	 */
	public int findRankCommutativeOnly(int n, int m, int p) {
		int[] sorted = { n, m, p };
		Arrays.sort(sorted);
		String key = sorted[0] + "x" + sorted[1] + "x" + sorted[2];
		List<FileEntry> candidates = fileIndex().get(key);
		if (candidates == null) return Recombination.SotaResolver.UNKNOWN_RANK;
		List<eu.solven.matmul.algebra.Field> allowedFields = field.fallbackChain();
		for (FileEntry e : candidates) {
			if (!allowedFields.contains(e.field)) continue;
			if (!e.commutative) continue;
			return e.rank;
		}
		return Recombination.SotaResolver.UNKNOWN_RANK;
	}

	/**
	 * Like {@link #find} but also returns the source file path so callers
	 * can read other metadata from the on-disk JSON (e.g. the file's own
	 * {@code lineage} field to recursively reconstruct deep lineage trees).
	 */
	public Optional<WithSource> findWithSource(int n, int m, int p) {
		int[] sorted = { n, m, p };
		Arrays.sort(sorted);
		String key = sorted[0] + "x" + sorted[1] + "x" + sorted[2];
		List<FileEntry> candidates = fileIndex().get(key);
		if (candidates == null) return Optional.empty();
		List<eu.solven.matmul.algebra.Field> allowedFields = field.fallbackChain();
		for (FileEntry e : candidates) {
			if (!allowedFields.contains(e.field)) continue;
			// Skip commutative entries — they're invalid for NC recombination
			// even when the field nominally matches. Keeps find() consistent
			// with findRank's fast-path filter.
			if (e.commutative) continue;
			NonCubicBilinearAlgorithm alg = parseCached(e.path);
			if (alg == null) continue;
			Optional<NonCubicBilinearAlgorithm> oriented = alg.orientAs(n, m, p);
			if (oriented.isPresent()) {
				return Optional.of(new WithSource(oriented.get(), e.path));
			}
		}
		return Optional.empty();
	}

	/** Pairs a found scheme with the file path it was loaded from. */
	public record WithSource(NonCubicBilinearAlgorithm alg, Path path) {}

	/**
	 * Return the on-disk file path for the best catalog entry of shape
	 * {@code ⟨n,m,p⟩}, WITHOUT parsing it. Used by the lineage replayer
	 * to discover stub files whose {@link #parseCached} returns null
	 * (stubs are skipped from the lookup but their files exist and carry
	 * lineage). Returns the first candidate matching the field-fallback
	 * chain; orientation-mismatched entries are NOT filtered here (since
	 * we have no algorithm to orient).
	 */
	public Optional<Path> findFile(int n, int m, int p) {
		int[] sorted = { n, m, p };
		Arrays.sort(sorted);
		String key = sorted[0] + "x" + sorted[1] + "x" + sorted[2];
		List<FileEntry> candidates = fileIndex().get(key);
		if (candidates == null) return Optional.empty();
		List<eu.solven.matmul.algebra.Field> allowedFields = field.fallbackChain();
		for (FileEntry e : candidates) {
			if (allowedFields.contains(e.field)) return Optional.of(e.path);
		}
		return Optional.empty();
	}

	/**
	 * All on-disk files for {@code ⟨n,m,p⟩} in the field-fallback chain,
	 * <strong>rank-ascending</strong> (the {@code fileIndex} bucket is pre-sorted
	 * by filename-claimed rank). Unlike {@link #findFile}, which returns only the
	 * single lowest-rank candidate, this lets a caller fall through to the
	 * next-best file when the lowest is e.g. a non-replayable stub — so the best
	 * <em>usable</em> scheme is found rather than giving up on the first.
	 */
	public List<Path> findFiles(int n, int m, int p) {
		int[] sorted = { n, m, p };
		Arrays.sort(sorted);
		String key = sorted[0] + "x" + sorted[1] + "x" + sorted[2];
		List<FileEntry> candidates = fileIndex().get(key);
		if (candidates == null) return List.of();
		List<eu.solven.matmul.algebra.Field> allowedFields = field.fallbackChain();
		List<Path> out = new java.util.ArrayList<>();
		for (FileEntry e : candidates) {
			if (allowedFields.contains(e.field)) out.add(e.path);
		}
		return out;
	}

	/**
	 * Like {@link #findFiles} but EXCLUDES commutative-only entries — the
	 * rank-ascending files of shape {@code ⟨n,m,p⟩} that are valid as
	 * <strong>non-commutative</strong> leaves (same filter {@link #findRank}
	 * applies). Use when a leaf must lift to recursive matmul over a non-commutative
	 * ring: the global rank-best may be a commutative Waksman/Rosowski scheme (e.g.
	 * ⟨26,3,3⟩=159 Rosowski vs the NC 175) that {@link #findFile} would wrongly
	 * surface — the 2026-06-02 ⟨17,17,17⟩=2868 commutative-leaf footgun.
	 */
	public List<Path> findFilesNonCommutative(int n, int m, int p) {
		int[] sorted = { n, m, p };
		Arrays.sort(sorted);
		String key = sorted[0] + "x" + sorted[1] + "x" + sorted[2];
		List<FileEntry> candidates = fileIndex().get(key);
		if (candidates == null) return List.of();
		List<eu.solven.matmul.algebra.Field> allowedFields = field.fallbackChain();
		List<Path> out = new java.util.ArrayList<>();
		for (FileEntry e : candidates) {
			if (allowedFields.contains(e.field) && !e.commutative) out.add(e.path);
		}
		return out;
	}

	/**
	 * Resolve a scheme by {@link SchemeIO#contentHash} at shape ⟨n,m,p⟩ — the
	 * precise, collision-free alternative to the ambiguous shape-ref (which picks
	 * the rank-best and can return the wrong content, e.g. a bud-poor sibling).
	 * Hashes each candidate's stored form and matches {@code hash} (full or short
	 * prefix). An exact full-hash match is returned immediately (unique); on a
	 * short-prefix collision the first match wins (the caller may disambiguate
	 * further by property — pass a longer hash to make it exact). Empty if nothing
	 * matches, so the caller can fall back to the shape-ref.
	 */
	public Optional<WithSource> findByHash(int n, int m, int p, String hash) {
		WithSource prefixHit = null;
		for (Path path : findFiles(n, m, p)) {
			NonCubicBilinearAlgorithm alg = parseCached(path);
			if (alg == null) {
				continue;  // stub — no matrices to hash
			}
			Optional<NonCubicBilinearAlgorithm> oriented = alg.orientAs(n, m, p);
			if (oriented.isEmpty()) {
				continue;
			}
			// Match the content hash registered for this shape. MOST refs hash the STORED
			// form, but orientation-specific refs hash the base AFTER orienting it to the
			// frame the construction needs — e.g. a SerendipitousProduct base hashed as its
			// ⟨6,5,5⟩ product-frame while the file is stored ⟨5,5,6⟩. Check BOTH so such a
			// ref still resolves to the EXACT (bud-rich) base, not the rank-best bud-free
			// sibling — which is what made ⟨18,20,30⟩/⟨24,24,27⟩ unreplayable (NoSuchElement
			// in productViaBuds) and blocked the fill. Backward-compatible (stored still wins).
			String storedHash = SchemeIO.contentHash(alg);
			String orientedHash = SchemeIO.contentHash(oriented.get());
			boolean exact = storedHash.equals(hash) || orientedHash.equals(hash);
			boolean prefix = storedHash.startsWith(hash) || orientedHash.startsWith(hash);
			if (exact || prefix) {
				WithSource hit = new WithSource(oriented.get(), path);
				if (exact) {
					return Optional.of(hit);
				}
				if (prefixHit == null) {
					prefixHit = hit;  // first short-prefix match; keep scanning for an exact one
				}
			}
		}
		return Optional.ofNullable(prefixHit);
	}

	/**
	 * Clears the global index + parse cache. Use as the heavy hammer
	 * after schemes/ has been mutated in bulk (manual edits, file
	 * deletes, external generation). For single-file writes prefer
	 * {@link #onSchemeWritten} — surgical update, no re-walk.
	 */
	public static void invalidateCache() {
		INDEX_CACHE.clear();
		PARSE_CACHE.invalidateAll();
		// Drop the manifest-derived caches too — a bulk catalog mutation may have
		// regenerated docs/catalog.json, so the next buildIndex must re-read it.
		MANIFEST_BY_FILE = null;
		MANIFEST_FIELDS = null;
	}

	/**
	 * Surgical cache update: a new scheme file was just written to
	 * disk. If the index has been built for any registered root that
	 * contains {@code file}, splice the new entry into the existing
	 * shape-bucket (re-sorted by rank) so the next {@link #find} call
	 * sees it without a full directory re-walk.
	 *
	 * <p>Closure-loop materialisers (write-then-find-then-write-again
	 * patterns like {@code MaterializeClosureLoop}) rely on this — if
	 * the index were stale, an iteration's freshly-written scheme would
	 * be invisible to the next iteration's lookups.</p>
	 *
	 * <p>Called automatically by
	 * {@code SchemeIO.write(NonCubicBilinearAlgorithm, File[, Lineage.Node])}
	 * — callers writing through other paths can invoke this directly.</p>
	 */
	public static void onSchemeWritten(File file) {
		String name = file.getName();
		if (!name.endsWith(".json") || KNOWN_BROKEN_FILES.contains(name)) return;
		// Invalidate any stale parse-cache entry for this specific path
		// (in case it was loaded before the overwrite).
		PARSE_CACHE.invalidate(file.toPath());
		PARSE_CACHE.invalidate(file.toPath().toAbsolutePath());

		// CONTENT-DRIVEN (mirrors buildIndex): shape (n), rank (m), field
		// (fields[]) and commutativity are read from the just-written JSON —
		// never the filename, which is a pure label.
		int[] dims;
		int rank;
		eu.solven.matmul.algebra.Field fld;
		boolean commutative;
		try {
			JsonNode r = SchemeIO.parseJson(file);
			JsonNode n = r.get("n");
			if (n == null || !n.isArray() || n.size() != 3) return;
			dims = new int[] { n.get(0).asInt(), n.get(1).asInt(), n.get(2).asInt() };
			rank = r.has("m") ? r.get("m").asInt() : (r.has("rank") ? r.get("rank").asInt() : -1);
			if (rank < 0) return;
			fld = fieldFromTags(SchemeIO.fieldTags(r));
			commutative = r.has("commutative") && r.get("commutative").asBoolean();
		} catch (Exception e) {
			return; // unreadable / malformed → transparently absent
		}
		Arrays.sort(dims);
		String key = dims[0] + "x" + dims[1] + "x" + dims[2];
		Path absFile = file.toPath().toAbsolutePath();

		for (var e : INDEX_CACHE.entrySet()) {
			Path absRoot = e.getKey().toAbsolutePath();
			if (!absFile.startsWith(absRoot)) continue;
			Map<String, List<FileEntry>> index = e.getValue();
			List<FileEntry> bucket = index.computeIfAbsent(key, k -> new ArrayList<>());
			// Replace any existing entry for the same path (file was overwritten).
			bucket.removeIf(en -> en.path.toAbsolutePath().equals(absFile));
			bucket.add(new FileEntry(file.toPath(), fld, rank, commutative));
			bucket.sort(Comparator.comparingInt(FileEntry::rank)
					.thenComparing(en -> en.path().getFileName().toString()));
		}
	}
}
