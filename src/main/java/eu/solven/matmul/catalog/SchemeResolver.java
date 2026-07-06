package eu.solven.matmul.catalog;

import java.io.File;
import java.nio.file.Path;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Resolve a scheme FILE by its identity (shape + source), tolerant of the folder split
 * (known/derived/curated), the 2026 filename rename
 * ({@code {shape}-r{rank}-{note}-{hash7}.json}) and historical {@code _a/_b}-token drift.
 *
 * <p>Replaces hardcoded scheme paths and {@code BlockSplitSearch.resolvePoolFile}'s
 * stem-glob, both of which break when a file is moved or renamed. The SHAPE is taken
 * from the catalog index (content {@code n}); among same-shape files the source token
 * disambiguates (the rename preserves it as the filename {@code note}). Accepts an old
 * OR new hint path/stem — only the embedded {@code NxMxP} and the alpha prefix matter.</p>
 *
 * <p><b>Fail-loud discipline</b>: an unmatched non-blank source token falls back to the
 * first scheme at the shape — but WARNS, because the caller asked for a specific scheme
 * and is getting a different one (the 2026-06 purge/renames silently re-targeted two
 * ⟨17,17,17⟩ tests this way). Callers that verify a SPECIFIC scheme must use
 * {@link #byHintStrict} and skip when it returns {@code null}.</p>
 */
public final class SchemeResolver {
	private SchemeResolver() {}

	private static final Logger log = LoggerFactory.getLogger(SchemeResolver.class);

	private static final Pattern SHAPE = Pattern.compile("(\\d+)x(\\d+)x(\\d+)");
	private static volatile FieldAwareLookup defaultLookup;
	private static volatile FieldAwareLookup f2Lookup;
	private static volatile FieldAwareLookup f3Lookup;

	private static FieldAwareLookup lookup() {
		FieldAwareLookup lk = defaultLookup;
		if (lk == null) {
			// C is the widest CHAR-0 field (Z ⊂ Q ⊂ R ⊂ C) — but NOT all-covering:
			// F₂/F₃ are their own universes (no char-0 inclusion), so an F2-only file
			// (e.g. AlphaTensor ⟨4,4,4⟩=47) is absent from a C index. Candidate
			// enumeration must UNION the universes — see findFilesAllFields.
			lk = new FieldAwareLookup("C");
			defaultLookup = lk;
		}
		return lk;
	}

	/** Files at ⟨n,m,p⟩ across EVERY field universe: char-0 (C index) first, then the
	 *  F2-only / F3-only extras. Z-stamped schemes carry all six tags, so the de-dup
	 *  keeps their C-index entry only. Without this union, a hint naming an F2-only
	 *  scheme silently resolved to an arbitrary char-0 sibling (a Waksman commutative
	 *  file for the AlphaTensor ⟨4,4,4⟩=47 hint). */
	private static List<Path> findFilesAllFields(int n, int m, int p) {
		FieldAwareLookup f2 = f2Lookup, f3 = f3Lookup;
		if (f2 == null) { f2 = new FieldAwareLookup("F2"); f2Lookup = f2; }
		if (f3 == null) { f3 = new FieldAwareLookup("F3"); f3Lookup = f3; }
		java.util.LinkedHashSet<Path> out = new java.util.LinkedHashSet<>(lookup().findFiles(n, m, p));
		out.addAll(f2.findFiles(n, m, p));
		out.addAll(f3.findFiles(n, m, p));
		return new java.util.ArrayList<>(out);
	}

	/** Resolve a hint path/stem (old or new) to the actual file. Returns a {@code File}
	 *  at the hint path verbatim if the hint can't be parsed or nothing resolves — so the
	 *  caller's existing not-found handling still fires. */
	public static File byHint(String hint) {
		File exact = new File(hint);
		if (exact.exists()) {
			return exact;
		}
		ParsedHint ph = parseHint(exact);
		if (ph == null) {
			return exact;
		}
		File r = byShapeAndSource(ph.n, ph.m, ph.p, ph.token);
		return r != null ? r : exact;
	}

	/**
	 * Strict variant of {@link #byHint}: returns the file ONLY when it is the scheme the
	 * hint identifies — exact path, or a same-shape file whose name matches the hint's
	 * source token. Returns {@code null} when the token matches nothing (e.g. the scheme
	 * was purged from the catalog), instead of silently handing back an arbitrary
	 * same-shape file. Use this in tests that verify a specific import, paired with a
	 * JUnit assumption so the test SKIPS (not vacuously passes) when the subject is gone.
	 */
	public static File byHintStrict(String hint) {
		File exact = new File(hint);
		if (exact.exists()) {
			return exact;
		}
		ParsedHint ph = parseHint(exact);
		if (ph == null) {
			return null;
		}
		List<Path> cands = findFilesAllFields(ph.n, ph.m, ph.p);
		File matched = matchToken(cands, ph.token);
		if (matched == null) {
			log.warn("byHintStrict: no scheme at ⟨{},{},{}⟩ matches token '{}' (hint {}); returning null",
					ph.n, ph.m, ph.p, ph.token, hint);
		}
		return matched;
	}

	private record ParsedHint(int n, int m, int p, String token) {}

	private static ParsedHint parseHint(File exact) {
		String stem = exact.getName().replaceFirst("\\.json$", "");
		Matcher m = SHAPE.matcher(stem);
		if (!m.find()) {
			return null;
		}
		int n = Integer.parseInt(m.group(1));
		int mm = Integer.parseInt(m.group(2));
		int p = Integer.parseInt(m.group(3));
		// Source token: alpha prefix BEFORE the shape (old names: "strassen-2x2x2…"), or
		// the note AFTER it (new names: "2x2x2-r7-strassen-…"). Take whichever is non-empty.
		String before = stem.substring(0, m.start()).replaceAll("[-_]+$", "");
		String after = stem.substring(m.end()).replaceFirst("^[-_]+(r\\d+[-_]+)?", "")
				.replaceFirst("[-_]+[0-9a-f]{6,}$", "");  // drop trailing hash
		String token = !before.isBlank() ? before : after;
		return new ParsedHint(n, mm, p, token);
	}

	/** First candidate whose filename matches the token (boundary match preferred), else null. */
	private static File matchToken(List<Path> cands, String sourceToken) {
		if (sourceToken == null || sourceToken.isBlank()) {
			return null;
		}
		String tok = sourceToken.toLowerCase();
		// Prefer an exact note-token boundary match, then a looser contains.
		for (Path c : cands) {
			if (c.getFileName().toString().toLowerCase().contains("-" + tok + "-")) {
				return c.toFile();
			}
		}
		for (Path c : cands) {
			if (c.getFileName().toString().toLowerCase().contains(tok)) {
				return c.toFile();
			}
		}
		return null;
	}

	/** Resolve the file at {@code ⟨n,m,p⟩} whose source/note matches {@code sourceToken}
	 *  (e.g. "strassen", "alphatensor_Z"); falls back to the first scheme at that shape
	 *  when the token is blank or unmatched — WITH a warning in the unmatched case, since
	 *  the caller asked for a specific scheme and is getting a different one.
	 *  {@code null} if no scheme exists at the shape. */
	public static File byShapeAndSource(int n, int m, int p, String sourceToken) {
		List<Path> cands = findFilesAllFields(n, m, p);
		if (cands.isEmpty()) {
			return null;
		}
		File matched = matchToken(cands, sourceToken);
		if (matched != null) {
			return matched;
		}
		File fallback = cands.get(0).toFile();
		if (sourceToken != null && !sourceToken.isBlank()) {
			log.warn("byShapeAndSource: no scheme at ⟨{},{},{}⟩ matches token '{}' — "
					+ "falling back to {} (of {} candidates). If the caller verifies a "
					+ "specific scheme, switch it to byHintStrict.",
					n, m, p, sourceToken, fallback.getName(), cands.size());
		}
		return fallback;
	}
}
