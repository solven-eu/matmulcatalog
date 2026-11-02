package eu.solven.matmul.additions;

import java.util.ArrayList;
import java.util.List;
import java.util.Optional;

import eu.solven.matmul.NonCubicBilinearAlgorithm;

import tools.jackson.databind.JsonNode;
import tools.jackson.databind.json.JsonMapper;
import tools.jackson.databind.node.ArrayNode;
import tools.jackson.databind.node.ObjectNode;

/**
 * The JSON codec of a scheme's straight-line programs — the one place that knows
 * both on-disk encodings, so a stored program can be read back and <em>replayed</em>
 * (nothing re-verified a stored {@code slp} before 2026-11-02; it was write-only).
 *
 * <ul>
 *   <li><b>Our object layout</b> ({@code slp: {additions, scalar_mults, optimality,
 *       optimality_method, source?, a:{inputs, ops, forms}, b:…, out:…}}), written by
 *       {@code MaterialiseAdditionsSlp} for the greedy CSE bound and by the upstream
 *       importers for a contributor's program: {@link #sideJson} / {@link #readSide},
 *       {@link #toJson} / {@link #readScheme}.</li>
 *   <li><b>The Perminov / we*star "reduced" layout</b> ({@code u_fresh}, {@code v_fresh},
 *       {@code w_fresh} + {@code u}/{@code v}/{@code w} rows of {@code {index, value}}
 *       references): {@link #fromFreshVars}. {@code SchemeIO.readReduced} expands the
 *       same layout into dense matrices; this reader keeps the <em>program</em>.</li>
 * </ul>
 *
 * <p>Honesty tier (CLAUDE.md optimality discipline): every program carries an
 * {@code optimality} of {@value #OPTIMALITY_BOUND} unless produced by an exact
 * solver; a contributor's reduced count is a bound too — a better schedule is
 * never excluded, only the replay is certified.</p>
 */
public final class SlpJson {
	private SlpJson() {}

	public static final String OPTIMALITY_BOUND = "bound";
	public static final String OPTIMALITY_PROVEN = "proven";

	/** The three per-side programs of one bilinear scheme: the {@code r} left factors
	 *  over the A-variables, the {@code r} right factors over the B-variables, the
	 *  {@code n·p} outputs over the products ({@code SchemeAdditiveComplexity} order). */
	public record SchemeSlp(Slp a, Slp b, Slp out) {
		public int additions() {
			return a.additions() + b.additions() + out.additions();
		}

		public int scalarMults() {
			return a.scalarMults() + b.scalarMults() + out.scalarMults();
		}

		/** True iff every side replays to the scheme's factor matrices exactly. */
		public boolean reconstructs(NonCubicBilinearAlgorithm alg) {
			return a.reconstructs(SchemeAdditiveComplexity.columnsAsForms(alg.denseU(), alg.r))
					&& b.reconstructs(SchemeAdditiveComplexity.columnsAsForms(alg.denseV(), alg.r))
					&& out.reconstructs(alg.denseW());
		}
	}

	// ── our object layout ──────────────────────────────────────────────────────

	/** One side as {@code {inputs, ops:[[kind, target, x, y|c]…], forms:[slot|-1…]}}. */
	public static ObjectNode sideJson(JsonMapper mapper, Slp slp) {
		ObjectNode side = mapper.createObjectNode();
		side.put("inputs", slp.nInputs());
		ArrayNode ops = mapper.createArrayNode();
		for (Slp.Op o : slp.ops()) {
			ArrayNode op = mapper.createArrayNode();
			op.add(String.valueOf(o.kind()));
			op.add(o.target());
			op.add(o.x());
			if (o.kind() == '*') {
				op.add(o.c());
			} else {
				op.add(o.y());
			}
			ops.add(op);
		}
		side.set("ops", ops);
		ArrayNode forms = mapper.createArrayNode();
		for (int fr : slp.formResult()) {
			forms.add(fr);
		}
		side.set("forms", forms);
		return side;
	}

	/** Inverse of {@link #sideJson}. */
	public static Slp readSide(JsonNode side) {
		int nInputs = side.get("inputs").asInt();
		List<Slp.Op> ops = new ArrayList<>();
		int adds = 0, muls = 0;
		for (JsonNode op : side.get("ops")) {
			char kind = op.get(0).asString().charAt(0);
			int target = op.get(1).asInt();
			int x = op.get(2).asInt();
			if (kind == '*') {
				ops.add(new Slp.Op('*', target, x, -1, op.get(3).asDouble()));
				muls++;
			} else if (kind == '+' || kind == '-') {
				ops.add(new Slp.Op(kind, target, x, op.get(3).asInt(), 0.0));
				adds++;
			} else {
				throw new IllegalArgumentException("unknown slp op kind " + kind);
			}
		}
		JsonNode formsNode = side.get("forms");
		int[] forms = new int[formsNode.size()];
		for (int i = 0; i < forms.length; i++) {
			forms[i] = formsNode.get(i).asInt();
		}
		return new Slp(nInputs, List.copyOf(ops), forms, adds, muls);
	}

	/**
	 * The {@code slp} object of a scheme file.
	 *
	 * @param optimality {@link #OPTIMALITY_BOUND} or {@link #OPTIMALITY_PROVEN}
	 * @param method     how the program was obtained (free text, e.g. the reducer)
	 * @param source     who produced it (the catalog's attribution key), or null for our own
	 */
	public static ObjectNode toJson(JsonMapper mapper, SchemeSlp s, String optimality, String method, String source) {
		ObjectNode slp = mapper.createObjectNode();
		slp.put("additions", s.additions());
		slp.put("scalar_mults", s.scalarMults());
		slp.put("optimality", optimality);
		slp.put("optimality_method", method);
		if (source != null) {
			slp.put("source", source);
		}
		slp.set("a", sideJson(mapper, s.a()));
		slp.set("b", sideJson(mapper, s.b()));
		slp.set("out", sideJson(mapper, s.out()));
		return slp;
	}

	/** Reads our object layout; empty for a missing / array-form / malformed {@code slp}. */
	public static Optional<SchemeSlp> readScheme(JsonNode slp) {
		if (slp == null || !slp.isObject() || !slp.has("a") || !slp.has("b") || !slp.has("out")) {
			return Optional.empty();
		}
		return Optional.of(new SchemeSlp(readSide(slp.get("a")), readSide(slp.get("b")), readSide(slp.get("out"))));
	}

	// ── Perminov / we*star reduced layout ──────────────────────────────────────

	/** Does the file carry fresh-variable programs ({@code u_fresh} / {@code v_fresh} / {@code w_fresh})? */
	public static boolean isFreshVarLayout(JsonNode root) {
		return root.has("u_fresh") || root.has("v_fresh") || root.has("w_fresh");
	}

	/**
	 * Read the fresh-variable programs of a reduced file into replayable SLPs.
	 *
	 * <p>Index space per side: {@code 0 … dim-1} the inputs (A-entries, B-entries, or
	 * the products for {@code w}), {@code dim + k} the k-th fresh variable. {@code u}/{@code v}
	 * rows are the {@code r} products (one linear form each); {@code w} rows are the
	 * {@code n·p} outputs in Perminov's <em>column-major</em> C order ({@code p = j·n + i}),
	 * mapped to our row-major order here — the same convention
	 * {@code SchemeIO.readReduced} applies to the matrices, so {@link SchemeSlp#reconstructs}
	 * is the consistency proof between the two readers.</p>
	 *
	 * <p>Additions = one per extra term of every fresh variable and every row; a
	 * negation alone is a scalar multiplication (free in the usual count), never an
	 * addition — the same convention as the contributor's {@code complexity.reduced}
	 * and our {@code Verifier.additionCount}.</p>
	 */
	public static SchemeSlp fromFreshVars(JsonNode root) {
		JsonNode nArr = root.get("n");
		int n1 = nArr.get(0).asInt(), n2 = nArr.get(1).asInt(), n3 = nArr.get(2).asInt();
		int rank = root.get("m").asInt();
		int dimU = n1 * n2, dimV = n2 * n3, dimW = n1 * n3;
		Slp a = side(root.get("u_fresh"), root.get("u"), dimU, rank, null);
		Slp b = side(root.get("v_fresh"), root.get("v"), dimV, rank, null);
		int[] rowMap = new int[dimW];
		for (int pDron = 0; pDron < dimW; pDron++) {
			int i = pDron % n1, j = pDron / n1;
			rowMap[pDron] = i * n3 + j;
		}
		Slp out = side(root.get("w_fresh"), root.get("w"), rank, dimW, rowMap);
		return new SchemeSlp(a, b, out);
	}

	private static Slp side(JsonNode fresh, JsonNode rows, int nInputs, int nForms, int[] rowMap) {
		if (rows == null || !rows.isArray() || rows.size() != nForms) {
			throw new IllegalArgumentException(
					"expected " + nForms + " rows, got " + (rows == null ? "none" : rows.size()));
		}
		Builder bld = new Builder(nInputs);
		int nFresh = fresh == null ? 0 : fresh.size();
		int[] slotOfFresh = new int[nFresh];
		for (int k = 0; k < nFresh; k++) {
			slotOfFresh[k] = bld.sum(terms(fresh.get(k), nInputs, slotOfFresh, k));
		}
		int[] formResult = new int[nForms];
		for (int f = 0; f < nForms; f++) {
			List<long[]> terms = terms(rows.get(f), nInputs, slotOfFresh, nFresh);
			int target = rowMap == null ? f : rowMap[f];
			formResult[target] = terms.isEmpty() ? -1 : bld.sum(terms);
		}
		return bld.build(formResult);
	}

	/** {@code {slot, coefficient-bits}} pairs of one row; a fresh reference must point below {@code limit}. */
	private static List<long[]> terms(JsonNode row, int nInputs, int[] slotOfFresh, int limit) {
		List<long[]> out = new ArrayList<>();
		for (JsonNode obj : row) {
			int idx = obj.get("index").asInt();
			double val = obj.get("value").asDouble();
			if (val == 0.0) {
				continue;
			}
			int slot;
			if (idx < nInputs) {
				slot = idx;
			} else {
				int k = idx - nInputs;
				if (k >= limit) {
					throw new IllegalArgumentException("fresh reference " + idx + " is not defined yet");
				}
				slot = slotOfFresh[k];
			}
			out.add(new long[] { slot, Double.doubleToLongBits(val) });
		}
		return out;
	}

	/** Emits SSA ops; slots are handed out in op order, as {@link Slp} requires. */
	private static final class Builder {
		private final int nInputs;
		private final List<Slp.Op> ops = new ArrayList<>();
		private int adds, muls;

		Builder(int nInputs) {
			this.nInputs = nInputs;
		}

		private int emit(char kind, int x, int y, double c) {
			int target = nInputs + ops.size();
			ops.add(new Slp.Op(kind, target, x, y, c));
			if (kind == '*') {
				muls++;
			} else {
				adds++;
			}
			return target;
		}

		/** The slot holding Σ coef·slot over {@code terms}; one addition per extra term. */
		int sum(List<long[]> terms) {
			List<Integer> plus = new ArrayList<>(), minus = new ArrayList<>();
			List<long[]> other = new ArrayList<>();
			for (long[] t : terms) {
				double c = Double.longBitsToDouble(t[1]);
				if (c == 1.0) {
					plus.add((int) t[0]);
				} else if (c == -1.0) {
					minus.add((int) t[0]);
				} else {
					other.add(t);
				}
			}
			int acc;
			boolean negated = false;
			if (!plus.isEmpty()) {
				acc = plus.get(0);
				for (int i = 1; i < plus.size(); i++) {
					acc = emit('+', acc, plus.get(i), 0.0);
				}
				for (int m : minus) {
					acc = emit('-', acc, m, 0.0);
				}
			} else if (!minus.isEmpty()) {
				// A sum of negated terms: accumulate the positive sum, negate once at the end.
				acc = minus.get(0);
				for (int i = 1; i < minus.size(); i++) {
					acc = emit('+', acc, minus.get(i), 0.0);
				}
				negated = true;
			} else {
				long[] t0 = other.remove(0);
				acc = emit('*', (int) t0[0], -1, Double.longBitsToDouble(t0[1]));
			}
			for (long[] t : other) {
				int scaled = emit('*', (int) t[0], -1, Double.longBitsToDouble(t[1]));
				acc = negated ? emit('-', scaled, acc, 0.0) : emit('+', acc, scaled, 0.0);
				negated = false;
			}
			if (negated) {
				acc = emit('*', acc, -1, -1.0);
			}
			return acc;
		}

		Slp build(int[] formResult) {
			return new Slp(nInputs, List.copyOf(ops), formResult, adds, muls);
		}
	}
}
