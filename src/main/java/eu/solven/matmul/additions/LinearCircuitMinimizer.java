package eu.solven.matmul.additions;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Heuristic minimiser for a <em>linear straight-line program</em> (LSP) that
 * evaluates a set of linear forms over a shared variable set — the additive
 * structure behind a bilinear matmul scheme. Produces a complete, replayable
 * {@link Slp} (not just a count), so the construction can be reconstructed.
 *
 * <p>Pipeline per call:</p>
 * <ol>
 *   <li><b>Scalar pre-pass</b> — a term with coefficient other than ±1 (Waksman
 *       ±½, Smirnov ±2, …) emits one shared {@code '*'} (scalar-mul) op
 *       {@code s = c·v}, reused across forms; {@code s} then enters the CSE as a
 *       +1 atom.</li>
 *   <li><b>Common-subexpression elimination</b> (Paar 1997, signed): repeatedly
 *       extract the signed atom-pair co-occurring in the most forms as one shared
 *       temporary (one {@code +}/{@code -} op), until no pair recurs. A pair and
 *       its full negation share the temporary (the second use subtracts it).</li>
 *   <li><b>Assembly</b> — each form folds its remaining atoms into a single
 *       result slot; {@code formResult[f]} records it.</li>
 * </ol>
 *
 * <p>Finding the true minimum is NP-hard; this greedy is a strong, fast upper
 * bound (it reproduces Winograd ⟨2,2,2⟩'s canonical 15). The emitted SLP is
 * always exact — {@link Slp#reconstructs} verifies it against the input forms.</p>
 */
public final class LinearCircuitMinimizer {

	private LinearCircuitMinimizer() {}

	/** @param slp the replayable program; @param naiveAdditions the
	 *  independent-per-form baseline {@code Σ(nnz−1)} (= Verifier.additionCount). */
	public record Result(Slp slp, int naiveAdditions) {
		public int additions() {
			return slp.additions();
		}

		public int scalarMults() {
			return slp.scalarMults();
		}

		public int saved() {
			return naiveAdditions - slp.additions();
		}
	}

	/**
	 * Minimise the additions to compute every row of {@code forms}
	 * ({@code forms[i][v]} = coefficient of variable {@code v} in form {@code i}).
	 */
	public static Result minimize(double[][] forms) {
		int nForms = forms.length;
		int nVars = nForms == 0 ? 0 : forms[0].length;

		List<Slp.Op> ops = new ArrayList<>();
		int[] nextSlot = { nVars };
		int[] scalarMults = { 0 };
		Map<Long, Integer> scalarSlot = new HashMap<>(); // (var,coeff) -> slot

		// live[i] : atom-slot -> sign (±1). Atoms are SSA slots (base var, scalar
		// temp, or CSE temp). naive = Σ(nnz−1) over the original forms.
		List<Map<Integer, Integer>> live = new ArrayList<>(nForms);
		int naive = 0;
		for (int i = 0; i < nForms; i++) {
			Map<Integer, Integer> f = new HashMap<>();
			int nnz = 0;
			for (int v = 0; v < nVars; v++) {
				double c = forms[i][v];
				if (c == 0.0) {
					continue;
				}
				nnz++;
				if (c == 1.0) {
					f.merge(v, 1, Integer::sum);
				} else if (c == -1.0) {
					f.merge(v, -1, Integer::sum);
				} else {
					long sk = scalarKey(v, c);
					Integer slot = scalarSlot.get(sk);
					if (slot == null) {
						slot = nextSlot[0]++;
						ops.add(new Slp.Op('*', slot, v, -1, c));
						scalarMults[0]++;
						scalarSlot.put(sk, slot);
					}
					f.merge(slot, 1, Integer::sum);
				}
			}
			// merge() can leave a 0 sign if a var appeared twice with opposite signs
			f.values().removeIf(s -> s == 0);
			live.add(f);
			if (nnz > 0) {
				naive += nnz - 1;
			}
		}

		int additions = 0;

		// CSE: greedily extract the most-frequent signed pair.
		while (true) {
			Map<Long, Integer> count = new HashMap<>();
			Map<Long, int[]> canonical = new HashMap<>();
			for (Map<Integer, Integer> f : live) {
				if (f.size() < 2) {
					continue;
				}
				Integer[] atoms = f.keySet().toArray(new Integer[0]);
				for (int x = 0; x < atoms.length; x++) {
					for (int y = x + 1; y < atoms.length; y++) {
						int a = atoms[x], b = atoms[y];
						long key = pairKey(a, f.get(a), b, f.get(b));
						count.merge(key, 1, Integer::sum);
						canonical.putIfAbsent(key, canonPair(a, f.get(a), b, f.get(b)));
					}
				}
			}
			long bestKey = -1;
			int bestCount = 1;
			for (Map.Entry<Long, Integer> e : count.entrySet()) {
				if (e.getValue() > bestCount) {
					bestCount = e.getValue();
					bestKey = e.getKey();
				}
			}
			if (bestKey < 0) {
				break;
			}
			int[] cp = canonical.get(bestKey); // {a, sa=+1, b, sb}
			int a = cp[0], b = cp[2], sb = cp[3];
			int temp = nextSlot[0]++;
			ops.add(new Slp.Op(sb > 0 ? '+' : '-', temp, a, b, 0));
			additions++;
			for (Map<Integer, Integer> f : live) {
				Integer fa = f.get(a), fb = f.get(b);
				if (fa == null || fb == null) {
					continue;
				}
				if (fa == 1 && fb == sb) {
					f.remove(a);
					f.remove(b);
					f.put(temp, 1);
				} else if (fa == -1 && fb == -sb) {
					f.remove(a);
					f.remove(b);
					f.put(temp, -1);
				}
			}
		}

		// Assembly: fold each form's remaining atoms into one result slot.
		int[] formResult = new int[nForms];
		for (int i = 0; i < nForms; i++) {
			Map<Integer, Integer> f = live.get(i);
			if (f.isEmpty()) {
				formResult[i] = -1; // zero form
				continue;
			}
			List<Integer> pos = new ArrayList<>();
			List<Integer> neg = new ArrayList<>();
			for (Map.Entry<Integer, Integer> e : f.entrySet()) {
				(e.getValue() > 0 ? pos : neg).add(e.getKey());
			}
			int acc;
			List<Integer> rest;
			boolean negateAtEnd = false;
			if (!pos.isEmpty()) {
				acc = pos.get(0);
				for (int k = 1; k < pos.size(); k++) {
					int t = nextSlot[0]++;
					ops.add(new Slp.Op('+', t, acc, pos.get(k), 0));
					additions++;
					acc = t;
				}
				rest = neg;
				for (int nb : rest) {
					int t = nextSlot[0]++;
					ops.add(new Slp.Op('-', t, acc, nb, 0));
					additions++;
					acc = t;
				}
			} else {
				// all-negative form: build the positive sum, then negate once.
				acc = neg.get(0);
				for (int k = 1; k < neg.size(); k++) {
					int t = nextSlot[0]++;
					ops.add(new Slp.Op('+', t, acc, neg.get(k), 0));
					additions++;
					acc = t;
				}
				negateAtEnd = true;
			}
			if (negateAtEnd) {
				int t = nextSlot[0]++;
				ops.add(new Slp.Op('*', t, acc, -1, -1.0));
				scalarMults[0]++;
				acc = t;
			}
			formResult[i] = acc;
		}

		Slp slp = new Slp(nVars, ops, formResult, additions, scalarMults[0]);
		return new Result(slp, naive);
	}

	private static long scalarKey(int var, double coeff) {
		return (((long) var) << 32) ^ Double.hashCode(coeff);
	}

	/** Negation-normalised key for the unordered signed pair {(a,sa),(b,sb)}. */
	private static long pairKey(int a, int sa, int b, int sb) {
		int[] c = canonPair(a, sa, b, sb);
		return (((long) c[0]) << 34) ^ (((long) c[2]) << 2) ^ ((c[3] > 0) ? 1L : 0L);
	}

	/** Canonical orientation: atoms ascending, then flip both signs so the first
	 *  atom's sign is +1 (so a pair and its negation coincide). */
	private static int[] canonPair(int a, int sa, int b, int sb) {
		if (a > b) {
			int t = a; a = b; b = t;
			t = sa; sa = sb; sb = t;
		}
		if (sa < 0) {
			sa = -sa;
			sb = -sb;
		}
		return new int[] { a, sa, b, sb };
	}
}
