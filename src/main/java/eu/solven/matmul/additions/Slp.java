package eu.solven.matmul.additions;

import java.util.List;

/**
 * A linear straight-line program: the full, replayable sequence of operations
 * that evaluates a set of linear forms from input variables — registered so a
 * scheme's additive structure can be <em>reconstructed</em>, not just counted.
 *
 * <p>Value model: a flat value array {@code val[]} of length
 * {@code nInputs + ops.size()}. Slots {@code 0 .. nInputs-1} are the input
 * variables; op {@code i} writes slot {@code nInputs + i}. Each {@link Op} reads
 * earlier slots only (SSA / topologically ordered), so {@link #evaluate} is a
 * single left-to-right pass. {@code formResult[f]} is the slot holding form
 * {@code f}'s value (or {@code -1} for the zero form).</p>
 *
 * <p>Op kinds: {@code '+'} → {@code val[target]=val[x]+val[y]};
 * {@code '-'} → {@code val[x]-val[y]}; {@code '*'} → scalar
 * {@code c*val[x]} (used for coefficients other than ±1). Additions count
 * {@code +}/{@code -}; scalar-mults count {@code *}.</p>
 */
public record Slp(int nInputs, List<Op> ops, int[] formResult, int additions, int scalarMults) {

	/** One SSA instruction. For {@code '*'}, {@code y} is unused and {@code c} is
	 *  the scalar; for {@code '+'}/{@code '-'}, {@code c} is unused. */
	public record Op(char kind, int target, int x, int y, double c) {}

	/** Evaluate every form for the given input vector. */
	public double[] evaluate(double[] inputs) {
		double[] val = new double[nInputs + ops.size()];
		System.arraycopy(inputs, 0, val, 0, nInputs);
		for (Op o : ops) {
			val[o.target()] = switch (o.kind()) {
				case '+' -> val[o.x()] + val[o.y()];
				case '-' -> val[o.x()] - val[o.y()];
				case '*' -> o.c() * val[o.x()];
				default -> throw new IllegalStateException("bad op " + o.kind());
			};
		}
		double[] out = new double[formResult.length];
		for (int f = 0; f < formResult.length; f++) {
			out[f] = formResult[f] < 0 ? 0.0 : val[formResult[f]];
		}
		return out;
	}

	/** Human-readable program, one op per line (inputs named {@code varName0…}). */
	public java.util.List<String> render(String varName) {
		java.util.List<String> out = new java.util.ArrayList<>();
		for (Op o : ops) {
			String t = "t" + o.target();
			if (o.kind() == '*') {
				out.add(t + " = " + o.c() + "*" + slot(o.x(), varName));
			} else {
				out.add(t + " = " + slot(o.x(), varName) + " " + o.kind() + " " + slot(o.y(), varName));
			}
		}
		return out;
	}

	private String slot(int id, String varName) {
		return id < nInputs ? varName + id : "t" + id;
	}

	/**
	 * Self-check: confirm the SLP reproduces {@code forms} exactly (each form is a
	 * coefficient row over the inputs). Tests with the standard basis are enough —
	 * the program is linear, so agreement on the basis ⇒ agreement everywhere.
	 */
	public boolean reconstructs(double[][] forms) {
		if (forms.length != formResult.length) {
			return false;
		}
		for (int v = 0; v < nInputs; v++) {
			double[] e = new double[nInputs];
			e[v] = 1.0;
			double[] got = evaluate(e);
			for (int f = 0; f < forms.length; f++) {
				if (Math.abs(got[f] - forms[f][v]) > 1e-9) {
					return false;
				}
			}
		}
		return true;
	}
}
