package eu.solven.matmul.research;

import eu.solven.matmul.NonBilinearAlgorithm;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.papers.makarov1986.Makarov22;
import eu.solven.matmul.papers.makarov1986.Makarov22.WSpec;

/**
 * Single-flip search for the typo in Makarov 1986 ⟨3,3,3⟩=22.
 * Tries flipping the sign of one γ_k in one c_{i,j}, checks residual.
 * Reports any configuration with residual = 0.
 *
 * <p>Also tries: removing one γ_k from one c_{i,j} (without re-adding).
 * Skip: γ_k swapped between two c_{i,j}'s (combinatorially larger).</p>
 */
public final class MakarovSearch {
	private MakarovSearch() {}

	public static void main(String[] args) {
		WSpec[][] base = Makarov22.defaultOutputs();
		double defaultRes = Verifier.residualNonBilinear(Makarov22.build(base));
		System.out.printf("Baseline residual: %.4g%n", defaultRes);
		System.out.println();
		System.out.println("Single-sign-flip search (flip sign of one γ_k in one c_{i,j}):");
		int hits = 0;
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				for (var t : base[i][j].terms()) {
					int g = t.gammaIdx();
					WSpec flipped = base[i][j].flipSign(g);
					WSpec[][] mod = deepCopy(base);
					mod[i][j] = flipped;
					double r = Verifier.residualNonBilinear(Makarov22.build(mod));
					if (r < 1e-9) {
						System.out.printf("  HIT: flip sign of γ%d in c%d%d → residual %.4g%n",
								g, i + 1, j + 1, r);
						hits++;
					}
				}
			}
		}

		System.out.println();
		System.out.println("Single-removal search (remove one γ_k from one c_{i,j}):");
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				for (var t : base[i][j].terms()) {
					int g = t.gammaIdx();
					var newTerms = new java.util.ArrayList<>(base[i][j].terms());
					newTerms.removeIf(x -> x.gammaIdx() == g);
					WSpec[][] mod = deepCopy(base);
					mod[i][j] = new WSpec(newTerms);
					double r = Verifier.residualNonBilinear(Makarov22.build(mod));
					if (r < 1e-9) {
						System.out.printf("  HIT: drop γ%d from c%d%d → residual %.4g%n",
								g, i + 1, j + 1, r);
						hits++;
					}
				}
			}
		}

		System.out.println();
		System.out.println("Single-addition search (add ±γ_k to one c_{i,j}, k not currently present):");
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				java.util.Set<Integer> present = new java.util.HashSet<>();
				for (var t : base[i][j].terms()) present.add(t.gammaIdx());
				for (int g = 1; g <= 22; g++) {
					if (present.contains(g)) continue;
					for (int sign : new int[]{ +1, -1 }) {
						var newTerms = new java.util.ArrayList<>(base[i][j].terms());
						newTerms.add(new WSpec.Term(sign, g));
						WSpec[][] mod = deepCopy(base);
						mod[i][j] = new WSpec(newTerms);
						double r = Verifier.residualNonBilinear(Makarov22.build(mod));
						if (r < 1e-9) {
							System.out.printf("  HIT: add %sγ%d to c%d%d → residual %.4g%n",
									sign > 0 ? "+" : "-", g, i + 1, j + 1, r);
							hits++;
						}
					}
				}
			}
		}

		System.out.println();
		System.out.printf("Total single-edit hits: %d%n", hits);
		if (hits == 0) {
			System.out.println("No single-edit fix found. The typo (if any) involves >=2 edits, or");
			System.out.println("lies inside the γ_k definitions themselves (not the W output combos).");
		}
	}

	private static WSpec[][] deepCopy(WSpec[][] src) {
		WSpec[][] out = new WSpec[3][3];
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				out[i][j] = new WSpec(new java.util.ArrayList<>(src[i][j].terms()));
		return out;
	}
}
