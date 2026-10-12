package eu.solven.matmul.io;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.assertThatThrownBy;

import org.junit.jupiter.api.Test;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.verifiers.Verifier;

/**
 * FMM-Lille's tensor artifacts ({@code *_tensor.mpl}, Triad listing). The trap is the
 * third matrix: it is {@code p×n} (the file verifies with {@code Transpose(T[i][3])}),
 * so a parser that reads it as {@code n×p} yields a scheme of the right shape and rank
 * that is simply wrong — invisible on square or symmetric examples.
 */
public class TestMapleTensorParser {

	/** Render a scheme the way FMM does: {@code Triad([U (n×m), V (m×p), W (p×n)])}. */
	private static String triads(NonCubicBilinearAlgorithm a) {
		double[][] u = a.denseU(), v = a.denseV(), w = a.denseW();
		StringBuilder sb = new StringBuilder("A:=Matrix(" + a.n + ", " + a.m + ", [[A_1_1]]):\nTensor:=TriadSet([");
		for (int q = 0; q < a.r; q++) {
			if (q > 0) sb.append(", ");
			sb.append("Triad([");
			matrix(sb, a.n, a.m, q, (i, j) -> u[i * a.m + j]);
			sb.append(", ");
			matrix(sb, a.m, a.p, q, (j, k) -> v[j * a.p + k]);
			sb.append(", ");
			matrix(sb, a.p, a.n, q, (k, i) -> w[i * a.p + k]); // p×n: transposed
			sb.append("])");
		}
		return sb.append("]):\nmap(expand,A.B-add(x,i=1..").append(a.r).append("));\n").toString();
	}

	private interface Cell {
		double[] row(int r, int c);
	}

	private static void matrix(StringBuilder sb, int rows, int cols, int q, Cell cell) {
		sb.append("Matrix(").append(rows).append(", ").append(cols).append(", [");
		for (int r = 0; r < rows; r++) {
			sb.append(r > 0 ? "," : "").append('[');
			for (int c = 0; c < cols; c++) {
				double x = cell.row(r, c)[q];
				sb.append(c > 0 ? "," : "").append(x == 0.5 ? "1/2" : x == -0.5 ? "-1/2" : String.valueOf((long) x));
			}
			sb.append(']');
		}
		sb.append("])");
	}

	@Test
	public void a_non_square_scheme_round_trips_exactly() throws Exception {
		// ⟨2,3,4⟩: n ≠ p, so the p×n third matrix cannot be confused with n×p.
		NonCubicBilinearAlgorithm naive = NonCubicBilinearAlgorithm.naive(2, 3, 4);

		NonCubicBilinearAlgorithm parsed = MapleTensorParser.parse(triads(naive));

		assertThat(new int[] { parsed.n, parsed.m, parsed.p, parsed.r }).containsExactly(2, 3, 4, 24);
		assertThat(Verifier.isExactNonCubic(parsed)).isTrue();
	}

	@Test
	public void rational_entries_are_read_as_fractions() throws Exception {
		// Halve U and double V on one product: still exact, now with a 1/2 entry.
		NonCubicBilinearAlgorithm naive = NonCubicBilinearAlgorithm.naive(2, 3, 4);
		double[][] u = naive.denseU(), v = naive.denseV(), w = naive.denseW();
		for (double[] row : u) row[5] /= 2;
		for (double[] row : v) row[5] *= 2;
		NonCubicBilinearAlgorithm scaled = new NonCubicBilinearAlgorithm(2, 3, 4, u, v, w);

		NonCubicBilinearAlgorithm parsed = MapleTensorParser.parse(triads(scaled));

		assertThat(parsed.r).isEqualTo(24);
		assertThat(Verifier.isExactNonCubic(parsed)).isTrue();
		assertThat(java.util.Arrays.stream(parsed.denseU()).flatMapToDouble(java.util.Arrays::stream).anyMatch(x -> x == 0.5))
				.as("the 1/2 entry survived as a fraction").isTrue();
	}

	@Test
	public void a_placeholder_and_a_parametric_family_are_not_schemes() {
		String placeholder = "A:=Matrix(2, 2, [[A_1_1,A_1_2],[A_2_1,A_2_2]]):\nTensor:=Tensor:\n";
		assertThatThrownBy(() -> MapleTensorParser.parse(placeholder))
				.isInstanceOf(MapleTensorParser.UnsupportedArtifactException.class).hasMessageContaining("placeholder");

		String parametric = "Tensor:=TriadSet([Triad([Matrix(1, 1, [[1/l]]), Matrix(1, 1, [[l]]), Matrix(1, 1, [[1]])])]):";
		assertThatThrownBy(() -> MapleTensorParser.parse(parametric))
				.isInstanceOf(MapleTensorParser.UnsupportedArtifactException.class).hasMessageContaining("parametric");
	}
}
