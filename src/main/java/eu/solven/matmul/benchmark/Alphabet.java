package eu.solven.matmul.benchmark;

/**
 * Coefficient alphabet for a decomposition. Determines which encoder the
 * {@link BenchmarkRunner} can use.
 *
 * <ul>
 *   <li>{@link #Z2} — GF(2): scalars in {0, 1}, addition is XOR. Encoded by
 *       {@code Z2CnfEncoder} (implemented).</li>
 *   <li>{@link #Z3} — GF(3): scalars in {0, 1, 2}; mod-3 arithmetic. No encoder
 *       yet; rows are recorded as {@code NO_ENCODER} for now.</li>
 *   <li>{@link #Z} — integers with restricted ranges (typically
 *       {-1, 0, +1}). No encoder yet.</li>
 * </ul>
 */
public enum Alphabet {
	Z2, Z3, Z;
}
