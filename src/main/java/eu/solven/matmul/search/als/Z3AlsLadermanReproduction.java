package eu.solven.matmul.search.als;

import lombok.extern.slf4j.Slf4j;

import eu.solven.matmul.Verifier;

import eu.solven.matmul.BilinearAlgorithm;

/**
 * Smoke / sanity attempt at reproducing a Laderman-equivalent rank-23 algorithm
 * for ⟨3,3,3⟩ via random-restart Z/3-equivariant ALS.
 *
 * Targets {@link Verifier#trilinTensor(int)} — see {@link Z3Als#Z3Als(int, int, int)}
 * for why the trilinear form is the right tensor for raw Z/3 equivariance.
 * Any rank-23 solution found here converts to a {@link Verifier#matmulTensor}
 * decomposition by applying {@link Verifier#transposeW}.
 *
 * Structure: {@code 2 fixed + 7 orbits} (one of the published Z/3 partitions of
 * Laderman's algorithm). Other valid r=23 partitions: {@code 5 + 6, 8 + 5,
 * 11 + 4, ...} per [[RANK_3X3_SEARCH]] §4.1; if 2+7 plateaus, those are
 * straightforward to try by editing {@code f, g} below.
 *
 * Expected behaviour per the literature: a non-trivial fraction of random
 * inits converge at r=23. If after several dozen restarts no convergence is
 * seen, raise {@code maxIters}, then try other structures.
 *
 * Run with: {@code mvn -q compile && java -cp target/classes
 *           eu.solven.matmul.search.als.Z3AlsLadermanReproduction}
 */
@Slf4j
public class Z3AlsLadermanReproduction {

	public static void main(String[] args) {
		int n = 3;
		int f = 2;
		int g = 7;
		int r = f + 3 * g; // = 23
		int restarts = 30;
		int maxIters = 5000;
		double tol = 1e-10;
		long baseSeed = 42L;

		log.info(String.format("Z3Als reproduction attempt for ⟨3,3,3⟩ at r=%d%n", r));
		log.info(String.format("  structure  : %d fixed + %d orbits (Z/3-equivariant)%n", f, g));
		log.info(String.format("  target     : trilinTensor(3) — cyclic-symmetric trace(ABC)%n"));
		log.info(String.format("  restarts   : %d, maxIters: %d, tol: %.0e%n%n", restarts, maxIters, tol));

		Z3Als als = new Z3Als(n, f, g); // defaults to trilin
		long start = System.nanoTime();
		Z3Als.Result best = null;
		int firstConvergedSeed = -1;

		for (int s = 0; s < restarts; s++) {
			long t = System.nanoTime();
			Z3Als.Result result = als.fit(maxIters, tol, baseSeed + s);
			long ms = (System.nanoTime() - t) / 1_000_000;
			log.info(String.format("[seed %3d] residual=%.4e  iters=%5d  converged=%s  (%d ms)%n",
					baseSeed + s, result.residual, result.iterations, result.converged, ms));

			if (best == null || result.residual < best.residual) best = result;
			if (result.converged && firstConvergedSeed < 0) {
				firstConvergedSeed = (int) (baseSeed + s);
			}
		}

		long totalMs = (System.nanoTime() - start) / 1_000_000;
		log.info(String.format("%n=== DONE in %d ms ===%n", totalMs));
		log.info(String.format("Best residual across %d restarts: %.4e%n", restarts, best.residual));

		if (best.converged) {
			log.info(String.format("Converged at seed %d.%n", firstConvergedSeed));

			BilinearAlgorithm trilinAlg = best.algorithm;
			BilinearAlgorithm matmulAlg = Verifier.transposeW(trilinAlg);
			double resMatmul = Verifier.residual(matmulAlg);
			log.info(String.format("Converted to matmul-tensor decomposition; residual = %.4e%n", resMatmul));

			if (resMatmul < 1e-8) {
				log.info("");
				log.info("SUCCESS: this is a real-valued Laderman-equivalent rank-23 algorithm");
				log.info("for the ⟨3,3,3⟩ matmul tensor. Entries are floating-point and would");
				log.info("need rationalization (snap-to-rational + re-verify) to recover a");
				log.info("publishable closed-form algorithm.");
			} else {
				log.info("");
				log.info("WARNING: trilin residual is small but matmul residual is large.");
				log.info("Conversion may have a bug or the algorithm structure is degenerate.");
			}
		} else {
			log.info("Did not converge.");
			log.info("Next steps to try:");
			log.info("  - Increase restarts (current: " + restarts + ")");
			log.info("  - Increase maxIters (current: " + maxIters + ")");
			log.info("  - Try a different Z/3 structure: f ∈ {5, 8, 11, ...} with g = (23-f)/3");
		}
	}
}
