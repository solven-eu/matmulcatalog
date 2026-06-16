package eu.solven.matmul.docs.explore;

import java.util.List;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.Verifier;
import eu.solven.matmul.catalog.ProjectionSearch;
import eu.solven.matmul.papers.dis2009.PanTrilinearAggregation;
import lombok.extern.slf4j.Slf4j;

/**
 * One-off probe: can a STRUCTURED Pan-trilinear-aggregation cube ⟨n,n,n⟩ project
 * to ⟨n-1,n,n⟩ as well as (or better than) FMM? FMM gets ⟨29,30,30⟩=12588 by
 * projecting a ⟨30,30,30⟩=12710 Pan-TA cube — yet our rank-best cube is the
 * flat SZ ⟨30,30,30⟩=12688, which only projects to 12626. This tests the
 * hypothesis that the higher-rank STRUCTURED cube projects to a lower child
 * (the projection analogue of "bud-rich base beats rank-minimal base").
 *
 * <pre>{@code
 * mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.explore.ProbePanTAProjection -Dexec.args="30"
 * }</pre>
 */
@Slf4j
public final class ProbePanTAProjection {
	private ProbePanTAProjection() {}

	public static void main(String[] args) {
		int n = args.length > 0 ? Integer.parseInt(args[0]) : 30;
		long t0 = System.nanoTime();
		NonCubicBilinearAlgorithm cube = PanTrilinearAggregation.build(n);
		log.info("PanTA ⟨{},{},{}⟩ built: r={}  ({} ms)  exact={}  projectionMargin μ={}",
				n, n, n, cube.r, (System.nanoTime() - t0) / 1_000_000L,
				Verifier.passesRandomMatmulSpotCheck(cube),
				ProjectionSearch.projectionMargin(cube));

		// Project one axis n→n-1 (FMM's [[1,0],[0]]). bestFor exhaustively tries every
		// single-index drop and keeps the survivor-minimal one.
		long t1 = System.nanoTime();
		var hit = ProjectionSearch.bestFor(n - 1, n, n, List.of(cube), Long.MAX_VALUE, 1);
		long ms = (System.nanoTime() - t1) / 1_000_000L;
		if (hit.isEmpty()) {
			log.info("no projection found ({} ms)", ms);
			return;
		}
		var h = hit.get();
		log.info("PROJECTED ⟨{},{},{}⟩ = {}  (drop best single index; {} ms)  exact={}",
				n - 1, n, n, h.rank(), ms, Verifier.passesRandomMatmulSpotCheck(h.scheme()));
		log.info("Compare: FMM ⟨29,30,30⟩=12588 ; our flat-SZ projection = 12626 ; PanTA-projection = {}",
				h.rank());
	}
}
