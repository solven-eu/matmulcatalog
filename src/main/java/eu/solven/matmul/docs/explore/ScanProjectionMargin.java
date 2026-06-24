package eu.solven.matmul.docs.explore;

import eu.solven.matmul.recombination.Recombination;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Stream;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.ProjectionSearch;
import eu.solven.matmul.catalog.SchemeIO;
import lombok.extern.slf4j.Slf4j;

/**
 * Scan small committed schemes, compute the projection margin μ (death rate) of
 * each, and report how well each projects one index down — for the paper's
 * worked example of the (rank, projection-margin) trade-off. Output columns:
 * shape, source, rank R, μ, projected rank {@code R−μ}, and the current
 * catalog-best rank at the projected shape (so "acceptable, not optimal" is
 * visible). Sorted by μ descending.
 *
 * <pre>{@code
 * mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.explore.ScanProjectionMargin -Dexec.args="8"
 * }</pre>
 */
@Slf4j
public final class ScanProjectionMargin {
	private ScanProjectionMargin() {}

	private record Row(String shape, String source, int n, int m, int p, int r,
			int mu, char axis, int projShape0, int projShape1, int projShape2, int catalogBestProj) {}

	public static void main(String[] args) throws Exception {
		int maxDim = args.length > 0 ? Integer.parseInt(args[0]) : 8;
		FieldAwareLookup lookup = new FieldAwareLookup(Field.Q);
		Path root = Path.of("src/main/resources/schemes");
		List<Row> rows = new ArrayList<>();
		try (Stream<Path> walk = Files.walk(root)) {
			for (Path p : (Iterable<Path>) walk::iterator) {
				String name = p.getFileName().toString();
				if (!name.endsWith(".json")) continue;
				NonCubicBilinearAlgorithm a;
				try {
					var root2 = SchemeIO.parseJson(p.toFile());
					if (SchemeIO.isStub(root2) || SchemeIO.isNonBilinear(root2) || SchemeIO.isComplex(root2)) {
						continue;  // need explicit real matrices
					}
					a = SchemeIO.read(root2);
				} catch (Exception e) {
					continue;
				}
				if (Math.max(a.n, Math.max(a.m, a.p)) > maxDim) continue;
				if (a.n < 2 && a.m < 2 && a.p < 2) continue;

				// μ per axis, plus which axis achieves it (so we know which dim drops).
				int muN = axisOnly(a, 'n'), muM = axisOnly(a, 'm'), muP = axisOnly(a, 'p');
				int mu = Math.max(muN, Math.max(muM, muP));
				char axis = muN >= muM && muN >= muP ? 'n' : (muM >= muP ? 'm' : 'p');
				if (mu == 0) continue;  // projecting any index can't help

				int pn = a.n, pm = a.m, pp = a.p;
				if (axis == 'n') pn--; else if (axis == 'm') pm--; else pp--;
				if (pn < 1 || pm < 1 || pp < 1) continue;
				int best = lookup.findRank(pn, pm, pp);
				rows.add(new Row(a.n + "x" + a.m + "x" + a.p, sourceOf(name), a.n, a.m, a.p, a.r,
						mu, axis, pn, pm, pp, best >= Recombination.SotaResolver.UNKNOWN_RANK ? -1 : best));
			}
		}
		rows.sort(Comparator.comparingInt((Row x) -> x.mu).reversed());
		log.info("scheme                         R    μ  drop  ->  ⟨proj⟩    R−μ   catalogBest(proj)  verdict");
		for (Row x : rows) {
			int projRank = x.r - x.mu;
			String verdict = x.catalogBestProj < 0 ? "(no entry — fills gap)"
					: projRank < x.catalogBestProj ? "BEATS catalog"
					: projRank == x.catalogBestProj ? "= catalog best" : "above best by " + (projRank - x.catalogBestProj);
			log.info(String.format("%-22s %-8s %4d %3d   %c   -> ⟨%d,%d,%d⟩   %4d      %5s         %s",
					x.shape, x.source, x.r, x.mu, x.axis, x.projShape0, x.projShape1, x.projShape2,
					projRank, x.catalogBestProj < 0 ? "—" : Integer.toString(x.catalogBestProj), verdict));
		}
	}

	/** μ restricted to one axis, by zeroing the other two axes' contribution: we
	 *  reuse the public μ over a single-axis view by temporarily reporting only
	 *  that axis. Simpler: recompute via projectionMargin on a shape where the
	 *  other axes are dimension-1 is wrong; instead expose per-axis here. */
	private static int axisOnly(NonCubicBilinearAlgorithm a, char axis) {
		// projectionMargin already takes the max over axes; to get a per-axis value
		// we evaluate every single-index drop on that axis and take R − survivors.
		int dim = axis == 'n' ? a.n : axis == 'm' ? a.m : a.p;
		int best = 0;
		for (int i = 0; i < dim; i++) {
			int[] keepN = keepAllBut(a.n, axis == 'n' ? i : -1);
			int[] keepM = keepAllBut(a.m, axis == 'm' ? i : -1);
			int[] keepP = keepAllBut(a.p, axis == 'p' ? i : -1);
			NonCubicBilinearAlgorithm proj = eu.solven.matmul.catalog.Compose.project(a, keepN, keepM, keepP);
			best = Math.max(best, a.r - proj.r);
		}
		return best;
	}

	private static int[] keepAllBut(int dim, int drop) {
		int[] out = new int[drop < 0 ? dim : dim - 1];
		int w = 0;
		for (int i = 0; i < dim; i++) if (i != drop) out[w++] = i;
		return out;
	}

	private static String sourceOf(String name) {
		int dash = name.indexOf('-');
		return dash > 0 ? name.substring(0, Math.min(dash, 8)) : name;
	}
}
