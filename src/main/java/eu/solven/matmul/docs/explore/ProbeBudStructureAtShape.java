package eu.solven.matmul.docs.explore;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Stream;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.algebra.Field;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.catalog.SerendipitousBudProduct;
import eu.solven.matmul.search.flip.FlipObjectives;
import eu.solven.matmul.search.flip.FlipScheme;
import lombok.extern.slf4j.Slf4j;

/**
 * Throwaway probe: report the bud-structure of every catalog scheme at a given
 * shape — budScore (Σ independent class sizes ≥2 over U/V/W), the per-axis
 * independent class-size multisets, the greedy disjoint bud decomposition,
 * projection margin, and the per-axis self-serendipity savings σ (priced at
 * inner = own shape; a greedy bound like every bud figure).
 *
 * <pre>
 *   mvn -q -ntp exec:java \
 *       -Dexec.mainClass=eu.solven.matmul.docs.explore.ProbeBudStructureAtShape \
 *       -Dexec.args="--shape=2x9x10"
 * </pre>
 */
@Slf4j
public final class ProbeBudStructureAtShape {

	private ProbeBudStructureAtShape() {}

	public static void main(String[] args) throws IOException {
		String shape = Arrays.stream(args)
				.filter(a -> a.startsWith("--shape="))
				.map(a -> a.substring("--shape=".length()))
				.findFirst()
				.orElseThrow(() -> new IllegalArgumentException("--shape=NxMxP is required"));
		String[] d = shape.split("x");
		int n = Integer.parseInt(d[0]);
		int m = Integer.parseInt(d[1]);
		int p = Integer.parseInt(d[2]);

		// Filename prefilter only (fast); shape is then confirmed from CONTENT.
		List<File> files = new ArrayList<>();
		try (Stream<Path> walk = Files.walk(Path.of("src/main/resources/schemes"))) {
			walk.filter(f -> f.getFileName().toString().startsWith(shape + "-"))
					.filter(f -> f.toString().endsWith(".json"))
					.forEach(f -> files.add(f.toFile()));
		}
		log.info("{} candidate files at ⟨{},{},{}⟩", files.size(), n, m, p);

		FieldAwareLookup qLookup = new FieldAwareLookup(Field.Q);
		for (File f : files) {
			NonCubicBilinearAlgorithm a = SchemeIO.readBilinear(f);
			if (a.n != n || a.m != m || a.p != p) {
				log.warn("{}: filename label disagrees with content shape ⟨{},{},{}⟩ — skipped",
						f.getName(), a.n, a.m, a.p);
				continue;
			}
			FlipScheme s = FlipScheme.of(a);
			int[][] classes = SerendipitousBudProduct.independentClassSizes(a);
			SerendipitousBudProduct.BudDecomposition greedy = SerendipitousBudProduct.findBuds(a);
			long[] sigma = FlipObjectives.serendipitySavingByAxis(s, qLookup, n, m, p);
			log.info("{}", f.getName());
			log.info("  rank={} budScore={} margin={} σ(U,V,W)=({},{},{})",
					a.r, FlipObjectives.budScore(s), FlipObjectives.projectionMargin(s),
					sigma[0], sigma[1], sigma[2]);
			log.info("  independent classes ≥2: U={} V={} W={}",
					nonTrivial(classes[0]), nonTrivial(classes[1]), nonTrivial(classes[2]));
			log.info("  greedy disjoint buds: {} (trivial={})",
					greedy.buds().stream()
							.map(b -> b.type() + "[" + b.terms().length + "]")
							.toList(),
					greedy.trivial().length);
		}
	}

	private static List<Integer> nonTrivial(int[] sizes) {
		return Arrays.stream(sizes).filter(x -> x >= 2).boxed().toList();
	}
}
