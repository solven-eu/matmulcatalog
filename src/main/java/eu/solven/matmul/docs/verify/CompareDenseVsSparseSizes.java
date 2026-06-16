package eu.solven.matmul.docs.verify;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import java.util.TreeMap;

import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.catalog.MatrixJsonFormatter;
import eu.solven.matmul.catalog.SchemeIO;

/**
 * Measures, per cubic dim, the on-disk JSON byte size of the dense
 * {@code u/v/w} encoding versus the sparse {@code *_sparse} encoding, to
 * inform the {@code SPARSE_DIM_THRESHOLD} sparse toggle (currently 16; the
 * question is whether 8 is better).
 *
 * <p>For every cubic scheme {@code ⟨k,k,k⟩} whose factor matrices are present
 * on disk (sections 2..15 are dense; atoms at 16 too), we load it, re-serialise
 * it BOTH ways through the canonical {@code MatrixJsonFormatter}, and compare
 * byte counts. The crossover dim (where sparse first becomes smaller than
 * dense) is the natural value for the toggle.</p>
 */
public final class CompareDenseVsSparseSizes {

	private CompareDenseVsSparseSizes() {}

	private static final int MAX_PER_DIM = 40;

	public static void main(String[] args) {
		File root = new File("src/main/resources/schemes");
		// dim -> [count, sumDense, sumSparse, sumRank, sparseSmallerCount]
		TreeMap<Integer, long[]> byDim = new TreeMap<>();

		// sectionN dirs live under known/derived/curated since the 2026-06-09 split —
		// recurse to find them at any depth. byDim aggregates by maxDim, so processing
		// e.g. known/section8 + derived/section8 separately is fine (both add to dim 8).
		File[] sections;
		try (java.util.stream.Stream<java.nio.file.Path> w = java.nio.file.Files.walk(root.toPath())) {
			sections = w.filter(java.nio.file.Files::isDirectory)
					.filter(p -> p.getFileName().toString().startsWith("section"))
					.map(java.nio.file.Path::toFile).toArray(File[]::new);
		} catch (java.io.IOException e) { sections = null; }
		if (sections == null || sections.length == 0) {
			System.err.println("no section dirs under " + root.getAbsolutePath());
			return;
		}
		for (File sec : sections) {
			File[] files = sec.listFiles((d, n) -> n.endsWith(".json"));
			if (files == null) continue;
			List<File> cubic = new ArrayList<>();
			for (File f : files) {
				cubic.add(f);
			}
			int sampled = 0;
			for (File f : cubic) {
				if (sampled >= MAX_PER_DIM) break;
				NonCubicBilinearAlgorithm alg;
				try {
					alg = SchemeIO.read(f);
				} catch (Exception e) {
					continue; // stub / non-bilinear / unreadable — skip
				}
				if (alg == null) continue;
				if (alg.n != alg.m || alg.m != alg.p) continue; // cubic only
				int k = alg.n;
				int dense, sparse;
				try {
					dense = MatrixJsonFormatter.format(SchemeIO.toJson(alg, null, null)).length();
					sparse = MatrixJsonFormatter.format(SchemeIO.toJsonSparse(alg)).length();
				} catch (Exception e) {
					continue;
				}
				long[] acc = byDim.computeIfAbsent(k, kk -> new long[5]);
				acc[0]++;
				acc[1] += dense;
				acc[2] += sparse;
				acc[3] += alg.r;
				if (sparse < dense) acc[4]++;
				sampled++;
			}
		}

		System.out.printf("%nDense vs sparse JSON size for cubic schemes (sampled ≤%d per dim)%n", MAX_PER_DIM);
		System.out.printf("%-4s %6s %8s %12s %12s %8s %10s%n",
				"dim", "n", "avgRank", "denseBytes", "sparseBytes", "sp/dn", "sparseWins");
		for (var e : byDim.entrySet()) {
			long[] a = e.getValue();
			if (a[0] == 0) continue;
			double dn = (double) a[1] / a[0];
			double sp = (double) a[2] / a[0];
			System.out.printf(Locale.US, "%-4d %6d %8.0f %12.0f %12.0f %8.2f %6d/%-3d%n",
					e.getKey(), a[0], (double) a[3] / a[0], dn, sp, sp / dn, a[4], a[0]);
		}
	}
}
