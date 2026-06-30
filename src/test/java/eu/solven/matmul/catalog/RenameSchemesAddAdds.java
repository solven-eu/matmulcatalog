package eu.solven.matmul.catalog;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import eu.solven.matmul.ComplexNonCubicBilinearAlgorithm;
import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.verifiers.Verifier;
import eu.solven.matmul.catalog.SchemeIO;

/**
 * One-shot migration that suffixes each scheme filename with its addition
 * count: {@code _r{rank}_a{adds}}. Two schemes at the same rank can have very
 * different addition counts — surfacing both in the filename makes
 * improvements visible at a glance (Laderman r=23 has 98 adds;
 * dronperminov-c88 r=23 has 88 adds — same multiplication count, fewer
 * additions = better).
 *
 * <p>Idempotent: files already named {@code …_r{N}_a{M}.json} are skipped.</p>
 *
 * <p>Run via:</p>
 * <pre>
 *   mvn -q test-compile
 *   java -cp target/classes:target/test-classes:&dollar;CLASSPATH \
 *        eu.solven.matmul.catalog.RenameSchemesAddAdds
 * </pre>
 */
public class RenameSchemesAddAdds {

	private static final File ROOT = new File("src/main/resources/schemes");
	/** Pattern for filenames already in the new {@code _a{adds}} form — skip these. */
	private static final Pattern NEW_FORM = Pattern.compile(".*_(?:r|m)\\d+_a\\d+.*\\.json");
	/** Pattern that finds the {@code _r{rank}} (or {@code _m{rank}}) tail to insert {@code _a{adds}} after. */
	private static final Pattern RANK_TAIL = Pattern.compile("(.*_(?:r|m)\\d+)(.*\\.json)");

	public static void main(String[] args) throws IOException {
		List<File> files;
		try (var s = Files.walk(ROOT.toPath())) {
			files = s.filter(p -> p.toString().endsWith(".json"))
					.map(Path::toFile)
					.sorted()
					.toList();
		}

		int renamed = 0, skipped = 0, errors = 0;
		List<String> errorList = new ArrayList<>();
		for (File f : files) {
			if (NEW_FORM.matcher(f.getName()).matches()) {
				skipped++;
				continue;
			}
			try {
				int adds = computeAdds(f);
				Matcher m = RANK_TAIL.matcher(f.getName());
				if (!m.matches()) {
					errorList.add(f.getName() + ": no _r tail to suffix");
					errors++;
					continue;
				}
				String newName = m.group(1) + "_a" + adds + m.group(2);
				File dst = new File(f.getParentFile(), newName);
				Files.move(f.toPath(), dst.toPath(), StandardCopyOption.REPLACE_EXISTING);
				renamed++;
			} catch (Exception e) {
				errorList.add(f.getName() + ": " + e.getMessage());
				errors++;
			}
		}

		System.out.printf("Rename: renamed=%d already-OK=%d errors=%d (total %d)%n",
				renamed, skipped, errors, files.size());
		if (!errorList.isEmpty()) {
			System.out.println("Errors:");
			errorList.forEach(s -> System.out.println("  " + s));
		}
	}

	private static int computeAdds(File f) throws IOException {
		tools.jackson.databind.JsonNode root = SchemeIO.parseJson(f);
		if (SchemeIO.isComplex(root)) {
			return Verifier.additionCount(SchemeIO.readComplex(root));
		}
		if (SchemeIO.isReduced(root)) {
			return Verifier.additionCount(SchemeIO.readReduced(root));
		}
		return Verifier.additionCount(SchemeIO.read(root));
	}
}
