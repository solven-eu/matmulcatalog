package eu.solven.matmul.docs.migrate;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;

import eu.solven.matmul.catalog.MatrixJsonFormatter;
import eu.solven.matmul.catalog.SchemeIO;

import tools.jackson.databind.JsonNode;
import tools.jackson.databind.node.ObjectNode;

/**
 * One-shot migration: strip the retired {@code "z2"} boolean from every scheme
 * JSON.
 *
 * <p>The {@code z2} flag was always written {@code false} (all 11.5k files
 * carried {@code z2:false} — even F₂-valid integer schemes, whose
 * {@code fields[]} already lists {@code F2}). It was both redundant with the
 * unified {@code fields[]} / {@code fields_not[]} model and routinely wrong, so
 * the writers stopped emitting it and {@link SchemeIO#isZ2} now decides F₂-only
 * membership purely from {@code fields[]}. This driver removes the stale key
 * from the on-disk catalog so the JSON stops asserting a field fact it never
 * tracked correctly.</p>
 *
 * <p>Files are re-emitted through the canonical {@link MatrixJsonFormatter}, so
 * the diff is exactly the dropped {@code z2} line — coefficient tokens
 * (integers, exact {@code "p/q"} rationals, {@code [re,im]} complex pairs) are
 * preserved verbatim because the formatter is representation-preserving and
 * idempotent on already-canonical files.</p>
 *
 * <p>Default DRY-RUN. {@code --execute} writes.</p>
 *
 * <pre>mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.migrate.DropZ2Key [-Dexec.args=--execute]</pre>
 */
public final class DropZ2Key {
	private DropZ2Key() {}

	public static void main(String[] args) throws Exception {
		boolean execute = List.of(args).contains("--execute");
		Path root = Path.of("src/main/resources/schemes");

		List<Path> files;
		try (var s = Files.walk(root)) {
			files = s.filter(p -> p.toString().endsWith(".json")).collect(Collectors.toList());
		}
		System.out.println("Scanning " + files.size() + " scheme files (mode="
				+ (execute ? "EXECUTE" : "DRY-RUN") + ")…");

		AtomicInteger had = new AtomicInteger(), absent = new AtomicInteger();
		AtomicInteger rewritten = new AtomicInteger(), errors = new AtomicInteger();
		long t0 = System.currentTimeMillis();
		int done = 0;
		for (Path f : files) {
			done++;
			try {
				JsonNode parsed = SchemeIO.parseJson(f.toFile());
				if (!(parsed instanceof ObjectNode obj)) { errors.incrementAndGet(); continue; }
				if (!obj.has("z2")) { absent.incrementAndGet(); continue; }
				had.incrementAndGet();
				if (execute) {
					obj.remove("z2");
					Files.writeString(f, MatrixJsonFormatter.format(obj));
					rewritten.incrementAndGet();
				}
			} catch (Exception e) {
				errors.incrementAndGet();
				System.out.println("[error] " + f.getFileName() + ": " + e);
			}
			if (done % 2000 == 0) {
				long ms = System.currentTimeMillis() - t0;
				System.out.printf("[progress] %d/%d processed (%d had z2), %dms elapsed%n",
						done, files.size(), had.get(), ms);
			}
		}

		System.out.println("\n=== " + (execute ? "STRIPPED" : "PLAN") + " ===");
		System.out.println("had z2 key (removed):  " + had.get());
		System.out.println("no z2 key (untouched): " + absent.get());
		System.out.println("rewritten:             " + rewritten.get());
		System.out.println("errors (skipped):      " + errors.get());
		if (!execute) System.out.println("\n(DRY-RUN — pass --execute to write)");
	}
}
