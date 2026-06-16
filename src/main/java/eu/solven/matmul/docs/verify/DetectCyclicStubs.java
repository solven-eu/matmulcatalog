package eu.solven.matmul.docs.verify;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Deque;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.catalog.SchemeIO;

import tools.jackson.databind.JsonNode;

/**
 * Detect (and optionally delete) derived-stub lineage CYCLES — the class of
 * corruption that makes {@code LineageReplayer} recurse forever (now a catchable
 * {@code IllegalStateException} thanks to the cycle guard, previously an
 * uncatchable {@code StackOverflowError}).
 *
 * <p>The graph is over FILES, with each lineage ref resolved through the real
 * orientation-aware {@link FieldAwareLookup#findFile} — this is load-bearing:
 * the canonical ⟨26,29,28⟩ corruption references shape ⟨26,28,29⟩ (a different
 * ordering) which {@code findFile} aliases back to the SAME file, so a pure
 * string graph over shapes would miss it.</p>
 *
 * <p>A stub is <b>bad</b> if it is on a cycle OR can transitively reach one — its
 * replay would throw. We report the cycle SCCs and the full bad set, and with
 * {@code --delete} remove them (the gap sweep then refills the shapes cleanly).</p>
 *
 * <pre>mvn -q -ntp exec:java -Dexec.mainClass=eu.solven.matmul.docs.verify.DetectCyclicStubs [-Dexec.args=--delete]</pre>
 */
public final class DetectCyclicStubs {
	private DetectCyclicStubs() {}

	private static final Pattern SHAPE = Pattern.compile("(\\d+)x(\\d+)x(\\d+)");

	public static void main(String[] args) throws Exception {
		boolean delete = List.of(args).contains("--delete");
		FieldAwareLookup lookup = new FieldAwareLookup("Q");
		Path root = Path.of("src/main/resources/schemes");

		List<Path> stubs;
		try (var s = Files.walk(root)) {
			stubs = s.filter(p -> p.getFileName().toString().startsWith("derived_recursive-")
					&& p.toString().endsWith(".json")).collect(Collectors.toList());
		}
		System.out.println("Scanning " + stubs.size() + " derived stubs for lineage cycles…");

		// file (canonical abs path) -> list of dependency files (stubs only)
		Map<String, List<String>> graph = new HashMap<>();
		Map<String, Path> pathOf = new HashMap<>();
		Map<String, Optional<String>> resolveCache = new HashMap<>();

		int processed = 0;
		for (Path stub : stubs) {
			String key = stub.toAbsolutePath().normalize().toString();
			pathOf.put(key, stub);
			Set<String> deps = new LinkedHashSet<>();
			for (int[] sh : referencedShapes(stub)) {
				if (sh[0] == 1 || sh[1] == 1 || sh[2] == 1) continue; // naive / elementary leaf
				String depKey = resolveCache.computeIfAbsent(sh[0] + "x" + sh[1] + "x" + sh[2], k -> {
					Optional<Path> f = lookup.findFile(sh[0], sh[1], sh[2]);
					return f.map(pp -> pp.toAbsolutePath().normalize().toString());
				}).orElse(null);
				if (depKey != null && depKey.contains("derived_recursive-")) deps.add(depKey);
			}
			graph.put(key, new ArrayList<>(deps));
			if (++processed % 2000 == 0) {
				System.out.println("[progress] " + processed + "/" + stubs.size() + " graphed");
			}
		}

		// Tarjan-free cycle detection: iterative DFS with grey/black colours; on a
		// back-edge to a grey node, walk the active stack to extract the cycle.
		Set<String> onCycle = new HashSet<>();
		List<List<String>> cycles = new ArrayList<>();
		Map<String, Integer> colour = new HashMap<>(); // 0 white,1 grey,2 black
		for (String start : graph.keySet()) {
			if (colour.getOrDefault(start, 0) != 0) continue;
			Deque<String> stack = new ArrayDeque<>();
			Deque<Integer> idx = new ArrayDeque<>();
			Deque<String> path = new ArrayDeque<>(); // active grey path (insertion order)
			stack.push(start); idx.push(0);
			while (!stack.isEmpty()) {
				String u = stack.peek();
				int i = idx.pop();
				if (i == 0) { colour.put(u, 1); path.addLast(u); }
				List<String> adj = graph.getOrDefault(u, List.of());
				if (i < adj.size()) {
					idx.push(i + 1);
					String v = adj.get(i);
					int c = colour.getOrDefault(v, 0);
					if (c == 0) { stack.push(v); idx.push(0); }
					else if (c == 1) { // back-edge → cycle from v..u
						List<String> cyc = new ArrayList<>();
						boolean collecting = false;
						for (String node : path) {
							if (node.equals(v)) collecting = true;
							if (collecting) { cyc.add(node); onCycle.add(node); }
						}
						if (!cyc.isEmpty()) cycles.add(cyc);
					}
				} else {
					colour.put(u, 2); path.pollLast(); stack.pop();
				}
			}
		}

		// Bad set = onCycle ∪ everything that can reach a cycle node (replay would throw).
		Set<String> bad = new HashSet<>(onCycle);
		// reverse edges
		Map<String, List<String>> rev = new HashMap<>();
		for (var e : graph.entrySet())
			for (String v : e.getValue()) rev.computeIfAbsent(v, k -> new ArrayList<>()).add(e.getKey());
		Deque<String> q = new ArrayDeque<>(onCycle);
		while (!q.isEmpty()) {
			String u = q.poll();
			for (String pre : rev.getOrDefault(u, List.of()))
				if (bad.add(pre)) q.add(pre);
		}

		System.out.println("\n=== RESULT ===");
		System.out.println("cycles found: " + cycles.size());
		int show = Math.min(cycles.size(), 15);
		for (int i = 0; i < show; i++) {
			System.out.println("  cycle " + (i + 1) + ": " + cycles.get(i).stream()
					.map(k -> new File(k).getName().replace("derived_recursive-", "").replace(".json", ""))
					.collect(Collectors.joining(" → ")));
		}
		if (cycles.size() > show) System.out.println("  … " + (cycles.size() - show) + " more");
		System.out.println("files ON a cycle:            " + onCycle.size());
		System.out.println("files reaching a cycle (bad): " + bad.size() + "  (of " + stubs.size() + " stubs)");

		if (delete) {
			int deleted = 0;
			for (String k : bad) {
				try { Files.deleteIfExists(pathOf.get(k)); deleted++; } catch (Exception e) {
					System.out.println("  delete failed: " + pathOf.get(k) + " — " + e);
				}
			}
			System.out.println("\nDELETED " + deleted + " corrupt stubs. Re-sweep by gap to refill.");
		} else {
			System.out.println("\n(dry run — pass --delete to remove the " + bad.size() + " bad stubs)");
		}
	}

	/** All ⟨n,m,p⟩ shapes a stub's lineage references (collected from every "ref" string). */
	private static List<int[]> referencedShapes(Path stub) {
		List<int[]> out = new ArrayList<>();
		try {
			JsonNode root = SchemeIO.parseJson(stub.toFile());
			JsonNode lin = root.get("lineage");
			if (lin != null) collectRefs(lin, out);
		} catch (Exception e) {
			// unreadable lineage — treat as no deps (a separate corruption class)
		}
		return out;
	}

	private static void collectRefs(JsonNode node, List<int[]> out) {
		if (node == null) return;
		if (node.isObject()) {
			JsonNode ref = node.get("ref");
			if (ref != null && ref.isString()) {
				Matcher m = SHAPE.matcher(ref.asString());
				if (m.find()) out.add(new int[] { Integer.parseInt(m.group(1)),
						Integer.parseInt(m.group(2)), Integer.parseInt(m.group(3)) });
			}
			for (Map.Entry<String, JsonNode> e : node.properties()) collectRefs(e.getValue(), out);
		} else if (node.isArray()) {
			for (JsonNode c : node) collectRefs(c, out);
		}
	}
}
