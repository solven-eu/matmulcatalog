package eu.solven.matmul.docs.migrate;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import eu.solven.matmul.CanonicalShape;
import eu.solven.matmul.Shape;
import eu.solven.matmul.docs.migrate.ImportUpstreamSchemeRepos.Candidate;
import eu.solven.matmul.docs.migrate.UpstreamSchemeRepos.UpstreamSchemeRepo;

import tools.jackson.databind.JsonNode;
import tools.jackson.databind.json.JsonMapper;
import tools.jackson.databind.node.ObjectNode;

/**
 * The <b>rank digest</b> of the repositories in {@link UpstreamSchemeRepos}: for every
 * repository, the lowest rank it publishes per format (sorted shape) and the file that
 * carries it — read off the repository's file listing by the registry pattern, at
 * <em>every</em> dimension.
 *
 * <p><b>Why it exists.</b> A repository with an {@link UpstreamSchemeRepo#explicitMaxDim}
 * is only imported up to that dimension; what it publishes above is held by no catalog
 * file. Without this digest those schemes would be invisible here: no gap, no alert, a
 * catalog that reads as up to date while an explicit, verified scheme sits upstream below
 * our rank. The digest makes them a source of the upstream scan
 * ({@code ScanUpstreamSources}, kind {@code base-candidate}) — a claim read from a
 * filename, labelled as such, never a rank of ours.</p>
 *
 * <p>Written by {@link ImportUpstreamSchemeRepos} on every non-dry run (the 3-hourly sync);
 * repositories not listed in a run keep their previous entry. No timestamp: an unchanged
 * upstream rewrites the file byte for byte.</p>
 */
public final class UpstreamRepoRanks {
	private UpstreamRepoRanks() {}

	public static final Path FILE = Path.of("references/catalogs/upstream-repo-ranks.json");

	private static final JsonMapper MAPPER = JsonMapper.builder().build();

	/** One repository's best rank per format. */
	public record RepoRanks(String key, String repo, String source, int explicitMaxDim,
			Map<CanonicalShape, Integer> ranks, Map<CanonicalShape, String> paths) {
	}

	/** Lowest rank per sorted shape among the paths the repository's pattern names. */
	public static RepoRanks of(UpstreamSchemeRepo repo, List<String> upstreamPaths) {
		Map<CanonicalShape, Integer> ranks = new TreeMap<>(ORDER);
		Map<CanonicalShape, String> paths = new TreeMap<>(ORDER);
		for (String path : upstreamPaths) {
			Candidate c = repo.parse(path).orElse(null);
			if (c == null) continue;
			CanonicalShape sh = Shape.of(c.n(), c.m(), c.p()).canonical();
			Integer held = ranks.get(sh);
			if (held == null || c.rank() < held) {
				ranks.put(sh, c.rank());
				paths.put(sh, path);
			}
		}
		return new RepoRanks(repo.key(), repo.repo(), repo.source(), repo.explicitMaxDim(), ranks, paths);
	}

	private static final java.util.Comparator<CanonicalShape> ORDER = java.util.Comparator
			.comparingInt(CanonicalShape::n).thenComparingInt(CanonicalShape::m).thenComparingInt(CanonicalShape::p);

	/** Replace the entries of the listed repositories in {@code file}, keeping the others. */
	public static void merge(Path file, Map<UpstreamSchemeRepo, List<String>> listings) throws IOException {
		Map<String, RepoRanks> all = new TreeMap<>(read(file));
		listings.forEach((repo, paths) -> all.put(repo.key(), of(repo, paths)));
		ObjectNode root = MAPPER.createObjectNode();
		root.put("generated_by", ImportUpstreamSchemeRepos.class.getName());
		root.put("note", "Lowest rank per format (sorted shape) that each UpstreamSchemeRepos repository publishes, read"
				+ " from its file listing by the registry pattern. A rank here is the upstream's claim as named by its"
				+ " file - verified by this catalog only where the scheme was imported (max dimension <="
				+ " explicit_max_dim). Read by ScanUpstreamSources.");
		ObjectNode repos = root.putObject("repos");
		for (RepoRanks r : all.values()) {
			ObjectNode o = repos.putObject(r.key());
			o.put("repo", r.repo());
			o.put("source", r.source());
			if (r.explicitMaxDim() != Integer.MAX_VALUE) o.put("explicit_max_dim", r.explicitMaxDim());
			ObjectNode ranks = o.putObject("ranks");
			r.ranks().forEach((sh, rank) -> {
				ObjectNode e = ranks.putObject(sh.toFilenameToken());
				e.put("rank", rank);
				e.put("path", r.paths().get(sh));
			});
		}
		Files.createDirectories(file.toAbsolutePath().getParent());
		Files.writeString(file, MAPPER.writerWithDefaultPrettyPrinter().writeValueAsString(root) + "\n",
				StandardCharsets.UTF_8);
	}

	/** The digest by repository key; an absent file reads as empty. */
	public static Map<String, RepoRanks> read(Path file) throws IOException {
		Map<String, RepoRanks> out = new LinkedHashMap<>();
		if (!Files.isRegularFile(file)) return out;
		JsonNode repos = MAPPER.readTree(Files.readString(file, StandardCharsets.UTF_8)).path("repos");
		for (Map.Entry<String, JsonNode> e : repos.properties()) {
			Map<CanonicalShape, Integer> ranks = new TreeMap<>(ORDER);
			Map<CanonicalShape, String> paths = new TreeMap<>(ORDER);
			for (Map.Entry<String, JsonNode> s : e.getValue().path("ranks").properties()) {
				String[] d = s.getKey().split("x");
				CanonicalShape sh = Shape.of(Integer.parseInt(d[0]), Integer.parseInt(d[1]), Integer.parseInt(d[2]))
						.canonical();
				ranks.put(sh, s.getValue().path("rank").asInt());
				paths.put(sh, s.getValue().path("path").asString(""));
			}
			out.put(e.getKey(), new RepoRanks(e.getKey(), e.getValue().path("repo").asString(""),
					e.getValue().path("source").asString(""),
					e.getValue().path("explicit_max_dim").asInt(Integer.MAX_VALUE), ranks, paths));
		}
		return out;
	}
}
