package eu.solven.matmul.docs.verify;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;

import tools.jackson.databind.JsonNode;
import tools.jackson.databind.ObjectMapper;

/**
 * The manual FMM artifact audit ({@code references/fmm-artifact-audit.json}): formats
 * whose FMM <em>index</em> rank is not backed by FMM's own published tensor artifact
 * ({@code index_only}: the artifact's rank is ≥ ours; {@code placeholder}: no tensor
 * data; {@code artifact_invalid}: the artifact does not compute matmul). Such a row is
 * an upstream-unverified claim, not an actionable gap — see
 * {@code references/fmm-artifact-audit.md} for the method.
 *
 * <p>An entry is a statement about ONE index value and <b>expires</b> when the index
 * moves: it carries the {@code index_rank} it was audited against and applies only
 * while the digest still shows that rank (a legacy entry without {@code index_rank}
 * is always honoured). Shared by {@link FmmCrossCheck} (the committed report) and
 * {@code ScanUpstreamSources} (the weekly alert), so neither counts the same audited
 * row as a gap.</p>
 */
public final class FmmArtifactAudit {

	public static final File DEFAULT = new File("references/fmm-artifact-audit.json");

	/** One audited format: its class and the index value it was audited at ({@code null} = legacy). */
	public record Entry(String cls, Integer indexRank) {}

	private final Map<String, Entry> entries;

	private FmmArtifactAudit(Map<String, Entry> entries) {
		this.entries = entries;
	}

	/** Empty audit (no file → every WORSE row is a gap). */
	public static FmmArtifactAudit none() {
		return new FmmArtifactAudit(Map.of());
	}

	public static FmmArtifactAudit load() throws IOException {
		return load(DEFAULT);
	}

	public static FmmArtifactAudit load(File f) throws IOException {
		if (!f.isFile()) return none();
		JsonNode root = new ObjectMapper().readTree(f);
		Map<String, Entry> out = new LinkedHashMap<>();
		root.get("shapes").properties().forEach(e -> {
			JsonNode v = e.getValue();
			out.put(e.getKey(), new Entry(v.get("class").asString(),
					v.has("index_rank") ? v.get("index_rank").asInt() : null));
		});
		return new FmmArtifactAudit(out);
	}

	public int size() {
		return entries.size();
	}

	/**
	 * The audit class covering {@code ⟨n,m,p⟩} at FMM index rank {@code indexRank}, if a
	 * still-valid entry exists. Empty when the format is not audited OR its entry
	 * expired (index moved) — the caller treats the row as an ordinary gap.
	 */
	public Optional<String> coveringClass(int n, int m, int p, int indexRank) {
		Entry e = entries.get(key(n, m, p));
		if (e == null) return Optional.empty();
		if (e.indexRank() != null && e.indexRank() != indexRank) return Optional.empty();
		return Optional.of(e.cls());
	}

	/** The audited index value of an entry that no longer matches, for the expiry log. */
	public Optional<Integer> expiredAt(int n, int m, int p, int indexRank) {
		Entry e = entries.get(key(n, m, p));
		if (e == null || e.indexRank() == null || e.indexRank() == indexRank) return Optional.empty();
		return Optional.of(e.indexRank());
	}

	/** Every audited format key ({@code "nxmxp"}, sorted dims) — for the curation guard. */
	public List<String> keys() {
		return new ArrayList<>(entries.keySet());
	}

	static String key(int n, int m, int p) {
		int[] d = { n, m, p };
		Arrays.sort(d);
		return d[0] + "x" + d[1] + "x" + d[2];
	}
}
