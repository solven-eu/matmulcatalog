package eu.solven.matmul.benchmark;

import lombok.extern.slf4j.Slf4j;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

/**
 * Run a (problem × strategy) matrix and emit two roll-up files alongside the
 * per-run artifacts:
 *
 * <ul>
 *   <li>{@code target/benchmarks/index.csv} — one row per run, full precision,
 *       machine-readable.</li>
 *   <li>{@code target/benchmarks/index.md} — pivot table (problem rows × strategy
 *       columns) with rounded timings, for at-a-glance human review.</li>
 * </ul>
 *
 * <p>CLI usage:</p>
 * <pre>
 *   java eu.solven.matmul.benchmark.BenchmarkSuite [tier]
 *     tier ∈ {trivial, fast, slow, research, all}    default: fast
 * </pre>
 */
@Slf4j
public final class BenchmarkSuite {

	public static void main(String[] args) throws IOException {
		String tier = args.length > 0 ? args[0] : "fast";
		Path outputRoot = Paths.get("target", "benchmarks");
		Files.createDirectories(outputRoot);

		List<BenchmarkProblem> problems = pickProblems(tier);
		List<BenchmarkStrategy> strategies = defaultStrategies();

		log.info(String.format("BenchmarkSuite: tier=%s, %d problems × %d strategies = %d runs%n",
				tier, problems.size(), strategies.size(),
				problems.size() * strategies.size()));
		log.info(String.format("output: %s%n%n", outputRoot.toAbsolutePath()));

		BenchmarkRunner runner = new BenchmarkRunner(outputRoot);
		List<RunRecord> records = new ArrayList<>();
		for (BenchmarkProblem problem : problems) {
			for (BenchmarkStrategy strategy : strategies) {
				log.info(String.format("running %s × %s ... ", problem.id, strategy.id));
				System.out.flush();
				BenchmarkResult res = runner.run(problem, strategy);
				log.info(String.format("%s (wall=%dms)%n", res.verdict, res.wallMs));
				records.add(new RunRecord(problem, strategy, res));
			}
		}

		writeCsv(outputRoot.resolve("index.csv").toFile(), records);
		writeMd(outputRoot.resolve("index.md").toFile(), problems, strategies, records);

		log.info(String.format("%nwrote %s%n", outputRoot.resolve("index.csv")));
		log.info(String.format("wrote %s%n", outputRoot.resolve("index.md")));
	}

	private static List<BenchmarkProblem> pickProblems(String tier) {
		switch (tier.toLowerCase()) {
		case "trivial":
			return BenchmarkCatalog.trivial();
		case "fast":
			List<BenchmarkProblem> fast = new ArrayList<>(BenchmarkCatalog.trivial());
			fast.addAll(BenchmarkCatalog.fast());
			return fast;
		case "slow":
			List<BenchmarkProblem> slow = new ArrayList<>(BenchmarkCatalog.trivial());
			slow.addAll(BenchmarkCatalog.fast());
			slow.addAll(BenchmarkCatalog.slow());
			return slow;
		case "research":
		case "all":
			List<BenchmarkProblem> all = new ArrayList<>(BenchmarkCatalog.trivial());
			all.addAll(BenchmarkCatalog.fast());
			all.addAll(BenchmarkCatalog.slow());
			all.addAll(BenchmarkCatalog.research());
			return all;
		default:
			throw new IllegalArgumentException("unknown tier: " + tier);
		}
	}

	public static List<BenchmarkStrategy> defaultStrategies() {
		// `nolex_sat4j` is intentionally NOT included: without symmetry breaking
		// SAT4J hangs on dense ⟨2,2,2⟩ r=7 (and any harder problem). It works
		// for trivial-tier upper-triangular cases only; if you want to measure
		// the impact of the lex-ordering rule specifically, run a custom suite.
		List<BenchmarkStrategy> ss = new ArrayList<>();
		ss.add(new BenchmarkStrategy("lex_sat4j", true, false, SolverChoice.SAT4J));
		ss.add(new BenchmarkStrategy("lex_kissat", true, false, SolverChoice.KISSAT));
		ss.add(new BenchmarkStrategy("lex_breakid_kissat", true, true, SolverChoice.KISSAT));
		return ss;
	}

	// ───────────────────────────────────────────────────────────────────────────
	// Roll-up writers
	// ───────────────────────────────────────────────────────────────────────────

	private static void writeCsv(File f, List<RunRecord> records) throws IOException {
		try (PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(f)))) {
			pw.println("problem,alphabet,format,variant,rank,strategy,columnLex,useBreakId,solver,"
					+ "verdict,wallMs,encodeMs,breakIdMs,solveMs,"
					+ "finalVars,finalClauses,preBreakIdVars,preBreakIdClauses,message");
			for (RunRecord r : records) {
				pw.printf("%s,%s,%s,%s,%d,%s,%b,%b,%s,%s,%d,%d,%d,%d,%d,%d,%d,%d,%s%n",
						r.problem.id, r.problem.alphabet, r.problem.formatLabel(),
						r.problem.variantLabel, r.problem.rank,
						r.strategy.id, r.strategy.columnLex, r.strategy.useBreakId,
						r.strategy.solver,
						r.result.verdict, r.result.wallMs, r.result.encodeMs,
						r.result.breakIdMs, r.result.solveMs,
						r.result.finalVars, r.result.finalClauses,
						r.result.preBreakIdVars, r.result.preBreakIdClauses,
						r.result.message == null ? "" : csvEscape(r.result.message));
			}
		}
	}

	private static String csvEscape(String s) {
		if (s.contains(",") || s.contains("\"") || s.contains("\n")) {
			return "\"" + s.replace("\"", "\"\"") + "\"";
		}
		return s;
	}

	private static void writeMd(File f, List<BenchmarkProblem> problems,
			List<BenchmarkStrategy> strategies, List<RunRecord> records) throws IOException {
		Map<String, RunRecord> byCell = new LinkedHashMap<>();
		for (RunRecord r : records) {
			byCell.put(r.problem.id + "@" + r.strategy.id, r);
		}

		try (PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(f)))) {
			pw.println("# Benchmark roll-up");
			pw.println();
			pw.println("Each cell: `verdict @ wall-time`. `vars/cls` shows CNF size handed to the solver.");
			pw.println();

			// Section 1: verdict × time pivot
			pw.println("## Verdicts and timings");
			pw.println();
			pw.print("| problem | format | r | variant |");
			for (BenchmarkStrategy s : strategies) pw.printf(" %s |", s.id);
			pw.println();
			pw.print("|---|---|---|---|");
			for (int i = 0; i < strategies.size(); i++) pw.print("---|");
			pw.println();
			for (BenchmarkProblem p : problems) {
				pw.printf("| `%s` | %s | %d | %s |", p.id, p.formatLabel(), p.rank, p.variantLabel);
				for (BenchmarkStrategy s : strategies) {
					RunRecord r = byCell.get(p.id + "@" + s.id);
					pw.printf(" %s |", formatCell(r));
				}
				pw.println();
			}
			pw.println();

			// Section 2: CNF size impact
			pw.println("## CNF sizes (impact of strategy on encoding)");
			pw.println();
			pw.println("| problem | strategy | pre-BreakID vars | pre-BreakID cls | final vars | final cls | growth |");
			pw.println("|---|---|---|---|---|---|---|");
			for (RunRecord r : records) {
				if (r.result.verdict == Verdict.NO_ENCODER || r.result.verdict == Verdict.UNAVAILABLE) continue;
				double growth = r.result.preBreakIdClauses == 0 ? 0.0
						: (double) r.result.finalClauses / r.result.preBreakIdClauses;
				pw.printf(Locale.ROOT, "| `%s` | %s | %d | %d | %d | %d | %.2f× |%n",
						r.problem.id, r.strategy.id,
						r.result.preBreakIdVars, r.result.preBreakIdClauses,
						r.result.finalVars, r.result.finalClauses, growth);
			}
		}
	}

	private static String formatCell(RunRecord r) {
		if (r == null) return "—";
		switch (r.result.verdict) {
		case SAT:
		case UNSAT:
			return String.format("**%s** %s", r.result.verdict, fmtTime(r.result.wallMs));
		case TIMEOUT:
			return String.format("_TIMEOUT_ >%s", fmtTime(r.result.wallMs));
		case NO_ENCODER:
		case UNAVAILABLE:
		case ERROR:
			return String.format("_%s_", r.result.verdict);
		default:
			return r.result.verdict.toString();
		}
	}

	private static String fmtTime(long ms) {
		if (ms < 1000) return ms + "ms";
		if (ms < 60_000) return String.format(Locale.ROOT, "%.1fs", ms / 1000.0);
		return String.format(Locale.ROOT, "%dm%02ds", ms / 60_000, (ms / 1000) % 60);
	}

	private static final class RunRecord {
		final BenchmarkProblem problem;
		final BenchmarkStrategy strategy;
		final BenchmarkResult result;

		RunRecord(BenchmarkProblem problem, BenchmarkStrategy strategy, BenchmarkResult result) {
			this.problem = problem;
			this.strategy = strategy;
			this.result = result;
		}
	}
}
