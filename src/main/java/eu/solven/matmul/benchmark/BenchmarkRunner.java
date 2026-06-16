package eu.solven.matmul.benchmark;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.util.List;

import eu.solven.matmul.f2.sat.BreakIdBridge;
import eu.solven.matmul.f2.sat.Cnf;
import eu.solven.matmul.f2.sat.SatMatmulPipeline;
import eu.solven.matmul.f2.sat.SolverOutcome;
import eu.solven.matmul.f2.sat.Z2CnfEncoder;

/**
 * Runs a single ({@link BenchmarkProblem}, {@link BenchmarkStrategy}) pair and
 * persists artifacts under {@code <outputRoot>/<problem-id>/<strategy-id>/}:
 *
 * <ul>
 *   <li>{@code pre_solver.cnf} — DIMACS from the encoder, before BreakID.</li>
 *   <li>{@code solver_in.cnf} — DIMACS actually consumed by the solver (= pre_solver
 *       when BreakID is off; the BreakID output otherwise). Always written for
 *       reproducibility.</li>
 *   <li>{@code result.json} — verdict + timing + CNF size metrics.</li>
 *   <li>{@code decoded.txt} — the {U, V, W} factors when SAT (human-readable).</li>
 * </ul>
 *
 * <p>One run = one directory. Re-running overwrites; we don't timestamp by
 * default because the most-recent result is usually what matters; uncomment
 * the timestamp suffix in {@link #runDir} to keep a history.</p>
 */
public final class BenchmarkRunner {

	private final Path outputRoot;

	public BenchmarkRunner(Path outputRoot) {
		this.outputRoot = outputRoot;
	}

	public BenchmarkResult run(BenchmarkProblem problem, BenchmarkStrategy strategy) {
		Path dir = runDir(problem, strategy);
		try {
			Files.createDirectories(dir);
		} catch (IOException e) {
			return new BenchmarkResult(Verdict.ERROR, 0, 0, 0, 0, 0, 0, 0, 0,
					"failed to create output dir " + dir + ": " + e.getMessage());
		}

		// Alphabet gate. Z3/Z encoders not implemented yet.
		if (problem.alphabet != Alphabet.Z2) {
			BenchmarkResult res = BenchmarkResult.noEncoder("no encoder for " + problem.alphabet);
			writeResultJson(dir, problem, strategy, res);
			return res;
		}

		// Solver availability gate.
		if (!strategy.solver.isAvailable()) {
			BenchmarkResult res = BenchmarkResult.unavailable(strategy.solver + " binary not on PATH");
			writeResultJson(dir, problem, strategy, res);
			return res;
		}
		if (strategy.useBreakId && !BreakIdBridge.isAvailable()) {
			BenchmarkResult res = BenchmarkResult.unavailable("BreakID binary not on PATH");
			writeResultJson(dir, problem, strategy, res);
			return res;
		}

		// 1. Encode.
		long encodeStart = System.currentTimeMillis();
		Z2CnfEncoder encoder = new Z2CnfEncoder(
				problem.dimU(), problem.dimV(), problem.dimW(),
				problem.rank, problem.target, strategy.columnLex, problem.forceZeroSlots);
		int preBreakIdVars = encoder.getVarCount();
		List<int[]> preBreakIdClauses = encoder.getClauses();
		long encodeMs = System.currentTimeMillis() - encodeStart;

		File preSolverFile = dir.resolve("pre_solver.cnf").toFile();
		try {
			Cnf.writeDimacs(preSolverFile, preBreakIdVars, preBreakIdClauses);
		} catch (IOException e) {
			return error(dir, problem, strategy, encodeMs, 0, 0, preBreakIdVars,
					preBreakIdClauses.size(), preBreakIdVars, preBreakIdClauses.size(),
					"failed writing pre_solver.cnf: " + e.getMessage());
		}

		// 2. Optional BreakID pass.
		int finalVars = preBreakIdVars;
		List<int[]> finalClauses = preBreakIdClauses;
		long breakIdMs = 0;
		File solverInFile = dir.resolve("solver_in.cnf").toFile();
		if (strategy.useBreakId) {
			long breakIdStart = System.currentTimeMillis();
			try {
				Cnf.ReadResult augmented = BreakIdBridge.preprocess(preBreakIdVars, preBreakIdClauses);
				finalVars = augmented.varCount;
				finalClauses = augmented.clauses;
				breakIdMs = System.currentTimeMillis() - breakIdStart;
				Cnf.writeDimacs(solverInFile, finalVars, finalClauses);
			} catch (IOException e) {
				return error(dir, problem, strategy, encodeMs,
						System.currentTimeMillis() - breakIdStart, 0,
						preBreakIdVars, preBreakIdClauses.size(),
						preBreakIdVars, preBreakIdClauses.size(),
						"BreakID failed: " + e.getMessage());
			}
		} else {
			try {
				Files.copy(preSolverFile.toPath(), solverInFile.toPath(),
						StandardCopyOption.REPLACE_EXISTING);
			} catch (IOException e) {
				return error(dir, problem, strategy, encodeMs, 0, 0, finalVars,
						finalClauses.size(), preBreakIdVars, preBreakIdClauses.size(),
						"failed writing solver_in.cnf: " + e.getMessage());
			}
		}

		// 3. Solve.
		long solveStart = System.currentTimeMillis();
		SolverOutcome outcome = strategy.solver.solveBounded(finalVars, finalClauses, strategy.timeoutSec);
		long solveMs = System.currentTimeMillis() - solveStart;
		long wallMs = encodeMs + breakIdMs + solveMs;

		Verdict verdict;
		String message = null;
		switch (outcome.state) {
		case UNSAT:
			verdict = Verdict.UNSAT;
			break;
		case TIMEOUT:
			verdict = Verdict.TIMEOUT;
			message = "solver exceeded wall-time budget of " + strategy.timeoutSec + "s";
			break;
		case ERROR:
			verdict = Verdict.ERROR;
			message = outcome.message;
			break;
		case SAT:
			double[][][] uvw = encoder.decodeRaw(outcome.assignment.get());
			boolean ok = SatMatmulPipeline.verifyZ2NonCubic(uvw, problem.target);
			if (!ok) {
				BenchmarkResult res = new BenchmarkResult(Verdict.ERROR, wallMs, encodeMs,
						breakIdMs, solveMs, finalVars, finalClauses.size(),
						preBreakIdVars, preBreakIdClauses.size(),
						"SAT result failed verification");
				writeResultJson(dir, problem, strategy, res);
				return res;
			}
			verdict = Verdict.SAT;
			writeDecoded(dir, uvw);
			break;
		default:
			throw new IllegalStateException("unknown solver state " + outcome.state);
		}

		BenchmarkResult res = new BenchmarkResult(verdict, wallMs, encodeMs, breakIdMs,
				solveMs, finalVars, finalClauses.size(), preBreakIdVars,
				preBreakIdClauses.size(), message);
		writeResultJson(dir, problem, strategy, res);
		return res;
	}

	private Path runDir(BenchmarkProblem problem, BenchmarkStrategy strategy) {
		return outputRoot.resolve(problem.id).resolve(strategy.id);
	}

	private BenchmarkResult error(Path dir, BenchmarkProblem problem, BenchmarkStrategy strategy,
			long encodeMs, long breakIdMs, long solveMs,
			int finalVars, int finalClauses, int preBreakIdVars, int preBreakIdClauses,
			String message) {
		BenchmarkResult res = new BenchmarkResult(Verdict.ERROR,
				encodeMs + breakIdMs + solveMs, encodeMs, breakIdMs, solveMs,
				finalVars, finalClauses, preBreakIdVars, preBreakIdClauses, message);
		writeResultJson(dir, problem, strategy, res);
		return res;
	}

	private static void writeResultJson(Path dir, BenchmarkProblem problem,
			BenchmarkStrategy strategy, BenchmarkResult res) {
		File f = dir.resolve("result.json").toFile();
		try (PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(f)))) {
			pw.println("{");
			pw.printf("  \"problem\": \"%s\",%n", problem.id);
			pw.printf("  \"strategy\": \"%s\",%n", strategy.id);
			pw.printf("  \"alphabet\": \"%s\",%n", problem.alphabet);
			pw.printf("  \"format\": \"%s\",%n", problem.formatLabel());
			pw.printf("  \"variant\": \"%s\",%n", problem.variantLabel);
			pw.printf("  \"rank\": %d,%n", problem.rank);
			pw.printf("  \"columnLex\": %b,%n", strategy.columnLex);
			pw.printf("  \"useBreakId\": %b,%n", strategy.useBreakId);
			pw.printf("  \"solver\": \"%s\",%n", strategy.solver);
			pw.printf("  \"verdict\": \"%s\",%n", res.verdict);
			pw.printf("  \"wallMs\": %d,%n", res.wallMs);
			pw.printf("  \"encodeMs\": %d,%n", res.encodeMs);
			pw.printf("  \"breakIdMs\": %d,%n", res.breakIdMs);
			pw.printf("  \"solveMs\": %d,%n", res.solveMs);
			pw.printf("  \"finalVars\": %d,%n", res.finalVars);
			pw.printf("  \"finalClauses\": %d,%n", res.finalClauses);
			pw.printf("  \"preBreakIdVars\": %d,%n", res.preBreakIdVars);
			pw.printf("  \"preBreakIdClauses\": %d,%n", res.preBreakIdClauses);
			pw.printf("  \"message\": %s%n",
					res.message == null ? "null" : "\"" + escapeJson(res.message) + "\"");
			pw.println("}");
		} catch (IOException ignored) {
			// best-effort
		}
	}

	private static String escapeJson(String s) {
		return s.replace("\\", "\\\\").replace("\"", "\\\"").replace("\n", "\\n");
	}

	private static void writeDecoded(Path dir, double[][][] uvw) {
		File f = dir.resolve("decoded.txt").toFile();
		try (PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(f)))) {
			pw.println("# U (dimU × r)");
			writeFactor(pw, uvw[0]);
			pw.println("# V (dimV × r)");
			writeFactor(pw, uvw[1]);
			pw.println("# W (dimW × r)");
			writeFactor(pw, uvw[2]);
		} catch (IOException ignored) {
			// best-effort
		}
	}

	private static void writeFactor(PrintWriter pw, double[][] m) {
		for (int i = 0; i < m.length; i++) {
			StringBuilder sb = new StringBuilder();
			for (int k = 0; k < m[i].length; k++) {
				if (k > 0) sb.append(' ');
				sb.append((int) Math.round(m[i][k]));
			}
			pw.println(sb);
		}
	}
}
