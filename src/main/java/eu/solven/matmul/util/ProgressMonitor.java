package eu.solven.matmul.util;

import java.lang.management.GarbageCollectorMXBean;
import java.lang.management.ManagementFactory;
import java.lang.management.MemoryMXBean;
import java.lang.management.MemoryUsage;
import java.util.List;
import java.util.Locale;
import java.util.function.LongSupplier;
import java.util.function.LongUnaryOperator;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Reusable progress reporter for long-running drivers.
 *
 * <p>Per the project convention (see {@code CLAUDE.md} §"Long-running
 * procedures"), any procedure that takes more than ~30s wall-clock should
 * emit periodic {@code [progress] X processed (… counters …) Yms elapsed}
 * lines with an ETA — and that ETA <strong>must adapt to workload
 * non-uniformity</strong>. A flat per-item throughput is wrong for cubic
 * shape sweeps where the last batches can be 30× heavier than the first.</p>
 *
 * <h2>Workload-aware ETA</h2>
 *
 * <p>The caller may supply a {@link LongUnaryOperator} mapping each
 * item-index to its predicted work weight (e.g. for cubic ⟨n,n,n⟩ pass
 * {@code i -> (long) shapes.get(i).n³}). The monitor then extrapolates
 * remaining time as</p>
 *
 * <pre>remaining_ms = elapsed_ms · (weight_left / weight_done)</pre>
 *
 * <p>If no weight function is supplied, weight defaults to 1 per item and
 * the ETA degrades gracefully to the naive linear formula.</p>
 *
 * <h2>Heap statistics</h2>
 *
 * <p>Each progress line includes used / committed / max heap (MiB) and the
 * delta of GC count + GC time since {@link #start()}. The numbers are
 * read from {@link MemoryMXBean} / {@link GarbageCollectorMXBean} and are
 * cheap (no JVM stop-the-world).</p>
 *
 * <h2>Usage</h2>
 *
 * <pre>{@code
 * ProgressMonitor m = ProgressMonitor.builder()
 *     .total(shapes.size())
 *     .weight(i -> {
 *         int[] s = shapes.get((int) i);
 *         return (long) s[0] * s[1] * s[2];
 *     })
 *     .logEveryMillis(2_000)
 *     .logger(log)
 *     .label("FMM compare")
 *     .build();
 * m.start();
 * for (int i = 0; i < shapes.size(); i++) {
 *     doWork(shapes.get(i));
 *     m.tick(i);
 * }
 * m.done();
 * }</pre>
 *
 * <p>Thread-safe: {@link #tick(long)} may be called concurrently. Only
 * one thread will emit a log line per interval (loser-yields).</p>
 */
public final class ProgressMonitor {

	/** Construct via {@link #builder()}. */
	public static Builder builder() {
		return new Builder();
	}

	private final long total;
	private final LongUnaryOperator weight;
	private final long totalWeight;
	private final long logEveryMillis;
	private final Logger logger;
	private final String label;
	/** Pluggable clock for tests; defaults to {@link System#nanoTime}. */
	private final LongSupplier nanoClock;
	private final MemoryMXBean memoryBean;
	private final List<GarbageCollectorMXBean> gcBeans;

	// Mutable state — guarded by `this` for the log-emission section.
	private long t0Nanos;
	private long lastLogNanos;
	private long baseGcCount;
	private long baseGcTimeMs;
	private long processedCount;
	private long processedWeight;

	private ProgressMonitor(Builder b) {
		this.total = b.total;
		this.weight = b.weight != null ? b.weight : i -> 1L;
		this.logEveryMillis = b.logEveryMillis;
		this.logger = b.logger != null ? b.logger : LoggerFactory.getLogger(ProgressMonitor.class);
		this.label = b.label != null ? b.label : "progress";
		this.nanoClock = b.nanoClock != null ? b.nanoClock : System::nanoTime;
		this.memoryBean = b.memoryBean != null ? b.memoryBean : ManagementFactory.getMemoryMXBean();
		this.gcBeans = b.gcBeans != null ? b.gcBeans : ManagementFactory.getGarbageCollectorMXBeans();

		// Precompute total weight if we have a finite total. For total<0
		// (unknown), the ETA is suppressed.
		long sum = 0L;
		if (b.total >= 0 && b.weight != null) {
			for (long i = 0; i < b.total; i++) {
				sum += b.weight.applyAsLong(i);
			}
		} else if (b.total >= 0) {
			sum = b.total;
		}
		this.totalWeight = sum;
	}

	/** Snapshot the clock + GC baseline. Must be called before {@link #tick}. */
	public synchronized void start() {
		this.t0Nanos = nanoClock.getAsLong();
		this.lastLogNanos = t0Nanos;
		this.baseGcCount = totalGcCount();
		this.baseGcTimeMs = totalGcTime();
		this.processedCount = 0;
		this.processedWeight = 0;
	}

	/**
	 * Record that item {@code itemIndex} has just completed. Emits a log
	 * line at most once per {@link Builder#logEveryMillis(long)}.
	 *
	 * @param itemIndex 0-based index of the item that completed.
	 */
	public void tick(long itemIndex) {
		long w = weight.applyAsLong(itemIndex);
		synchronized (this) {
			processedCount++;
			processedWeight += w;
			long now = nanoClock.getAsLong();
			boolean enough = (now - lastLogNanos) / 1_000_000L >= logEveryMillis;
			boolean isLast = total >= 0 && processedCount >= total;
			if (enough || isLast) {
				lastLogNanos = now;
				emit(now);
			}
		}
	}

	/** Force a final progress line at the end of the run. */
	public synchronized void done() {
		emit(nanoClock.getAsLong());
	}

	/**
	 * Format the line that would be emitted right now. Exposed for tests
	 * that want to assert on the exact string without binding to a logger.
	 */
	public synchronized String format() {
		return buildLine(nanoClock.getAsLong());
	}

	private void emit(long now) {
		logger.info("{}", buildLine(now));
	}

	private String buildLine(long now) {
		long elapsedMs = (now - t0Nanos) / 1_000_000L;
		StringBuilder sb = new StringBuilder(160);
		sb.append("[progress][").append(label).append("] ");
		if (total >= 0) {
			sb.append(processedCount).append('/').append(total);
		} else {
			sb.append(processedCount);
		}
		sb.append(" processed");
		// Weighted ETA.
		if (total >= 0 && totalWeight > 0 && processedWeight > 0) {
			double pctByWeight = 100.0 * processedWeight / totalWeight;
			long remainingMs = (long) (elapsedMs
					* ((double) (totalWeight - processedWeight) / processedWeight));
			sb.append(String.format(Locale.ROOT, " (%.1f%% by weight, ~%ds remaining)",
					pctByWeight, remainingMs / 1000L));
		}
		sb.append(' ').append(elapsedMs).append("ms elapsed");
		sb.append("  ").append(formatHeap());
		return sb.toString();
	}

	/** Format heap + GC delta. Visible for testing. */
	String formatHeap() {
		MemoryUsage h = memoryBean.getHeapMemoryUsage();
		long used = h.getUsed() >> 20;       // MiB
		long committed = h.getCommitted() >> 20;
		long max = h.getMax() >> 20;          // may be -1 on some JVMs
		long gcCount = totalGcCount() - baseGcCount;
		long gcTime = totalGcTime() - baseGcTimeMs;
		String maxStr = max < 0 ? "?" : (max + "M");
		return String.format(Locale.ROOT, "heap=%dM/%dM/%s gc=%d/+%dms",
				used, committed, maxStr, gcCount, gcTime);
	}

	private long totalGcCount() {
		long s = 0;
		for (GarbageCollectorMXBean b : gcBeans) {
			long c = b.getCollectionCount();
			if (c > 0) s += c;
		}
		return s;
	}

	private long totalGcTime() {
		long s = 0;
		for (GarbageCollectorMXBean b : gcBeans) {
			long c = b.getCollectionTime();
			if (c > 0) s += c;
		}
		return s;
	}

	/** Builder. Use {@link ProgressMonitor#builder()}. */
	public static final class Builder {
		private long total = -1;
		private LongUnaryOperator weight;
		private long logEveryMillis = 2_000L;
		private Logger logger;
		private String label;
		private LongSupplier nanoClock;
		private MemoryMXBean memoryBean;
		private List<GarbageCollectorMXBean> gcBeans;

		/** Total number of items, or negative for "unknown". */
		public Builder total(long total) { this.total = total; return this; }

		/**
		 * Per-item work weight (e.g. {@code i -> n*m*p} for shape ⟨n,m,p⟩).
		 * Used to compute a workload-aware ETA. Default: 1 per item.
		 */
		public Builder weight(LongUnaryOperator weight) { this.weight = weight; return this; }

		/** Throttle: emit at most one line per N ms. Default: 2000. */
		public Builder logEveryMillis(long ms) { this.logEveryMillis = ms; return this; }

		public Builder logger(Logger logger) { this.logger = logger; return this; }

		/** Human label appearing in {@code [progress][label]}. */
		public Builder label(String label) { this.label = label; return this; }

		/** Inject a synthetic clock (nanoseconds) for tests. */
		public Builder nanoClock(LongSupplier clock) { this.nanoClock = clock; return this; }

		/** Inject memory bean for tests. */
		public Builder memoryBean(MemoryMXBean bean) { this.memoryBean = bean; return this; }

		/** Inject GC beans for tests. */
		public Builder gcBeans(List<GarbageCollectorMXBean> beans) { this.gcBeans = beans; return this; }

		public ProgressMonitor build() { return new ProgressMonitor(this); }
	}
}
