package eu.solven.matmul.util;

import static org.assertj.core.api.Assertions.assertThat;

import java.lang.management.GarbageCollectorMXBean;
import java.lang.management.MemoryMXBean;
import java.lang.management.MemoryUsage;
import java.util.List;
import java.util.concurrent.atomic.AtomicLong;
import javax.management.ObjectName;

import org.junit.jupiter.api.Test;

/**
 * Sanity-checks for {@link ProgressMonitor}: weighted ETA math and the
 * heap-stat formatting. Uses a synthetic nano clock so timings are
 * deterministic, and a stub {@link MemoryMXBean} so the heap-format
 * assertion isn't flaky.
 */
public class TestProgressMonitor {

	/**
	 * If item 0 has weight 1 and item 1 has weight 9, then after item 0
	 * completes at t=1s, the ETA must reflect that 90% of work remains
	 * — i.e. ~9s remaining, NOT ~1s as a naive count-based ETA would say.
	 */
	@Test
	public void weighted_eta_reflects_remaining_work_not_remaining_count() {
		AtomicLong clockNs = new AtomicLong(0L);
		ProgressMonitor m = ProgressMonitor.builder()
				.total(2)
				.weight(i -> i == 0 ? 1L : 9L)
				.logEveryMillis(0)
				.label("test")
				.nanoClock(clockNs::get)
				.memoryBean(stubMemoryBean(0, 0, 0))
				.gcBeans(List.of())
				.build();
		m.start();
		// Advance to 1s and tick item 0 (weight 1).
		clockNs.set(1_000_000_000L);
		m.tick(0);

		String line = m.format();
		// 10% weight done, 90% remaining; elapsed=1s → remaining≈9s.
		assertThat(line).contains("1/2 processed")
				.contains("10.0% by weight")
				.contains("~9s remaining")
				.contains("1000ms elapsed");
	}

	/**
	 * With uniform weight the ETA must equal {@code elapsed · remaining/done}.
	 * Sanity-check the degenerate case where the workload-aware formula
	 * collapses to the naive throughput formula.
	 */
	@Test
	public void uniform_weight_eta_matches_naive_formula() {
		AtomicLong clockNs = new AtomicLong(0L);
		ProgressMonitor m = ProgressMonitor.builder()
				.total(10)
				.logEveryMillis(0)
				.nanoClock(clockNs::get)
				.memoryBean(stubMemoryBean(0, 0, 0))
				.gcBeans(List.of())
				.build();
		m.start();
		for (int i = 0; i < 2; i++) {
			clockNs.set((i + 1) * 500_000_000L); // 0.5s per item
			m.tick(i);
		}
		// 2/10 done, 1s elapsed → 8 left → ~4s remaining.
		String line = m.format();
		assertThat(line).contains("2/10 processed")
				.contains("20.0% by weight")
				.contains("~4s remaining");
	}

	/** Heap stat must be formatted as {@code heap=usedM/committedM/maxM gc=count/+timems}. */
	@Test
	public void heap_format_includes_used_committed_max_and_gc_delta() {
		ProgressMonitor m = ProgressMonitor.builder()
				.total(1)
				.memoryBean(stubMemoryBean(100L << 20, 200L << 20, 400L << 20))
				.gcBeans(List.of(stubGcBean(3, 17)))
				.nanoClock(() -> 0L)
				.build();
		m.start();
		assertThat(m.formatHeap()).isEqualTo("heap=100M/200M/400M gc=0/+0ms");
	}

	/** When max is -1 (unbounded), the format should print '?'. */
	@Test
	public void heap_format_handles_unbounded_max() {
		ProgressMonitor m = ProgressMonitor.builder()
				.total(1)
				.memoryBean(stubMemoryBean(50L << 20, 80L << 20, -1L))
				.gcBeans(List.of())
				.nanoClock(() -> 0L)
				.build();
		m.start();
		assertThat(m.formatHeap()).isEqualTo("heap=50M/80M/? gc=0/+0ms");
	}

	/** GC counts after {@code start} should be reported as a delta, not the absolute. */
	@Test
	public void gc_counts_are_delta_from_start() {
		AtomicLong gcCount = new AtomicLong(5);
		AtomicLong gcTime = new AtomicLong(123);
		GarbageCollectorMXBean bean = stubGcBean(gcCount, gcTime);
		ProgressMonitor m = ProgressMonitor.builder()
				.total(1)
				.memoryBean(stubMemoryBean(0, 0, 0))
				.gcBeans(List.of(bean))
				.nanoClock(() -> 0L)
				.build();
		m.start();
		gcCount.set(8);
		gcTime.set(200);
		assertThat(m.formatHeap()).contains("gc=3/+77ms");
	}

	/** With total<0 (unknown), no ETA / pct is reported. */
	@Test
	public void unknown_total_skips_eta() {
		AtomicLong clockNs = new AtomicLong(0L);
		ProgressMonitor m = ProgressMonitor.builder()
				.total(-1)
				.logEveryMillis(0)
				.nanoClock(clockNs::get)
				.memoryBean(stubMemoryBean(0, 0, 0))
				.gcBeans(List.of())
				.build();
		m.start();
		clockNs.set(500_000_000L);
		m.tick(0);
		String line = m.format();
		assertThat(line).contains("1 processed")
				.contains("500ms elapsed")
				.doesNotContain("remaining")
				.doesNotContain("by weight");
	}

	// ---- Stubs ----

	private static MemoryMXBean stubMemoryBean(long used, long committed, long max) {
		MemoryUsage usage = new MemoryUsage(0L, used, committed, max < 0 ? -1L : max);
		return new MemoryMXBean() {
			@Override public MemoryUsage getHeapMemoryUsage() { return usage; }
			@Override public MemoryUsage getNonHeapMemoryUsage() { return usage; }
			@Override public int getObjectPendingFinalizationCount() { return 0; }
			@Override public void gc() {}
			@Override public boolean isVerbose() { return false; }
			@Override public void setVerbose(boolean v) {}
			@Override public ObjectName getObjectName() { return null; }
		};
	}

	private static GarbageCollectorMXBean stubGcBean(long count, long timeMs) {
		return stubGcBean(new AtomicLong(count), new AtomicLong(timeMs));
	}

	private static GarbageCollectorMXBean stubGcBean(AtomicLong count, AtomicLong timeMs) {
		return new GarbageCollectorMXBean() {
			@Override public long getCollectionCount() { return count.get(); }
			@Override public long getCollectionTime() { return timeMs.get(); }
			@Override public String getName() { return "stub"; }
			@Override public boolean isValid() { return true; }
			@Override public String[] getMemoryPoolNames() { return new String[0]; }
			@Override public ObjectName getObjectName() { return null; }
		};
	}
}
