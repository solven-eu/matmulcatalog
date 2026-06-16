package eu.solven.matmul;

import static org.assertj.core.api.Assertions.assertThat;

import org.junit.jupiter.api.Test;

import lombok.Builder;
import lombok.Getter;
import lombok.RequiredArgsConstructor;
import lombok.Value;
import lombok.extern.slf4j.Slf4j;

/**
 * Smoke test: verifies the Lombok annotation processor is wired and the
 * generated methods (getters, builder, slf4j {@code log}) are available
 * at compile + runtime.
 */
@Slf4j
public class TestLombokWiring {

	@Value
	@Builder
	static class Sample {
		int n;
		String label;
	}

	@RequiredArgsConstructor
	@Getter
	static class Counter {
		private final int start;
	}

	@Test
	public void lombok_generated_methods_work() {
		Sample s = Sample.builder().n(3).label("hi").build();
		assertThat(s.getN()).isEqualTo(3);
		assertThat(s.getLabel()).isEqualTo("hi");

		Counter c = new Counter(7);
		assertThat(c.getStart()).isEqualTo(7);

		log.info("Lombok @Slf4j logger is wired: {}", s);
	}
}
