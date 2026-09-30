package eu.solven.matmul.io;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.math.BigInteger;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.charset.StandardCharsets;
import java.util.HashMap;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import eu.solven.matmul.NonCubicBilinearAlgorithm;

/**
 * Reader for the {@code .npz} sparse-rational scheme archives published by
 * Khoruzhii–Serafin–Gelß–Pokutta (<a href="https://github.com/khoruzhii/lita">khoruzhii/lita</a>,
 * format id {@code fmm.q.csr.v1}). No Python/NumPy dependency: the archive is a
 * plain zip of {@code .npy} files, parsed here directly.
 *
 * <p>Layout (per their README): for each axis {@code u|v|w}, CSR arrays
 * {@code {axis}_indptr} (int64, R+1), {@code {axis}_indices} (int32, nnz),
 * {@code {axis}_numerators} / {@code {axis}_denominators} (int64, nnz); rows are
 * the R products, columns the ROW-MAJOR flattened matrix entries
 * ({@code A_i}, {@code B_j}, {@code C_k}); {@code T_ijk = Σ_q U_qi V_qj W_qk} with
 * {@code C_k = Σ_ij T_ijk A_i B_j}. That is exactly this catalog's in-memory
 * convention ({@code U[i·m+j][q]}, {@code V[j·p+l][q]}, {@code W[i·p+l][q]}, all
 * row-major — {@code SchemeIO.write} does the on-disk col-major-W conversion), so
 * the factors are transposed into {@code (flattened index) × (product)} and handed
 * to the public constructor unchanged. {@code metadata_json} (a 0-d {@code <U…}
 * UCS-4 string) carries {@code tensor:[n,m,p]}, {@code rank}, {@code coefficient_field}.</p>
 */
public final class NpzSchemeReader {
	private NpzSchemeReader() {}

	/** The parsed scheme plus the coefficient-denominator facts the field stamping needs. */
	public record Loaded(NonCubicBilinearAlgorithm alg, long maxDenominator, BigInteger lcmDenominator,
			String coefficientField, String format) {}

	public static Loaded read(File npz) throws IOException {
		Map<String, Npy> arrays = new HashMap<>();
		try (ZipFile zip = new ZipFile(npz)) {
			var entries = zip.entries();
			while (entries.hasMoreElements()) {
				ZipEntry e = entries.nextElement();
				if (!e.getName().endsWith(".npy")) continue;
				try (InputStream in = zip.getInputStream(e)) {
					arrays.put(e.getName().substring(0, e.getName().length() - 4), parseNpy(in.readAllBytes()));
				}
			}
		}
		Npy meta = req(arrays, "metadata_json");
		String json = meta.asString();
		int[] dims = intsAfter(json, "\"tensor\"");
		int rank = intsAfter(json, "\"rank\"")[0];
		String field = strAfter(json, "\"coefficient_field\"");
		String format = strAfter(json, "\"format\"");
		if (dims.length != 3) throw new IOException("metadata tensor must have 3 dims, got " + json);
		int n = dims[0], m = dims[1], p = dims[2];

		long[] maxDen = { 1 };
		BigInteger[] lcm = { BigInteger.ONE };
		double[][] U = axis(arrays, "u", rank, n * m, maxDen, lcm);
		double[][] V = axis(arrays, "v", rank, m * p, maxDen, lcm);
		double[][] W = axis(arrays, "w", rank, n * p, maxDen, lcm);
		return new Loaded(new NonCubicBilinearAlgorithm(n, m, p, U, V, W), maxDen[0], lcm[0], field, format);
	}

	/** CSR (products × flattened-entries) → dense (flattened-entries × products). */
	private static double[][] axis(Map<String, Npy> arrays, String name, int rank, int cols,
			long[] maxDen, BigInteger[] lcm) throws IOException {
		long[] indptr = req(arrays, name + "_indptr").asLongs();
		long[] indices = req(arrays, name + "_indices").asLongs();
		long[] num = req(arrays, name + "_numerators").asLongs();
		long[] den = req(arrays, name + "_denominators").asLongs();
		if (indptr.length != rank + 1) {
			throw new IOException(name + "_indptr length " + indptr.length + " != rank+1 " + (rank + 1));
		}
		double[][] out = new double[cols][rank];
		for (int q = 0; q < rank; q++) {
			for (long t = indptr[q]; t < indptr[q + 1]; t++) {
				int col = (int) indices[(int) t];
				if (col < 0 || col >= cols) throw new IOException(name + ": column " + col + " out of " + cols);
				long d = den[(int) t];
				if (d == 0) throw new IOException(name + ": zero denominator at " + t);
				if (d < 0) throw new IOException(name + ": negative denominator at " + t);
				maxDen[0] = Math.max(maxDen[0], d);
				BigInteger bd = BigInteger.valueOf(d);
				lcm[0] = lcm[0].divide(lcm[0].gcd(bd)).multiply(bd);
				out[col][q] += (double) num[(int) t] / (double) d;
			}
		}
		return out;
	}

	private static Npy req(Map<String, Npy> arrays, String key) throws IOException {
		Npy a = arrays.get(key);
		if (a == null) throw new IOException("npz lacks array '" + key + "' (has " + arrays.keySet() + ")");
		return a;
	}

	// ── .npy ────────────────────────────────────────────────────────────────

	/** A parsed {@code .npy}: dtype descriptor, shape, and the raw little-endian payload. */
	record Npy(String descr, int[] shape, byte[] data) {
		long[] asLongs() throws IOException {
			ByteBuffer bb = ByteBuffer.wrap(data).order(ByteOrder.LITTLE_ENDIAN);
			int count = 1;
			for (int s : shape) count *= s;
			long[] out = new long[count];
			switch (descr) {
				case "<i8" -> { for (int i = 0; i < count; i++) out[i] = bb.getLong(); }
				case "<i4" -> { for (int i = 0; i < count; i++) out[i] = bb.getInt(); }
				case "<i2" -> { for (int i = 0; i < count; i++) out[i] = bb.getShort(); }
				case "|i1", "<i1" -> { for (int i = 0; i < count; i++) out[i] = bb.get(); }
				case "<u4" -> { for (int i = 0; i < count; i++) out[i] = bb.getInt() & 0xFFFFFFFFL; }
				case "<u8" -> { for (int i = 0; i < count; i++) out[i] = bb.getLong(); }
				default -> throw new IOException("unsupported integer dtype " + descr);
			}
			return out;
		}

		/** A 0-d {@code <U{len}} (UCS-4 little-endian) string, NUL-padding stripped. */
		String asString() throws IOException {
			Matcher m = Pattern.compile("^<U(\\d+)$").matcher(descr);
			if (!m.matches()) throw new IOException("expected a <U dtype, got " + descr);
			int len = Integer.parseInt(m.group(1));
			ByteBuffer bb = ByteBuffer.wrap(data).order(ByteOrder.LITTLE_ENDIAN);
			StringBuilder sb = new StringBuilder(len);
			for (int i = 0; i < len && bb.remaining() >= 4; i++) {
				int cp = bb.getInt();
				if (cp == 0) break;
				sb.appendCodePoint(cp);
			}
			return sb.toString();
		}
	}

	private static final byte[] MAGIC = { (byte) 0x93, 'N', 'U', 'M', 'P', 'Y' };

	static Npy parseNpy(byte[] bytes) throws IOException {
		for (int i = 0; i < MAGIC.length; i++) {
			if (bytes[i] != MAGIC[i]) throw new IOException("not a .npy (bad magic)");
		}
		int major = bytes[6] & 0xFF;
		int headerLen, headerStart;
		if (major == 1) {
			headerLen = (bytes[8] & 0xFF) | ((bytes[9] & 0xFF) << 8);
			headerStart = 10;
		} else {
			headerLen = (bytes[8] & 0xFF) | ((bytes[9] & 0xFF) << 8) | ((bytes[10] & 0xFF) << 16) | ((bytes[11] & 0xFF) << 24);
			headerStart = 12;
		}
		String header = new String(bytes, headerStart, headerLen, StandardCharsets.UTF_8);
		Matcher d = Pattern.compile("'descr':\\s*'([^']+)'").matcher(header);
		if (!d.find()) throw new IOException("npy header without descr: " + header);
		if (header.contains("'fortran_order': True")) throw new IOException("fortran_order arrays unsupported");
		Matcher s = Pattern.compile("'shape':\\s*\\(([^)]*)\\)").matcher(header);
		if (!s.find()) throw new IOException("npy header without shape: " + header);
		String[] parts = s.group(1).trim().isEmpty() ? new String[0] : s.group(1).split(",");
		int[] shape = java.util.Arrays.stream(parts).map(String::trim).filter(x -> !x.isEmpty())
				.mapToInt(Integer::parseInt).toArray();
		int dataStart = headerStart + headerLen;
		ByteArrayOutputStream bo = new ByteArrayOutputStream(bytes.length - dataStart);
		bo.write(bytes, dataStart, bytes.length - dataStart);
		return new Npy(d.group(1), shape, bo.toByteArray());
	}

	// ── tiny JSON helpers (the metadata is a flat, well-formed object) ──────

	/** Integers of {@code "key": 123} or {@code "key": [1, 2, 3]}. */
	private static int[] intsAfter(String json, String key) throws IOException {
		Matcher m = Pattern.compile(Pattern.quote(key) + "\\s*:\\s*(\\[[^\\]]*\\]|-?\\d+)").matcher(json);
		if (!m.find()) throw new IOException("metadata lacks " + key + ": " + json);
		Matcher d = Pattern.compile("-?\\d+").matcher(m.group(1));
		java.util.List<Integer> out = new java.util.ArrayList<>();
		while (d.find()) out.add(Integer.parseInt(d.group()));
		return out.stream().mapToInt(Integer::intValue).toArray();
	}

	private static String strAfter(String json, String key) {
		Matcher m = Pattern.compile(Pattern.quote(key) + "\\s*:\\s*\"([^\"]*)\"").matcher(json);
		return m.find() ? m.group(1) : null;
	}
}
