package eu.solven.matmul.docs.explore;

import eu.solven.matmul.catalog.FieldAwareLookup;

/** Print {@code findRank} + the winning file for shapes, over a field: {@code --field=Q 31x31x31 30x31x31 …}. */
public final class ProbeFindRank {
	private ProbeFindRank() {}

	public static void main(String[] args) {
		String field = "Q";
		for (String a : args) if (a.startsWith("--field=")) field = a.substring(8);
		FieldAwareLookup lk = new FieldAwareLookup(field);
		for (String a : args) {
			if (a.startsWith("--")) continue;
			String[] d = a.split("x");
			int n = Integer.parseInt(d[0]), m = Integer.parseInt(d[1]), p = Integer.parseInt(d[2]);
			var ws = lk.findWithSource(n, m, p);
			System.out.printf("%s findRank(%s)=%d  file=%s%n", a, field, lk.findRank(n, m, p),
					ws.map(w -> w.path().getFileName().toString()).orElse("-"));
		}
	}
}
