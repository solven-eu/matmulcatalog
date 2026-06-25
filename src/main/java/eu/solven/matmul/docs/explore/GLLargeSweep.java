package eu.solven.matmul.docs.explore;
import java.io.File; import java.util.*;
import eu.solven.matmul.NonCubicBilinearAlgorithm; import eu.solven.matmul.SymmetryTransforms;
import eu.solven.matmul.algebra.*; import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.recombination.*;
import eu.solven.matmul.recombination.AnalyticalMaskSearch.SchemeSupports;
import eu.solven.matmul.recombination.Recombination.SotaResolver;
import eu.solven.matmul.recombination.RecombinationMultisetOrbit.Result;
import eu.solven.matmul.search.SearchBudget;
/**
 * LARGE-shape payoff test for the GL-orbit frontier: across the 17..BAND band (sorted n≤m≤p), does
 * the FULL GL-frontier expansion of the base pool STRICTLY beat the live {@code findBestStrategy}
 * (which only expands by AXIS_FLIP + hand-curated cousins)? Each GL support is costed by B&B with an
 * {@code upTo(poolRank)} budget so it only surfaces strict beats — the systematic re-discovery of the
 * ⟨17,17,17⟩=2930 class. Bounds, not proven optima (single-level recombination over a bounded base set).
 */
public final class GLLargeSweep {
  record Sup(String base, SchemeSupports s) {}
  public static void main(String[] a) throws Exception {
    int BAND = a.length>0? Integer.parseInt(a[0]) : 26;
    eu.solven.matmul.catalog.FieldAwareLookup look = new eu.solven.matmul.catalog.FieldAwareLookup(Field.R);
    // sota = FULL on-disk catalog (NOT Recombination.catalogResolver, which uses a small hardcoded
    // KnownAlgorithmCatalog and falls back to CUBIC a*b*c for missing shapes — e.g. ⟨9,9,9⟩→729).
    // Memoised: findRank walks the catalog index per call and the B&B hammers it.
    java.util.Map<Long,Integer> memo = new java.util.concurrent.ConcurrentHashMap<>();
    SotaResolver sota = (x,y,z) -> memo.computeIfAbsent((long)x*1000000L+y*1000L+z, k-> {
      int r = look.findRank(x,y,z); return r==SotaResolver.UNKNOWN_RANK? x*y*z : r; });
    // Default: just the ⟨2,2,2⟩ base (the proven 2930 vehicle) — 36 supports, tiny 2-part
    // allocation trees → fast. Pass arg[1]="multi" for the wider (slower) base set.
    // "duo" = ⟨2,2,2⟩+⟨2,2,3⟩ (both have two size-2 axes → small allocation trees, fast).
    // "multi" adds the 3-axis bases ⟨2,3,3⟩/⟨3,3,3⟩ — MUCH slower (3-part composition explosion).
    boolean multi = a.length>1 && a[1].equals("multi");
    boolean duo = a.length>1 && a[1].equals("duo");
    String[] bases = multi? new String[]{
      "src/main/resources/schemes/known/section2/2x2x2-r7-strassen-db11bcc.json",
      "src/main/resources/schemes/known/section3/2x2x3-r11-alphatensor_Z-682e003.json",
      "src/main/resources/schemes/known/section3/2x3x3-r15-alphatensor_Z-497eea7.json",
      "src/main/resources/schemes/derived/section3/3x3x3-r23-derived-b173cf2.json"}
      : duo? new String[]{
      "src/main/resources/schemes/known/section2/2x2x2-r7-strassen-db11bcc.json",
      "src/main/resources/schemes/known/section3/2x2x3-r11-alphatensor_Z-682e003.json"}
      : new String[]{"src/main/resources/schemes/known/section2/2x2x2-r7-strassen-db11bcc.json"};
    List<Sup> sups = new ArrayList<>();
    for (String f : bases) {
      NonCubicBilinearAlgorithm b = SchemeIO.read(new File(f));
      int md = Math.max(b.n, Math.max(b.m, b.p));
      // EXACT enumerate for md≤3 (complete frontier — sampled both under- and over-counts:
      // ⟨2,2,3⟩ sampled=15 vs exact=21, ⟨2,3,3⟩ sampled=226 vs exact=170). Sampled only for md≥4.
      Result orbit = md<=3? RecombinationMultisetOrbit.enumerate(b,2)
                          : RecombinationMultisetOrbit.enumerateSampled(b, 60_000, 2);
      int before = sups.size();
      String lbl = "⟨"+b.n+","+b.m+","+b.p+"⟩="+b.r;
      for (String key : orbit.dominanceFrontier()) {
        int[][][] xyz = orbit.representativeTransforms.get(key); if (xyz==null) continue;
        NonCubicBilinearAlgorithm mem = RecombinationMultisetOrbit.materialise(b, xyz[0],xyz[1],xyz[2]);
        for (NonCubicBilinearAlgorithm o : SymmetryTransforms.s3Orbit(mem)) sups.add(new Sup(lbl, SchemeSupports.extract(o)));
      }
      System.out.printf("[base] %s frontier=%d -> %d supports (total %d)%n", lbl, orbit.dominanceFrontier().size(), sups.size()-before, sups.size());
    }
    // sorted targets 17..BAND
    List<int[]> T = new ArrayList<>();
    for (int n=17;n<=BAND;n++) for(int m=n;m<=BAND;m++) for(int p=m;p<=BAND;p++) T.add(new int[]{n,m,p});
    System.out.printf("%n[sweep] %d supports x %d targets (band 17..%d, sorted)%n%n", sups.size(), T.size(), BAND);
    long t0=System.nanoTime(); int wins=0, reDerived=0; long workDone=0, workTotal=0;
    for (int[] t: T) workTotal += (long)t[0]*t[1]*t[2];
    for (int ti=0; ti<T.size(); ti++) {
      int[] t = T.get(ti); int N=t[0],M=t[1],P=t[2];
      int catR = look.findRank(N,M,P);
      long pool0 = catR==SotaResolver.UNKNOWN_RANK? Long.MAX_VALUE : catR;
      long gl=Long.MAX_VALUE; String glBase="?";
      for (Sup su: sups) {
        if (N<su.s.n||M<su.s.m||P<su.s.p) continue;
        long r = AllocationOptimizer.optimize(su.s, sota, N,M,P, SearchBudget.upTo(pool0), null).rank();
        if (r<gl){ gl=r; glBase=su.base; }
      }
      if (gl<pool0){ wins++; System.out.printf("  <<< WIN ⟨%d,%d,%d⟩  catalog=%d  GL=%d (-%d)  via %s%n", N,M,P, pool0, gl, pool0-gl, glBase); }
      else if (gl==pool0 && pool0!=Long.MAX_VALUE) reDerived++; // GL single-base recomb re-derives catalog SOTA
      workDone += (long)N*M*P;
      if ((ti+1)%15==0 || ti==T.size()-1){
        double el=(System.nanoTime()-t0)/1e9; double frac=workDone/(double)workTotal;
        System.out.printf("[progress] %d/%d targets, %d wins, %d re-derived, %.0fs elapsed, ~%.0fs remaining%n",
          ti+1, T.size(), wins, reDerived, el, frac>0? el*(1-frac)/frac : 0);
      }
    }
    System.out.printf("%n[done] GL-frontier vs catalog over %d band targets: %d STRICT WINS, %d re-derived SOTA (gl==catalog)%n", T.size(), wins, reDerived);
  }
}
