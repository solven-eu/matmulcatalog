package eu.solven.matmul.docs.explore;
import java.io.File; import java.util.*;
import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.SymmetryTransforms;
import eu.solven.matmul.algebra.*;
import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.recombination.*;
import eu.solven.matmul.recombination.AnalyticalMaskSearch.SchemeSupports;
import eu.solven.matmul.recombination.Recombination.SotaResolver;
import eu.solven.matmul.isotropy.RecombinationMultisetOrbit.Result;
import eu.solven.matmul.search.SearchBudget;
import eu.solven.matmul.isotropy.RecombinationMultisetOrbit;
public final class ProbeRecombBnB {
  public static void main(String[] a) throws Exception {
    SotaResolver sota=Recombination.catalogResolver(Algebra.nonCommutative(Field.R));
    String[] files={
      "src/main/resources/schemes/known/section2/2x2x2-r7-strassen-db11bcc.json",
      "src/main/resources/schemes/known/section3/2x2x3-r11-alphatensor_Z-682e003.json",
      "src/main/resources/schemes/derived/section4/2x2x4-r14-derived-a446d11.json",
      "src/main/resources/schemes/known/section3/2x3x3-r15-alphatensor_Z-497eea7.json",
      "src/main/resources/schemes/known/section4/3x3x4-r29-alphatensor_Z-a5da05f.json"};
    // SEARCH-recomb supports = each base's native, all orientations.
    // FRONTIER-recomb supports = each base's frontier members, materialised, all orientations.
    List<SchemeSupports> searchSup=new ArrayList<>(), frontSup=new ArrayList<>();
    for (String f: files){
      NonCubicBilinearAlgorithm b=SchemeIO.read(new File(f));
      int md=Math.max(b.n,Math.max(b.m,b.p));
      for (NonCubicBilinearAlgorithm ob: SymmetryTransforms.s3Orbit(b)) searchSup.add(SchemeSupports.extract(ob));
      Result orbit = md<=3? RecombinationMultisetOrbit.enumerate(b,2): RecombinationMultisetOrbit.enumerateSampled(b,40_000,2);
      for (String key: orbit.dominanceFrontier()){
        int[][][] xyz=orbit.representativeTransforms.get(key); if(xyz==null)continue;
        NonCubicBilinearAlgorithm orbited=RecombinationMultisetOrbit.materialise(b,xyz[0],xyz[1],xyz[2]);
        for (NonCubicBilinearAlgorithm ob: SymmetryTransforms.s3Orbit(orbited)) frontSup.add(SchemeSupports.extract(ob));
      }
    }
    System.out.printf("search-recomb supports=%d  frontier-recomb supports=%d%n%n", searchSup.size(), frontSup.size());
    int[][] shapes={{12,12,12},{16,16,16},{20,20,20},{18,19,20}};
    BlockSplitSearch.findBestStrategy(8,8,8, BlockSplitSearch.defaultPool(), sota, false); // warmup
    for (int[] t: shapes){
      long s0=System.nanoTime(); long sr=Long.MAX_VALUE;
      for (SchemeSupports s: searchSup) sr=Math.min(sr, AllocationOptimizer.optimize(s,sota,t[0],t[1],t[2],SearchBudget.EXACT,null).rank());
      long st=System.nanoTime()-s0;
      long f0=System.nanoTime(); long fr=Long.MAX_VALUE;
      for (SchemeSupports s: frontSup) fr=Math.min(fr, AllocationOptimizer.optimize(s,sota,t[0],t[1],t[2],SearchBudget.EXACT,null).rank());
      long ft=System.nanoTime()-f0;
      System.out.printf("⟨%d,%d,%d⟩  searchRecomb=%d (%.0fms)  frontierRecomb=%d (%.0fms)  Δrank=%+d  frontier %.1fx %s%n",
        t[0],t[1],t[2], sr, st/1e6, fr, ft/1e6, fr-sr, ft/(double)st, ft>st?"SLOWER":"faster");
    }
  }
}
