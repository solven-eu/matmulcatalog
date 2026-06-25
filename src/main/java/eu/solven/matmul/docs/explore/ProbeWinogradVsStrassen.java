package eu.solven.matmul.docs.explore;
import java.io.File; import java.util.*;
import eu.solven.matmul.NonCubicBilinearAlgorithm; import eu.solven.matmul.SymmetryTransforms;
import eu.solven.matmul.algebra.*; import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.recombination.*;
import eu.solven.matmul.recombination.AnalyticalMaskSearch.SchemeSupports;
import eu.solven.matmul.recombination.Recombination.SotaResolver;
import eu.solven.matmul.recombination.RecombinationMultisetOrbit.Result;
import eu.solven.matmul.search.SearchBudget;
/** ISOLATE one base's own GL orbit: does any frontier member STRICTLY beat the native support
 *  (both with B&B, both over s3 orientations) on composed-mul count? */
public final class ProbeWinogradVsStrassen {
  public static void main(String[] a) throws Exception {
    String path = a.length>0? a[0] : "src/main/resources/schemes/known/section2/2x2x2-r7-strassen-db11bcc.json";
    NonCubicBilinearAlgorithm base = SchemeIO.read(new File(path));
    SotaResolver sota = Recombination.catalogResolver(Algebra.nonCommutative(Field.R));
    int md=Math.max(base.n,Math.max(base.m,base.p));
    Result orbit = md<=3? RecombinationMultisetOrbit.enumerate(base,2): RecombinationMultisetOrbit.enumerateSampled(base,80_000,2);
    // native (+ s3 orientations)
    List<SchemeSupports> nat=new ArrayList<>();
    for (NonCubicBilinearAlgorithm o: SymmetryTransforms.s3Orbit(base)) nat.add(SchemeSupports.extract(o));
    // every frontier member (+ s3 orientations)
    List<SchemeSupports> front=new ArrayList<>();
    for (String key: orbit.dominanceFrontier()){
      int[][][] xyz=orbit.representativeTransforms.get(key); if(xyz==null)continue;
      NonCubicBilinearAlgorithm mem=RecombinationMultisetOrbit.materialise(base,xyz[0],xyz[1],xyz[2]);
      for (NonCubicBilinearAlgorithm o: SymmetryTransforms.s3Orbit(mem)) front.add(SchemeSupports.extract(o));
    }
    System.out.printf("base %s ⟨%d,%d,%d⟩=%d  frontier=%d  native-supports=%d  frontier-supports=%d%n%n",
      new File(path).getName(), base.n,base.m,base.p,base.r, orbit.dominanceFrontier().size(), nat.size(), front.size());
    int wins=0, tested=0;
    for (int N=base.n; N<=base.n*8; N+=base.n) for (int M=base.m; M<=base.m*8; M+=base.m) for (int P=base.p; P<=base.p*8; P+=base.p){
      long bn=Long.MAX_VALUE; for (SchemeSupports s: nat) bn=Math.min(bn, AllocationOptimizer.optimize(s,sota,N,M,P,SearchBudget.EXACT,null).rank());
      long bf=Long.MAX_VALUE; for (SchemeSupports s: front) bf=Math.min(bf, AllocationOptimizer.optimize(s,sota,N,M,P,SearchBudget.EXACT,null).rank());
      tested++;
      if (bf<bn){ wins++; System.out.printf("  ⟨%d,%d,%d⟩  native=%d  orbit-member=%d  <<< MEMBER beats native by %d%n",N,M,P,bn,bf,bn-bf); }
    }
    System.out.printf("%nORBIT-MEMBER strictly beats NATIVE (same base, both B&B+orientations): %d / %d targets%n", wins, tested);
  }
}
