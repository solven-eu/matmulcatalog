package eu.solven.matmul.docs.explore;
import java.io.File; import java.util.*;
import eu.solven.matmul.NonCubicBilinearAlgorithm; import eu.solven.matmul.SymmetryTransforms;
import eu.solven.matmul.algebra.*; import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.recombination.*;
import eu.solven.matmul.recombination.AnalyticalMaskSearch.SchemeSupports;
import eu.solven.matmul.recombination.Recombination.SotaResolver;
import eu.solven.matmul.recombination.RecombinationMultisetOrbit.Result;
import eu.solven.matmul.search.SearchBudget;
/** For ONE base: does the FULL GL-orbit frontier beat the live pool's AXIS_FLIP expansion?
 *  Both costed by B&B. AXIS_FLIP = 8 row-reversal variants (what defaultPool uses).
 *  GL frontier = full change-of-basis dominance. Sweep INCLUDES odd/unbalanced targets. */
public final class ProbeGlVsAxisFlip {
  static long bnb(SchemeSupports s, SotaResolver sota,int N,int M,int P){ return AllocationOptimizer.optimize(s,sota,N,M,P,SearchBudget.EXACT,null).rank(); }
  public static void main(String[] a) throws Exception {
    String path = a.length>0?a[0]:"src/main/resources/schemes/known/section2/2x2x2-r7-strassen-db11bcc.json";
    NonCubicBilinearAlgorithm base = SchemeIO.read(new File(path));
    SotaResolver sota = Recombination.catalogResolver(Algebra.nonCommutative(Field.R));
    // AXIS_FLIP pool (what live search uses) + s3 axis perms (also in live pool)
    List<SchemeSupports> flip=new ArrayList<>();
    for (var pv: SymmetryTransforms.internalOrbitWithPerms(base, SymmetryTransforms.InternalOrbitMode.AXIS_FLIP, 8))
      for (NonCubicBilinearAlgorithm o: SymmetryTransforms.s3Orbit(pv.alg())) flip.add(SchemeSupports.extract(o));
    // FULL GL frontier members + s3 perms
    int md=Math.max(base.n,Math.max(base.m,base.p));
    Result orbit = md<=3? RecombinationMultisetOrbit.enumerate(base,2): RecombinationMultisetOrbit.enumerateSampled(base,80_000,2);
    List<SchemeSupports> gl=new ArrayList<>();
    for (String key: orbit.dominanceFrontier()){
      int[][][] xyz=orbit.representativeTransforms.get(key); if(xyz==null)continue;
      NonCubicBilinearAlgorithm mem=RecombinationMultisetOrbit.materialise(base,xyz[0],xyz[1],xyz[2]);
      for (NonCubicBilinearAlgorithm o: SymmetryTransforms.s3Orbit(mem)) gl.add(SchemeSupports.extract(o));
    }
    System.out.printf("base ⟨%d,%d,%d⟩=%d  AXIS_FLIP-supports=%d  GL-frontier-supports=%d%n%n",base.n,base.m,base.p,base.r,flip.size(),gl.size());
    int wins=0,tested=0; int lo=base.n*3, hi=base.n*9;
    for (int N=lo;N<=hi;N++) for(int M=lo;M<=hi;M++) for(int P=lo;P<=hi;P++){
      long bfl=Long.MAX_VALUE; for(SchemeSupports s:flip) bfl=Math.min(bfl,bnb(s,sota,N,M,P));
      long bgl=Long.MAX_VALUE; for(SchemeSupports s:gl) bgl=Math.min(bgl,bnb(s,sota,N,M,P));
      tested++;
      if(bgl<bfl){wins++; System.out.printf("  ⟨%d,%d,%d⟩  axisflip=%d  GLfrontier=%d  <<< GL beats AXIS_FLIP by %d%n",N,M,P,bfl,bgl,bfl-bgl);}
    }
    System.out.printf("%nGL-frontier strictly beats AXIS_FLIP pool: %d / %d targets%n",wins,tested);
  }
}
