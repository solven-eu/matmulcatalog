package eu.solven.matmul.docs.explore;
import java.io.File;
import eu.solven.matmul.NonCubicBilinearAlgorithm; import eu.solven.matmul.SymmetryTransforms;
import eu.solven.matmul.algebra.*; import eu.solven.matmul.catalog.SchemeIO;
import eu.solven.matmul.recombination.*;
import eu.solven.matmul.recombination.AnalyticalMaskSearch.SchemeSupports;
import eu.solven.matmul.recombination.Recombination.SotaResolver;
import eu.solven.matmul.search.SearchBudget;
/** Strassen native vs Winograd native (same GL orbit), B&B + s3 orientations, at UNBALANCED cubics. */
public final class ProbeWinoUnbalanced {
  static long best(NonCubicBilinearAlgorithm b, SotaResolver s, int N,int M,int P){
    long r=Long.MAX_VALUE;
    for (NonCubicBilinearAlgorithm o: SymmetryTransforms.s3Orbit(b))
      r=Math.min(r, AllocationOptimizer.optimize(SchemeSupports.extract(o),s,N,M,P,SearchBudget.EXACT,null).rank());
    return r;
  }
  public static void main(String[] a) throws Exception {
    SotaResolver sota = Recombination.catalogResolver(Algebra.nonCommutative(Field.R));
    NonCubicBilinearAlgorithm str = SchemeIO.read(new File("src/main/resources/schemes/known/section2/2x2x2-r7-strassen-db11bcc.json"));
    NonCubicBilinearAlgorithm win = SchemeIO.read(new File("src/main/resources/schemes/known/section2/2x2x2-r7-winograd_1971-511df05.json"));
    int[][] T={{17,17,17},{15,15,15},{13,13,13},{11,11,11},{9,9,9},{17,16,16},{15,16,17},{13,14,15},{9,8,8},{11,10,10}};
    System.out.printf("%-12s %10s %10s %s%n","target","strassen","winograd","winner");
    for (int[] t: T){
      long rs=best(str,sota,t[0],t[1],t[2]), rw=best(win,sota,t[0],t[1],t[2]);
      String w = rw<rs? "WINOGRAD by "+(rs-rw) : rs<rw? "strassen by "+(rw-rs):"tie";
      System.out.printf("⟨%d,%d,%d⟩%s %10d %10d  %s%n",t[0],t[1],t[2]," ".repeat(Math.max(0,5-(""+t[0]+t[1]+t[2]).length())),rs,rw,w);
    }
  }
}
