package eu.solven.matmul.docs.explore;
import eu.solven.matmul.algebra.*;
import eu.solven.matmul.catalog.FieldAwareLookup;
import eu.solven.matmul.recombination.BlockSplitSearch;
import eu.solven.matmul.recombination.Recombination;
import eu.solven.matmul.recombination.Recombination.SotaResolver;
/** Are the GL-frontier single-base wins also wins vs the FULL pool / catalog SOTA? */
public final class ProbeGlVsFullPool {
  public static void main(String[] a){
    SotaResolver sota = Recombination.catalogResolver(Algebra.nonCommutative(Field.R));
    var pool = BlockSplitSearch.defaultPool();
    FieldAwareLookup look = new FieldAwareLookup(Field.R);
    int[][] gl = {{7,7,7,278},{9,9,9,648},{7,7,10,428},{7,7,11,477},{7,9,9,493},{9,9,7,493},{7,7,9,379},{6,9,9,428}};
    System.out.printf("%-12s %8s %10s %10s %s%n","shape","GLrecomb","fullPool","catalogR","verdict");
    for (int[] t: gl){
      long fp = BlockSplitSearch.findBestStrategy(t[0],t[1],t[2],pool,sota,false).map(s->s.rank()).orElse(Long.MAX_VALUE);
      int cr = look.findRank(t[0],t[1],t[2]);
      String v = t[3]<fp? "GL beats fullPool by "+(fp-t[3]) : t[3]<cr? "GL beats catalog by "+(cr-t[3]) : "subsumed (pool/catalog ≤ GL)";
      System.out.printf("⟨%d,%d,%d⟩%s %8d %10d %10d  %s%n",t[0],t[1],t[2]," ".repeat(Math.max(0,4-(""+t[0]+t[1]+t[2]).length())),t[3],fp,cr,v);
    }
  }
}
