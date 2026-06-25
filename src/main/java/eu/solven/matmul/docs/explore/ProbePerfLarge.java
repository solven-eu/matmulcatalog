package eu.solven.matmul.docs.explore;
import java.io.File; import java.nio.file.Path; import java.util.*;
import eu.solven.matmul.algebra.*;
import eu.solven.matmul.recombination.*;
import eu.solven.matmul.recombination.BlockSplitSearch.NamedBase;
import eu.solven.matmul.recombination.Recombination.SotaResolver;
public final class ProbePerfLarge {
  public static void main(String[] a) throws Exception {
    SotaResolver sota = Recombination.catalogResolver(Algebra.nonCommutative(Field.R));
    // pool restricted to baseDim<=4
    List<NamedBase> pool = new ArrayList<>();
    for (NamedBase nb : BlockSplitSearch.defaultPool()) { var b=nb.base();
      if (Math.max(b.n,Math.max(b.m,b.p))<=4) pool.add(nb); }
    System.out.println("pool (baseDim<=4): "+pool.size()+" base-variants");
    // frontier sidecars (baseDim<=4)
    File dir=new File("src/main/resources/frontiers");
    List<RecombFrontierIO.Loaded> loaded=new ArrayList<>();
    for (File f: dir.listFiles((d,n)->n.endsWith(".json"))) {
      var L=RecombFrontierIO.read(Path.of(f.getPath()));
      if (Math.max(L.dims()[0],Math.max(L.dims()[1],L.dims()[2]))<=4) loaded.add(L);
    }
    System.out.println("frontier sidecars (baseDim<=4): "+loaded.size());
    int[][] shapes={{12,12,12},{16,16,16},{20,20,20},{24,24,24},{18,19,20}};
    // warmup
    BlockSplitSearch.findBestStrategy(8,8,8,pool,sota,false);
    for (var c:loaded) FrontierRecombination.bestRankOverOrientations(8,8,8,c.dims(),c.frontier(),sota);
    for (int[] t:shapes){
      long s0=System.nanoTime();
      long sr=BlockSplitSearch.findBestStrategy(t[0],t[1],t[2],pool,sota,false).map(s->s.rank()).orElse(-1L);
      long st=System.nanoTime()-s0;
      long f0=System.nanoTime();
      long fr=Long.MAX_VALUE;
      for (var c:loaded) fr=Math.min(fr,FrontierRecombination.bestRankOverOrientations(t[0],t[1],t[2],c.dims(),c.frontier(),sota));
      long ft=System.nanoTime()-f0;
      System.out.printf("⟨%d,%d,%d⟩  search=%d (%.0fms)  frontier=%d (%.0fms)  speedup=%.1fx  Δrank=%+d%n",
        t[0],t[1],t[2], sr, st/1e6, fr, ft/1e6, st/(double)ft, fr-sr);
    }
  }
}
