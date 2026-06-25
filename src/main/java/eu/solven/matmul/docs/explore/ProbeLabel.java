package eu.solven.matmul.docs.explore;
import java.util.*;
import eu.solven.matmul.algebra.*;
import eu.solven.matmul.recombination.*;
import eu.solven.matmul.recombination.BlockSplitSearch.NamedBase;
public final class ProbeLabel {
  public static void main(String[] a){
    var sota=Recombination.catalogResolver(Algebra.nonCommutative(Field.R));
    List<NamedBase> pool=new ArrayList<>();
    for (NamedBase nb:BlockSplitSearch.defaultPool()){var b=nb.base(); if(Math.max(b.n,Math.max(b.m,b.p))<=4)pool.add(nb);}
    for (int[] t: new int[][]{{16,16,16},{20,20,20},{24,24,24}}){
      var best=BlockSplitSearch.findBestStrategy(t[0],t[1],t[2],pool,sota,false);
      System.out.printf("⟨%d,%d,%d⟩ = %s%n", t[0],t[1],t[2], best.map(s->s.rank()+"  ["+s.label()+"]").orElse("none"));
      // also: what does catalog directly say for this shape, and for its halves?
      System.out.printf("   catalog R⟨%d,%d,%d⟩=%d  R⟨%d,%d,%d⟩=%d%n", t[0],t[1],t[2], sota.getRank(t[0],t[1],t[2]),
        t[0]/2,t[1]/2,t[2]/2, sota.getRank(t[0]/2,t[1]/2,t[2]/2));
    }
  }
}
