package eu.solven.matmul.docs.verify;
import java.io.File; import java.nio.file.*; import java.util.*;
import eu.solven.matmul.catalog.*;
/** Scan all scheme files: readLineage each and report those whose IN-MEMORY lineage contains an
 *  unresolvable Atom("@ref?:Ln") — i.e. a STORED cyclic/dedup-broken lineage (the cycle source). */
public final class FindFallbackRefSources {
  static boolean hasFallback(Lineage.Node n){
    if (n instanceof Lineage.Atom a) return a.ref()!=null && a.ref().startsWith("@ref?:");
    for (Lineage.Node c: Lineage.childrenOf(n)) if (hasFallback(c)) return true; return false;
  }
  public static void main(String[] x) throws Exception {
    int[] scanned={0},hits={0};
    Files.walk(Path.of("src/main/resources/schemes")).filter(p->p.toString().endsWith(".json")).forEach(p->{
      try { var root=SchemeIO.parseJson(p.toFile()); var lin=SchemeIO.readLineage(root).orElse(null);
        scanned[0]++;
        if (lin!=null && hasFallback(lin)){ hits[0]++; System.out.println("CORRUPT-SOURCE "+p.getFileName()); }
      } catch(Exception e){}
    });
    System.out.println("scanned="+scanned[0]+" corrupt-sources="+hits[0]);
  }
}
