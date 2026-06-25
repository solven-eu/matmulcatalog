package eu.solven.matmul.docs.migrate;
import java.io.File; import java.util.*;
import eu.solven.matmul.catalog.*; import eu.solven.matmul.algebra.Field;
/** Re-stamp a scheme's fields[] by recomputing fieldNamesFromLineage over Q (the rational base
 *  field of our pipeline) — fixes R-floor under-claims (e.g. [R,C] → [Q,R,C]). Args: file paths. */
public final class ReStampFieldsQ {
  public static void main(String[] args) throws Exception {
    FieldAwareLookup q=new FieldAwareLookup(Field.Q);
    for (String p: args){
      File f=new File(p);
      var root=SchemeIO.parseJson(f);
      var lin=SchemeIO.readLineage(root).orElse(null);
      if (lin==null){ System.out.println("no lineage: "+f.getName()); continue; }
      List<String> before=SchemeIO.fieldTags(root);
      List<String> after=q.fieldNamesFromLineage(lin);
      if (after.isEmpty()){ System.out.println("empty fields computed, skip: "+f.getName()); continue; }
      SchemeIO.updateFields(f, Map.of("fields", after), List.of(), true);
      System.out.printf("%s : %s -> %s%n", f.getName(), before, after);
    }
  }
}
