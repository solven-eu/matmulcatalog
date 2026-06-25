package eu.solven.matmul.docs.explore;
import eu.solven.matmul.*; import eu.solven.matmul.catalog.*; import eu.solven.matmul.algebra.Field;
public class ProbeOrientedHash { public static void main(String[] a) throws Exception {
  FieldAwareLookup lk=new FieldAwareLookup(Field.R);
  // ref orientedShape, origHash, canonical
  Object[][] refs={{7,3,3,"00460b6"},{10,13,9,"16f3df1"}};
  for (Object[] rf: refs){
    int n=(int)rf[0],m=(int)rf[1],p=(int)rf[2]; String h=(String)rf[3];
    var byCanon=lk.findByHash(n,m,p,h); // resolves via orientAs of canonical
    if (byCanon.isEmpty()){ System.out.println(n+"x"+m+"x"+p+"@"+h+" : findByHash(origHash) EMPTY"); continue; }
    NonCubicBilinearAlgorithm oriented=byCanon.get().alg();
    String orientedHash=SchemeIO.contentHash(oriented);
    System.out.printf("%dx%dx%d : origHash=%s orientedHash=%s  same=%s  shape=⟨%d,%d,%d⟩%n",
      n,m,p,h,orientedHash.substring(0,7), h.equals(orientedHash.substring(0,7)), oriented.n,oriented.m,oriented.p);
    // does findByHash with the ORIENTED hash resolve precisely?
    var byOriented=lk.findByHash(n,m,p,orientedHash);
    System.out.println("   findByHash(orientedHash) resolves="+byOriented.isPresent());
  }
}}
