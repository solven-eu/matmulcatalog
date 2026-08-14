package eu.solven.matmul.docs.explore;
import eu.solven.matmul.NonCubicBilinearAlgorithm;
import eu.solven.matmul.papers.schwartzzwecher2025.TaNew25Construction;
public final class ProbeTaNew25 {
  public static void main(String[] a){
    int n0=4; NonCubicBilinearAlgorithm alg=TaNew25Construction.build(n0);
    double[][] U=alg.denseU(), V=alg.denseV(), W=alg.denseW(); int r=alg.r;
    int bad=0, frac=0;
    // per output-cell (x,y) mismatch counts, split by whether the CORRECT (i1,k2) path is diagonal.
    int[][] perCell=new int[n0][n0];
    long diagOut=0, offOut=0;
    for(int i1=0;i1<n0;i1++)for(int j1=0;j1<n0;j1++)for(int j2=0;j2<n0;j2++)for(int k2=0;k2<n0;k2++)for(int x=0;x<n0;x++)for(int y=0;y<n0;y++){
      double rec=0; for(int l=0;l<r;l++) rec+=U[i1*n0+j1][l]*V[j2*n0+k2][l]*W[x*n0+y][l];
      double exp=(j1==j2&&i1==x&&k2==y)?1:0;
      if(Math.abs(rec-exp)>1e-9){ bad++;
        if(Math.abs(rec-Math.rint(rec))>1e-6) frac++;
        perCell[x][y]++;
        if(x==y) diagOut++; else offOut++;
      }
    }
    System.out.printf("mismatches=%d (non-integer=%d)  outputDiag=%d outputOff=%d%n", bad, frac, diagOut, offOut);
    System.out.println("per output-cell (x,y) mismatch counts:");
    for(int x=0;x<n0;x++){StringBuilder sb=new StringBuilder("  ");for(int y=0;y<n0;y++)sb.append(String.format("%4d",perCell[x][y]));System.out.println(sb);}
  }
}
