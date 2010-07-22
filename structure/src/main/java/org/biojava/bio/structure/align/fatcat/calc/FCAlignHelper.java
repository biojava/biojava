/* This class is based on the original FATCAT implementation by
 * <pre>
 * Yuzhen Ye & Adam Godzik (2003)
 * Flexible structure alignment by chaining aligned fragment pairs allowing twists.
 * Bioinformatics vol.19 suppl. 2. ii246-ii255.
 * http://www.ncbi.nlm.nih.gov/pubmed/14534198
 * </pre>
 * 
 * Thanks to Yuzhen Ye and A. Godzik for granting permission to freely use and redistribute this code.
 *  
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 *
 * Created on Jun 17, 2009
 * Created by Andreas Prlic - RCSB PDB 
 * 
 */

package org.biojava.bio.structure.align.fatcat.calc;

public class FCAlignHelper
{
   
   int     M; //length of protein 1
   int     N; //length of protein 2
   double  g; //gap-create
   double  h; //gap-extend
   double  m; //g + h
   double[][]  sij;
   char[][]    trace; //trace-record
   char[][]    etrace; //trace-record
   char[][]    dtrace; //trace-record
   int     B1; //beginning position of protein 1 in alignment
   int     B2; //beginning position of protein 2 in alignment
   int     E1; //end position of protein 1 in alignment
   int     E2; //end position of protein 2 in alignment
   double  alignScore;
   double  identity;
   double  similarity;
   int[]   sapp;
   int[]     sapp0;
   int sappPos;
   int     last;

   char[]    seq1;
   char[]    seq2;
   char[]    aln1;
   char[]    aln2;
   char[]    mark;

   /** do an alignment given the provided matrix sij0
    * 
    * @param sij0 - the matrix to perform the calculations on.
    * @param M0
    * @param N0
    * @param g0
    * @param h0
    */
   public FCAlignHelper(double[][] sij0, int M0, int N0, double g0, double h0){
      init(M0, N0, g0, h0);
      int     i, j;
      for(i = 0; i < M; i ++) {
         for(j = 0; j < N; j ++){ 
            sij[i][j] = sij0[i][j];
            //System.out.println(i+"-"+j+":" +sij[i][j]);
         }         
      }
      
      doAlign();
   }


   private void init(int M0, int N0, double g0, double h0)
   {
      M = M0;
      N = N0;
      g = g0; //gap-create
      h = h0; //gap-extend
      m = g + h; //gap-create + gap-extend
      trace  = new char[M+1][N+1];
      etrace = new char[M+1][N+1];
      dtrace = new char[M+1][N+1];
      B1 = B2 = E1 = E2 = 0;
      alignScore = 0;
      last = 0;
      sapp = new int[M+N];
      sapp0 = sapp;
      sappPos = 0;
      sij = new double[M][N];
      seq1 = new char[M+1];
      seq2 = new char[N+1];
      aln1 = aln2 = mark = null;
      identity = similarity = 0;
   }

   //trace-back strategy
   //affine-gap penalty
   //local-model

   private void doAlign(){


      int     i, j;
      double  s, e, c, d, wa;
      double[]  CC = new double[N+1]; //note N + 1
      double[]  DD = new double[N+1];
      double  maxs = -100;
      char    trace_e, trace_d;

      //forward-phase
      CC[0] = 0;
      for(j = 1; j <= N; j ++)        {
         CC[j] = 0;
         DD[j] = -g;
      } //local-alignment, no terminal penalty
      for(i = 1; i <= M; i ++)        {
         CC[0] = c = s = 0;
         e = -g;
         for(j = 1; j <= N; j ++)        {
            trace_e = 'e';
            if ((c =   c   - m) > (e =   e   - h)) {
               e = c;  trace_e = 'E';
            }//insertion
            trace_d = 'd';
            if ((c = CC[j] - m) > (d = DD[j] - h)) {
               d = c;  trace_d = 'D';
            }//deletion
            //ie   CC[j]==CC[i-1][j]   DD[j]==DD[i-1][j]
            wa = sij[i - 1][j - 1]; //note i - 1, j - 1
            c = s + wa; //s==CC[i-1][j-1]
            trace[i][j] = 's';
            if (e > c) {
               c = e;
               trace[i][j] = trace_e;
            }
            if (d > c) {
               c = d;
               trace[i][j] = trace_d;
            }
            etrace[i][j] = trace_e;
            dtrace[i][j] = trace_d;
            s = CC[j]; //important for next replace
            CC[j] = c; //CC[i][j]
            DD[j] = d; //DD[i][j]
            if(c < 0)       {
               CC[j] = 0;
               DD[j] = -g;
               c = 0;
               e = -g;
               trace[i][j] = '0';
            } //local-N
            if(c > maxs)    {
               E1 = i;
               E2 = j;
               maxs = c;
            } //local-C
         }
      }
      alignScore = maxs;
      //printf("alignment score %f\n", alignScore);


      //trace-back
      if(trace[E1][E2] != 's')        {
         throw new RuntimeException("FCAlignHelper encoutered Exception: Not ending with substitution");

      }
      //Trace(maxs, E1, E2);
      trace('s', E1, E2);
      //printf("B1 %d B2 %d, E1 %d E2 %d\n", B1, B2, E1, E2);

      //check-alignment
      checkAlign();
      
      
   }


   /**
 trace-back, recorded in sapp, wrong method!
    */

   private void trace(char mod, int i, int j)
   {
      if(mod == '0' || i <= 0 || j <= 0)      {
         B1 = i + 1;
         B2 = j + 1;
      }
      if(mod == 's')  {
         trace(trace[i - 1][j - 1], i - 1, j - 1);
         rep();
      }
      else if(mod == 'D')     {
         trace(trace[i - 1][j], i - 1, j);
         del(1);
      }
      else if(mod == 'd')     {
         trace(dtrace[i - 1][j], i - 1, j);
         del(1);
      }
      else if(mod == 'E')     {
         trace(trace[i][j - 1], i, j - 1);
         ins(1);
      }
      else if(mod == 'e')     {
         trace(etrace[i][j - 1], i, j - 1);
         ins(1);
      }
   }

   //-----------------------------------------------------------------------------
   //record the alignment in sapp
   //deletion, sapp < 0, sequences in i, gaps in j
   //-----------------------------------------------------------------------------
   private void del(int k)
   {
      //if(last < 0)    last = sapp[-1] -= (k);
      //else            last = *sapp++ = -(k);
      
      if(last < 0)    last = sapp[sappPos-1]   -=  (k);
      else            last = sapp[(sappPos++)]  = -(k);
   }

   //Insertion, sapp > 0, gaps in i, sequences in j
   //-----------------------------------------------------------------------------
   private void ins(int k)
   {

      //if(last > 0)    last = sapp[-1] += k;
      //else            last = *sapp++ = (k);
      if(last > 0)    last = sapp[sappPos-1]   +=  k;
      else            last = sapp[(sappPos++)]  = (k);
   }

   //-----------------------------------------------------------------------------
   private void rep()
   {
      
      // last = *sapp++ = 0;
      last = sapp[(sappPos++)] = 0;
   }

   private void checkAlign(){

      if(sapp[0] != 0)        {
         System.err.println(String.format("warn: not a local-alignment result, first operation %d\n", sapp[0]));
      }
      double  sco = checkScore();
      if(Math.abs(sco - alignScore) > 1e-3)       {
         System.err.println(String.format("FCAlignHelper: warn: alignment scores are different %f(check) %f(align)\n", sco, alignScore));
      }
   }

   /**
    *  checkscore - return the score of the alignment stored in sapp
    */

   private double checkScore()
   {
      int     i, j, op, s;
      double  sco;

      sco = 0;
      op = 0;
      s = 0;

      i = B1;
      j = B2;
      while (i <= E1 && j <= E2) {
         op = sapp0[s ++];
         if (op == 0)    {
            sco += sij[i - 1][j - 1];
            //if (debug)
               //System.err.println(String.format("%d-%d %f\n", i - 1, j - 1, sij[i - 1][j - 1]));
            i ++;
            j ++;
         }
         else if (op > 0) {
            sco -= g+op*h;
            j = j+op;
         }
         else {
            sco -= g-op*h;
            i = i-op;
         }
      }
      return(sco);
   }

   /**
    * record the aligned pairs in alignList[][0], alignList[][1];
    * return the number of aligned pairs
    * @param alignList
    * @return the number of aligned pairs
    */
   public int getAlignPos(int[][] alignList)
   {
      int     i = B1;
      int     j = B2;
      int     s = 0;
      int     a = 0;
      int     op;
      while(i <= E1 && j <= E2)       {
         op = sapp0[s ++];
        
         if (op == 0)    {
            alignList[0][a] = i - 1; //i - 1
            alignList[1][a] = j - 1;
            a ++;
            i ++;
            j ++;
         }
         else if (op > 0) {
            j += op;
         }
         else {
            i -= op;
         }
      }
      return a;
   }

}


