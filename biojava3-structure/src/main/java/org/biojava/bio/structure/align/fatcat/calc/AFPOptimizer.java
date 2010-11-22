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

import java.util.List;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.model.AFP;
import org.biojava.bio.structure.align.model.AFPChain;

public class AFPOptimizer
{

   public static final boolean debug = FatCatAligner.debug;

   /**
 //optimize the alignment by dynamic programming
    */
   public static void optimizeAln(FatCatParameters params, AFPChain afpChain,Atom[] ca1, Atom[] ca2) throws StructureException
   {

      int minLen = afpChain.getMinLen();
      int fragLen = params.getFragLen();


      long optStart = System.currentTimeMillis();
      int     i, a, k, p1, p2, bk, b1, b2, e1, e2, a1, a2;
      int     iniLen;
      int[][]     iniSet = new int[2][minLen];
      int     maxi = 100;

      int[][][] optAln = afpChain.getOptAln();
      int[] optLen = afpChain.getOptLen();

      int maxTra = params.getMaxTra();
      double[] optRmsd = afpChain.getOptRmsd();

      int blockNum = afpChain.getBlockNum();

      if(optAln == null)      {        
         optAln     = new int[maxTra+1][2][minLen];
         afpChain.setOptAln(optAln);
         optLen     = new int[maxTra+1];
         afpChain.setOptLen(optLen);
         optRmsd    = new double[maxTra+1];
         afpChain.setOptRmsd(optRmsd);
      }

      List<AFP> afpSet = afpChain.getAfpSet();

      int optLength         = afpChain.getOptLength();
      int[] afpChainList    = afpChain.getAfpChainList();
      int[] block2Afp       = afpChain.getBlock2Afp();
      int[] blockSize       = afpChain.getBlockSize();


      if (debug)
         System.out.println("AFPOptimizer got blockNum: " +blockNum);
      //optimize each alignment defined by a block           
      b1 = b2 = e1 = e2 = optLength = 0;
      for(bk = 0; bk < blockNum; bk ++)       {
         //initial aligned position
         iniLen = 0;
         if(bk > 0)      {
            b1 = e1;
            b2 = e2;
         }
         if(bk < blockNum - 1)   {
            a1 = afpChainList[block2Afp[bk] + blockSize[bk] - 1]; //the last AFP in current block
            a2 = afpChainList[block2Afp[bk + 1]];  //the first AFP in next block
            e1 = (afpSet.get(a1).getP1() + fragLen +  afpSet.get(a2).getP1()) / 2;
            e2 = (afpSet.get(a1).getP2() + fragLen +  afpSet.get(a2).getP2()) / 2;
         } //use the middle point of the current and next AFPs. old (starting point of next AFP)
         else    {
            e1 = ca1.length;
            e2 = ca2.length;
         }


         for(i = block2Afp[bk]; i < block2Afp[bk] + blockSize[bk]; i ++) {
            a = afpChainList[i];
            p1 = afpSet.get(a).getP1();
            p2 = afpSet.get(a).getP2();
            for(k = 0; k < afpSet.get(a).getFragLen(); k ++)     {
               iniSet[0][iniLen] = p1 + k - b1; //note -b1
               iniSet[1][iniLen] = p2 + k - b2; //note -b2
               iniLen ++;

            }
         }
         //optimize the align by dynamic programming & constraint the optimization region
         if(debug) {
            System.err.println(String.format("optimize block %d (%d afp), region %d-%d(len %d), %d-%d(len %d)\n",
                  bk, blockSize[bk], b1, e1, e1-b1, b2, e2, e2-b2));
            
            System.err.println(" initial alignment Length: " + iniLen );
         }

         StructureAlignmentOptimizer opt = new StructureAlignmentOptimizer(b1,e1, ca1, b2,e2, ca2, iniLen, iniSet);
         opt.runOptimization(maxi);
         optRmsd[bk] = opt.optimizeResult(optLen,bk,optAln[bk]);

         //System.out.println(optRmsd[bk]);         
         // SALNOPT *opt = new SALNOPT(e1-b1, &pro1->caCod[3 * b1], e2-b2, &pro2->caCod[3 * b2], iniLen, iniSet, maxi);
         // optRmsd[bk] = opt->OptimizeResult(&optLen[bk], optAln[bk]);

         if(debug)
            System.out.println(String.format(" optimized len=%d, rmsd %f\n", optLen[bk], optRmsd[bk]));

         for(i = 0; i < optLen[bk]; i ++)        {
            optAln[bk][0][i] += b1; //restore the position
            optAln[bk][1][i] += b2; //restore the position
         }
         optLength += optLen[bk];


      }

      long optEnd = System.currentTimeMillis();
      if(debug)       System.out.println("complete AlignOpt " + (optEnd-optStart) +"\n");

      afpChain.setBlockNum(blockNum);
      afpChain.setOptLength(optLength);
      afpChain.setAfpChainList(afpChainList);
      afpChain.setBlock2Afp(block2Afp);
      afpChain.setBlockSize(blockSize);


   }

   /**
    * get the afp list and residue list for each block
    */

   public static void blockInfo(AFPChain afpChain)
   {
      int     i, j, k, a, n;

      int blockNum = afpChain.getBlockNum();

      int[] blockSize =afpChain.getBlockSize();
      int[] afpChainList = afpChain.getAfpChainList();
      int[] block2Afp = afpChain.getBlock2Afp();
      int[][][]blockResList = afpChain.getBlockResList();

      List<AFP>afpSet = afpChain.getAfpSet();
      int[] blockResSize = afpChain.getBlockResSize();

      for(i = 0; i < blockNum; i ++)  {
         n = 0;
         for(j = 0; j < blockSize[i]; j ++)      {
            //the index in afpChainList, not in the whole afp set
            a = afpChainList[block2Afp[i] + j];
            for(k = 0; k < afpSet.get(a).getFragLen(); k ++)     {
               blockResList[i][0][n] = afpSet.get(a).getP1() + k;
               blockResList[i][1][n] = afpSet.get(a).getP2() + k;
               n ++;
            }
         }
         blockResSize[i] = n;
      }

      afpChain.setBlockResSize(blockResSize);
      afpChain.setBlockSize(blockSize);
      afpChain.setAfpChainList(afpChainList);
      afpChain.setBlock2Afp(block2Afp);
      afpChain.setBlockResList(blockResList);
   }

   /**
    * to update the chaining score after block delete and merge processed
    * the blockScore value is important for significance evaluation
    */
   public static void updateScore(FatCatParameters params, AFPChain afpChain)
   {
      int     i, j, bknow, bkold, g1, g2;


      afpChain.setConn(0d);
      afpChain.setDVar(0d);      

      int blockNum = afpChain.getBlockNum();
      int alignScoreUpdate = 0;
      double[] blockScore = afpChain.getBlockScore();
      int[] blockGap = afpChain.getBlockGap();
      int[] blockSize =afpChain.getBlockSize();
      int[] afpChainList = afpChain.getAfpChainList();
      List<AFP>afpSet = afpChain.getAfpSet();
      int[] block2Afp = afpChain.getBlock2Afp();

      double torsionPenalty = params.getTorsionPenalty();


      bkold = 0;
      for(i = 0; i < blockNum; i ++)  {
         blockScore[i] = 0;
         blockGap[i] = 0;
         for(j = 0; j < blockSize[i]; j ++)      {
            bknow = afpChainList[block2Afp[i] + j];
            if(j == 0)      {
               blockScore[i] = afpSet.get(bknow).getScore();
            }
            else    {
               AFPChainer.afpPairConn(bkold, bknow, params, afpChain); //note: j, i
               Double conn = afpChain.getConn();
               blockScore[i] += afpSet.get(bknow).getScore() + conn;
               g1 = afpSet.get(bknow).getP1() - afpSet.get(bkold).getP1() - afpSet.get(bkold).getFragLen();
               g2 = afpSet.get(bknow).getP2() - afpSet.get(bkold).getP2() - afpSet.get(bkold).getFragLen();
               blockGap[i] += (g1 > g2)?g1:g2;
            }
            bkold = bknow;
         }
         alignScoreUpdate += blockScore[i];
      }
      if(blockNum >= 2)       {
         alignScoreUpdate += (double)(blockNum - 1) * torsionPenalty;
      }

      afpChain.setBlockGap(blockGap);
      afpChain.setAlignScoreUpdate(alignScoreUpdate);
      afpChain.setBlockScore(blockScore);
      afpChain.setBlockSize(blockSize);
      afpChain.setAfpChainList(afpChainList);
      afpChain.setBlock2Afp(block2Afp);
   }



}
