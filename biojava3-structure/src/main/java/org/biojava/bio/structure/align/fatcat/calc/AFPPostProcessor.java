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
import org.biojava.bio.structure.align.model.AFP;
import org.biojava.bio.structure.align.model.AFPChain;

/** does post processing after alignment chaingin
 * 
 * @author Andreas Prlic
 *
 */
public class AFPPostProcessor
{

   public static final boolean debug = FatCatAligner.debug;

   public static void postProcess(FatCatParameters params, AFPChain afpChain,Atom[] ca1, Atom[] ca2){

      int blockNum = afpChain.getBlockNum();
      afpChain.setBlockNumIni(blockNum);
      //PostProcess of chaining result

      afpChain.setBlockNumIni(blockNum);

      //split blocks (introduce twists) with high RMSD
      splitBlock(params,afpChain, ca1,ca2);
      blockNum = afpChain.getBlockNum();
      afpChain.setBlockNumSpt( blockNum);

      if ( debug){
         System.err.println("AFPPOstProcessor: postProcess blocknum = blocknumSpt:" + blockNum);
      }
      
      //redo: merge blocks with similar transformations & remove small blocks
      //if(blockNum >= 2)     ClustBlock();

      deleteBlock(params,afpChain,ca1,ca2);
      mergeBlock(params,afpChain,ca1,ca2);

      afpChain.setBlockNumClu(afpChain.getBlockNum());


   }

   /**
    * in some special cases, there is no maginificent twists in the
   final chaining result; however, their rmsd (original and after
   optimizing) are very large. Therefore, a post-process is taken
   to split the blocks further at the ralative bad connections (
   with relative high distance variation)
   to be tested:
     split or not according to optimized or initial chaining???
    */

   private static void splitBlock(FatCatParameters params, AFPChain afpChain, Atom[] ca1, Atom[] ca2)
   {
      if ( debug)
         System.err.println("AFPPostProcessor: splitBlock");
      int     i, a, bk, cut;
      double  maxs, maxt;
      int blockNum = afpChain.getBlockNum();
      int maxTra = params.getMaxTra();
      double badRmsd = params.getBadRmsd();

      int     blockNum0 = blockNum;

      double[] blockRmsd = afpChain.getBlockRmsd();
      int[] blockSize = afpChain.getBlockSize();
      int[] block2Afp = afpChain.getBlock2Afp();
      double[] afpChainTwiList = afpChain.getAfpChainTwiList();

      bk = 0;
      while(blockNum < maxTra + 1)    {
         maxs = 0;
         for(i = 0; i < blockNum; i ++)   {
            if(blockRmsd[i] > maxs && blockSize[i] > 2) { //according to the optimized alignment
               maxs = blockRmsd[i];
               bk = i;
            } //!(Note: optRmsd, not blockRmsd, according to the optimized alignment
         }
         if(maxs < badRmsd)      break;
         maxt = 0;
         cut = 0;
         for(i = 1; i < blockSize[bk]; i ++)     {
            a = i + block2Afp[bk];            
            if(afpChainTwiList[a] > maxt)   {
               maxt = afpChainTwiList[a];
               cut = i;
               
            }
         }
         if(debug)
            System.out.println(String.format("block %d original size %d rmsd %.3f maxt %.2f cut at %d\n", bk, blockSize[bk], maxs, maxt, cut));
         for(i = blockNum - 1; i > bk; i --)     {
            block2Afp[i + 1] = block2Afp[i];
            blockSize[i + 1] = blockSize[i];
            blockRmsd[i + 1] = blockRmsd[i];
         } //update block information
         block2Afp[bk + 1] = cut + block2Afp[bk];
         blockSize[bk + 1] = blockSize[bk] - cut;
         blockSize[bk] = cut;

         if(debug)
            System.out.println(String.format("  split into %d and %d sizes\n", blockSize[bk], blockSize[bk + 1]));


         int[] afpChainList = afpChain.getAfpChainList();
         //int[] subrange1    = getSubrange(afpChainList, block2Afp[bk + 1] );                           
         blockRmsd[bk + 1]  = AFPChainer.calAfpRmsd(blockSize[bk + 1],  afpChainList, block2Afp[bk + 1] , afpChain, ca1, ca2);

         //int[] subrange2    = getSubrange(afpChainList, block2Afp[bk] );   
         blockRmsd[bk]      = AFPChainer.calAfpRmsd(blockSize[bk],      afpChainList, block2Afp[bk], afpChain, ca1, ca2);

         //split a block at the biggest position
         blockNum ++;
         afpChain.setAfpChainList(afpChainList);
      }
      if(blockNum - blockNum0 > 0)    {
         if(debug)
            System.out.println(String.format("Split %d times:\n", blockNum - blockNum0));
         for(i = 0; i < blockNum; i ++)  {
            if(debug)
               System.out.println(String.format("  block %d size %d from %d rmsd %.3f\n", i, blockSize[i], block2Afp[i], blockRmsd[i]));
         }
      }

      
      afpChain.setBlockNum(blockNum);
      afpChain.setBlockSize(blockSize);
      afpChain.setBlockRmsd(blockRmsd);
      afpChain.setBlock2Afp(block2Afp);
     

   }

   /**
    * remove the artifical small rigid-body superimpose in the middle
    clust the similar superimpositions (caused by the small flexible
    region, which is detected as a seperate rigid superimposing region by adding
    two twists before and after it(artifically!)
    one possible solution: allowing long enough loops in the chaining process,
    which however increase the calculation complexity
    */
   private static void deleteBlock(FatCatParameters params, AFPChain afpChain,Atom[] ca1, Atom[] ca2)
   {
      int blockNum = afpChain.getBlockNum();
      List<AFP> afpSet = afpChain.getAfpSet();

      int[] afpChainList =afpChain.getAfpChainList();



      int[] block2Afp = afpChain.getBlock2Afp();
      int[] blockSize = afpChain.getBlockSize();

      double[] blockRmsd = afpChain.getBlockRmsd();

      int fragLen = params.getFragLen();

      //remove those blocks (both in terminals and in the middle) with only a AFP
      //but still keep those small blocks spaning large regions
      if(blockNum <= 1)       return;
      int     blockNumOld = blockNum;
      int     i, j, b1, b2, e1, e2, len;
      e1 = e2 = 0;
      for(i = 0; i < blockNum; i ++) {
         b1 = e1;
         b2 = e2;
         if(i < blockNum - 1)    {
            e1 = afpSet.get(afpChainList[block2Afp[i + 1]]).getP1();
            e2 = afpSet.get(afpChainList[block2Afp[i + 1]]).getP2();
         }
         else    {
            e1 = ca1.length;
            e2 = ca2.length;
         }
         if(blockSize[i] > 1)    continue;
         len = (e1 - b1) < (e2 - b2)?(e1 - b1):(e2 - b2);
         //if(i == blockNum - 1) blockNum --;
         if(len < 2 * fragLen)   {
            for(j = i; j < blockNum - 1; j ++)      {
               blockRmsd[j] = blockRmsd[j + 1];
               blockSize[j] = blockSize[j + 1];
               block2Afp[j] = block2Afp[j + 1];
            }
            blockNum --;
            i --;
         } //delete a block
      }
      if(blockNumOld > blockNum)
         if(debug)
            System.out.println(
                  String.format("Delete %d small blocks\n", blockNumOld - blockNum)
            );


      if (debug)
         System.err.println("deleteBlock: end blockNum:"+ blockNum);
      afpChain.setBlock2Afp(block2Afp);
      afpChain.setBlockSize(blockSize);
      afpChain.setAfpChainList(afpChainList);
      afpChain.setBlockNum(blockNum);
      afpChain.setBlockRmsd(blockRmsd);
   }


   /**
 //Merge consecutive blocks with similar transformation
    */
   private static  void mergeBlock(FatCatParameters params, AFPChain afpChain,Atom[] ca1,Atom[] ca2)
   {

      int blockNum = afpChain.getBlockNum();
      double badRmsd = params.getBadRmsd();

      int[] block2Afp = afpChain.getBlock2Afp();
      int[] blockSize = afpChain.getBlockSize();

      double[] blockRmsd = afpChain.getBlockRmsd();

      int afpChainTwiNum = afpChain.getAfpChainTwiNum();

      //clustering the neighbor blocks if their transformations are similar
      int     i, j, b1, b2, minb1, minb2;
      double  minrmsd;
      int     merge = 0;
      int     blockNumOld = blockNum;
      double[][]  rmsdlist = null;
      if(blockNum > 1)        {
         rmsdlist = new double[blockNumOld][blockNumOld];
         for(b1 = 0; b1 < blockNum - 1; b1 ++)   {
            for(b2 = b1 + 1; b2 < blockNum; b2 ++)  {
               rmsdlist[b1][b2] = combineRmsd(b1, b2,afpChain,ca1,ca2);
            }
         }
      }
      minb1 = 0;
      while(blockNum > 1)     {
         minrmsd = 1000;
         for(i = 0; i < blockNum - 1; i ++)      {
            j = i + 1; //only consider neighbor blocks
            if(minrmsd > rmsdlist[i][j])    {
               minrmsd = rmsdlist[i][j];
               minb1 = i;
            }
         }
         minb2 = minb1 + 1; //merge those most similar blocks
         //maxrmsd = (blockRmsd[minb1] > blockRmsd[minb2])?blockRmsd[minb1]:blockRmsd[minb2];
         if(minrmsd < badRmsd)   {
            if(debug)
               System.out.println(String.format("merge block %d (rmsd %.3f) and %d (rmsd %.3f), total rmsd %.3f\n",
                     minb1, blockRmsd[minb1], minb2, blockRmsd[minb2], minrmsd));
            blockSize[minb1] += blockSize[minb2];
            blockRmsd[minb1] = minrmsd;
            for(i = minb2; i < blockNum - 1; i ++)  {
               block2Afp[i] = block2Afp[i + 1];
               blockSize[i] = blockSize[i + 1];
               blockRmsd[i] = blockRmsd[i + 1];
            } //update block information
            afpChainTwiNum --;
            blockNum --;
            for(b1 = 0; b1 < blockNum - 1; b1 ++)   {
               for(b2 = b1 + 1; b2 < blockNum; b2 ++) {
                  if(b1 == minb1 || b2 == minb1)  {
                     rmsdlist[b1][b2] = combineRmsd(b1, b2, afpChain,ca1,ca2);
                  }
                  else if(b2 < minb1)     continue;
                  else if(b1 < minb1)     {
                     rmsdlist[b1][b2] = rmsdlist[b1][b2 + 1];
                  }
                  else    {
                     rmsdlist[b1][b2] = rmsdlist[b1 + 1][b2 + 1];
                  }
               }
            } //update the rmsdlist
            merge ++;
         } //merge two blocks
         else if(minrmsd >= 100) break;
         else    {
            rmsdlist[minb1][minb2] += 100;
         } //not merge, modify the rmsdlist so that this combination is not considered in next iteration
      }

      if(merge > 0)       {
         if(debug)
            System.out.println(String.format("Merge %d blocks, remaining %d blocks\n", merge, blockNum));
      }

      if (debug){
         System.err.println("AFPPostProcessor: mergeBlock end blocknum:" + blockNum);
      }
      afpChain.setBlock2Afp(block2Afp);
      afpChain.setBlockSize(blockSize);      
      afpChain.setBlockNum(blockNum);
      afpChain.setBlockRmsd(blockRmsd);
      afpChain.setAfpChainTwiNum(afpChainTwiNum);
   }


   /**
   return the rmsd of two blocks
    */
   private static double combineRmsd(int b1, int b2, AFPChain afpChain,Atom[] ca1,Atom[] ca2)
   {
      int     i;
      int     afpn = 0;

      int[] afpChainList =afpChain.getAfpChainList();

      int[] block2Afp = afpChain.getBlock2Afp();
      int[] blockSize = afpChain.getBlockSize();


      int[]   list = new int[blockSize[b1]+blockSize[b2]];
      for(i = block2Afp[b1]; i < block2Afp[b1] + blockSize[b1]; i ++) {
         list[afpn ++] = afpChainList[i];
      }
      for(i = block2Afp[b2]; i < block2Afp[b2] + blockSize[b2]; i ++) {
         list[afpn ++] = afpChainList[i];
      }
      double  rmsd = AFPChainer.calAfpRmsd(afpn, list,0, afpChain,ca1,ca2);

      afpChain.setBlock2Afp(block2Afp);
      afpChain.setBlockSize(blockSize);  
      afpChain.setAfpChainList(afpChainList);

      return rmsd;
   }


}
