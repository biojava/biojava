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
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.AFPTwister;
import org.biojava.bio.structure.align.model.AFP;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AFPAlignmentDisplay;


/** A class that does calculations on an AFPChain
 * 
 * @author Andreas Prlic
 *
 */
public class FatCatAligner 
{

   public static final boolean debug = false;
   public static final boolean printTimeStamps = false;

   AFPChain afpChain ;
   Group[] twistedGroups;
   
   
   
   public AFPChain getAfpChain()
   {
      return afpChain;
   }

  
   public Group[] getTwistedGroups()
   {
	   	  
      return twistedGroups;
   }
   
   public void setTwistedGroups(Group[] twistedGroups)
   {
	   this.twistedGroups = twistedGroups;
   }

   
   public  void align(Atom[] ca1, Atom[] ca2, boolean doRigid, FatCatParameters params) throws StructureException{

	  long tstart = System.currentTimeMillis();
    
      afpChain = new AFPChain();
      afpChain.setCa1Length(ca1.length);
      afpChain.setCa2Length(ca2.length);
 
      

      AFPCalculator.extractAFPChains(params, afpChain,ca1, ca2);

      long cend = System.currentTimeMillis();

      if (printTimeStamps)
         System.out.println("calculation took:" + (cend - tstart) + " ms.");

      
      AFPCalculator.sortAfps(afpChain,ca1,ca2);

      if ( printTimeStamps) {
         long send = System.currentTimeMillis();


         System.out.println("sorting  took:" + (send - cend) + " ms.");
      }
      
      if ( doRigid)
         this.twistedGroups = rChainAfp(params, afpChain,ca1,ca2);

      else {
         this.twistedGroups = chainAfp(params,afpChain,ca1,ca2);
      }
      
     // long start = System.currentTimeMillis();
		long end = System.currentTimeMillis();
		afpChain.setCalculationTime(end-tstart);
		if ( printTimeStamps)
			System.out.println("TOTAL calc time: " + (end -tstart) / 1000.0);
   
   }
   
   
   
   
   /** runs rigid chaining process
   *
   */
  private static Group[] rChainAfp(FatCatParameters params, AFPChain afpChain, Atom[] ca1, Atom[] ca2) throws StructureException{
     params.setMaxTra(0);
     afpChain.setMaxTra(0);
     return chainAfp(params,afpChain,ca1,ca2);
  }

  /**
   * run AFP chaining allowing up to maxTra flexible regions.
   * Input is original coordinates. 
   * 
   */
  private static Group[] chainAfp(FatCatParameters params,AFPChain afpChain, Atom[] ca1, Atom[] ca2) throws StructureException{
     
	// we don;t want to rotate input atoms, do we?
	  Atom[] ca2clone = StructureTools.cloneCAArray(ca2);
	  
     List<AFP> afpSet = afpChain.getAfpSet();
     
     if (debug)
        System.out.println("entering chainAfp");
     int afpNum = afpSet.size();

     if ( afpNum < 1)
        return new Group[0];

     long bgtime = System.currentTimeMillis();
     if(debug)    {
        System.out.println(String.format("total AFP %d\n", afpNum));
     }

     //run AFP chaining

     AFPChainer.doChainAfp(params,afpChain ,ca1,ca2);

     int afpChainLen = afpChain.getAfpChainLen();
     
     if(afpChainLen < 1)     {
        
        afpChain.setShortAlign(true);
        return new Group[0];
     } //very short alignment

     long chaintime = System.currentTimeMillis();
     if(debug)    {

        System.out.println("Afp chaining: time " + (chaintime-bgtime));
     }
     
     // do post processing
     
     AFPPostProcessor.postProcess(params, afpChain,ca1,ca2);
          
     // Optimize the final alignment 
     
     AFPOptimizer.optimizeAln(params, afpChain,ca1,ca2);

     AFPOptimizer.blockInfo( afpChain);

     AFPOptimizer.updateScore(params,afpChain);

     AFPAlignmentDisplay.getAlign(afpChain,ca1,ca2);
                                            
     Group[] twistedPDB = AFPTwister.twistPDB(afpChain, ca1, ca2clone);
     
     SigEva sig =  new SigEva();
     double probability = sig.calSigAll(params, afpChain);
     afpChain.setProbability(probability);
     double normAlignScore = sig.calNS(params,afpChain);
     afpChain.setNormAlignScore(normAlignScore);
    
     /*

  SIGEVA  sig;
  probability = sig.calSigAll(maxTra, sparse, pro1Len, pro2Len, alignScore, totalRmsdOpt, optLength, blockNum - 1);
  normAlignScore = sig.calNS(pro1Len, pro2Len, alignScore, totalRmsdOpt, optLength, blockNum - 1);

      */

     //if(maxTra == 0)       probability = sig.calSigRigid(pro1Len, pro2Len, alignScore, totalRmsdOpt, optLength);
     //else  probability = sig.calSigFlexi(pro1Len, pro2Len, alignScore, totalRmsdOpt, optLength, blockNum - 1);

     if(debug)    {

        long nowtime = System.currentTimeMillis();
        long diff = nowtime - chaintime;
        System.out.println("Alignment optimization: time "+ diff);

        System.out.println("score:      " + afpChain.getAlignScore());
        System.out.println("opt length: " + afpChain.getOptLength());
        System.out.println("opt rmsd:   "+ afpChain.getTotalRmsdOpt());

     }
     return twistedPDB;

  }
   
   



   
}
