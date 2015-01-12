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

package org.biojava.bio.structure.align.fatcat;

import org.biojava.bio.structure.Atom;


import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.StructureAlignment;
import org.biojava.bio.structure.align.fatcat.calc.FatCatAligner;
import org.biojava.bio.structure.align.fatcat.calc.FatCatParameters;
import org.biojava.bio.structure.align.model.AFPChain;


public class FatCat
{

   /**
    *  version history:
    *  1.1 - Added more parameters to the command line
    *  1.0 - Initial version
    */
   public static final String VERSION = "1.1";

   public static String newline = System.getProperty("line.separator");

   FatCatAligner aligner;




   /** See demo/FatCatDemo.java for an example how to run.
    * 
    * Launch FatCat from command line.
    * 
    * Parameters are:
    * 
    * @param argv
    */
   public static void main(String[] argv){
      FatCatUserArgumentProcessor processor = new FatCatUserArgumentProcessor();
      processor.process(argv);
   }

   @Override
   public String toString(){
      return "JFatCat v. " + VERSION;
   }


   public AFPChain alignRigid(Atom[] ca1, Atom[] ca2) throws StructureException{
      StructureAlignment fatCat = new FatCatRigid();
      return fatCat.align(ca1,ca2);
   }

   public AFPChain alignRigid(Atom[] ca1, Atom[] ca2, FatCatParameters params) throws StructureException{

      AFPChain afpChain = align(ca1,ca2,params,true);
      afpChain.setAlgorithmName(FatCatRigid.algorithmName);
      afpChain.setVersion(VERSION+"");
      return afpChain;
   }

   public AFPChain alignFlexible(Atom[] ca1, Atom[] ca2, FatCatParameters params) throws StructureException{

      AFPChain afpChain = align(ca1,ca2,params,false);
      afpChain.setAlgorithmName(FatCatFlexible.algorithmName);
      afpChain.setVersion(VERSION+"");
      return afpChain;
   }


   protected AFPChain align(Atom[] ca1, Atom[] ca2, FatCatParameters params, boolean doRigid) throws StructureException{

      aligner = new FatCatAligner();

      aligner.align(ca1, ca2, doRigid, params);

      return aligner.getAfpChain();


   }

   public FatCatAligner getFatCatAligner(){
      if ( aligner == null)
         aligner = new FatCatAligner();
      return aligner;
   }

   /** Display the results of an alignment. All input coordinate are not rotated. rotations are done automatically.
    * 
    * @param afpChain
    * @param ca1
    * @param ca2
    * @param hetatms
    * @param nucs
    * @param hetatms2
    * @param nucs2
    * @throws StructureException
    */
//   public StructureAlignmentJmol display(AFPChain afpChain, Atom[] ca1, Atom[] ca2,List<Group> hetatms, List<Group>nucs,List<Group> hetatms2, List<Group>nucs2 )
//   throws StructureException {
//
//
//      FatCatAligner ali = getFatCatAligner();
//      Group[] twistedGroups = ali.getTwistedGroups();
//
//      if ( twistedGroups == null) {	
//         System.out.println("twisting");
//         Atom[] ca3 = StructureTools.cloneCAArray(ca2);
//
//         // this can happen if the alignment got loaded from a flat file
//         twistedGroups = AFPTwister.twistOptimized(afpChain, ca1, ca3);
//         ali.setTwistedGroups(twistedGroups);
//      }
//
//      return DisplayAFP.display(afpChain,twistedGroups, ca1, ca2,hetatms, nucs, hetatms2, nucs2);
//
//   }


}
