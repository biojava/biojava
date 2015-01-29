/*
 *                    BioJava development code
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
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 * Created on Mar 15, 2010
 * Author: Andreas Prlic 
 *
 */

package demo;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.StructureAlignment;
import org.biojava.bio.structure.align.StructureAlignmentFactory;

import org.biojava.bio.structure.align.fatcat.FatCatRigid;
import org.biojava.bio.structure.align.fatcat.calc.FatCatParameters;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AtomCache;

public class DemoFATCAT
{

   public static void main(String[] args){

      //String name1 = "4hhb.A";
      //String name2 = "4hhb.B";

      String name1 = "1cdg.A";
      String name2 = "1tim.B";



      AtomCache cache = new AtomCache();

      Structure structure1 = null;
      Structure structure2 = null;

      try {

         StructureAlignment algorithm  = StructureAlignmentFactory.getAlgorithm(FatCatRigid.algorithmName);

         structure1 = cache.getStructure(name1);
         structure2 = cache.getStructure(name2);

         Atom[] ca1 = StructureTools.getAtomCAArray(structure1);
         Atom[] ca2 = StructureTools.getAtomCAArray(structure2);

         // the default parameters
         FatCatParameters params = new FatCatParameters();

         AFPChain afpChain = algorithm.align(ca1,ca2,params);            

         afpChain.setName1(name1);
         afpChain.setName2(name2);

         // flexible original results:
         System.out.println(afpChain.toFatcat(ca1,ca2));

         System.out.println(afpChain.toRotMat());
         //System.out.println(afpChain.toCE(ca1, ca2));

         //System.out.println(AFPChainXMLConverter.toXML(afpChain,ca1,ca2));

      } catch (Exception e) {
         e.printStackTrace();
         return;
      }
   }
}
