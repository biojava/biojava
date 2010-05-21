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
 * Created on Sep 9, 2009
 * Author: Andreas Prlic 
 *
 */

package org.biojava.bio.structure.align;

import java.io.IOException;


import org.biojava.bio.structure.Atom;

import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.TmpAtomCache;
import org.biojava.bio.structure.align.ce.CeCPMain;
import org.biojava.bio.structure.align.ce.CeMain;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.align.xml.AFPChainFlipper;
import org.biojava.bio.structure.align.xml.AFPChainXMLConverter;
import org.biojava.bio.structure.align.xml.AFPChainXMLParser;



import junit.framework.TestCase;

public class FlipAFPChainTest extends TestCase {

   public void testFlippingMultiBlock(){
      
      try {

         String name1 = "1tim.A";
         String name2 = "1cdg.A";
         
         flip(name1,name2, CeCPMain.algorithmName);
      } catch (Exception e){
         e.printStackTrace();
         fail(e.getMessage());
      }

   }

   public void testFlipping(){
      
      try {

         String name1 = "1cdg.A";
         String name2 = "1tim.A";
         flip(name1,name2, CeMain.algorithmName);
      } catch (Exception e){
         e.printStackTrace();
         fail(e.getMessage());
      }
   }

   private void flip(String name1, String name2, String algorithmName) throws StructureException, IOException{
      
      AtomCache cache = TmpAtomCache.cache;
      
      Structure s1 = cache.getStructure(name1);
      Structure s2 = cache.getStructure(name2);

      Atom[] ca1 = StructureTools.getAtomCAArray(s1);
      Atom[] ca2 = StructureTools.getAtomCAArray(s2);

      StructureAlignment algorithm = StructureAlignmentFactory.getAlgorithm(algorithmName );
      AFPChain afpChain = algorithm.align(ca1,ca2);
      afpChain.setName1(name1);
      afpChain.setName2(name2);



      String xml = AFPChainXMLConverter.toXML(afpChain, ca1, ca2);

      AFPChain newC    = AFPChainXMLParser.fromXML(xml, ca1, ca2);			
      AFPChain flipped = AFPChainFlipper.flipChain(newC);

      assertEquals(afpChain.getName1(), flipped.getName2());
      assertEquals(afpChain.getName2(),flipped.getName1());
      assertEquals(afpChain.getCa1Length(),flipped.getCa2Length());
      assertEquals(afpChain.getCa2Length(),flipped.getCa1Length());
      assertEquals(afpChain.getAlgorithmName(),flipped.getAlgorithmName());
      assertEquals(afpChain.getVersion(), flipped.getVersion());

      //System.out.println(AFPChainXMLConverter.toXML(flipped));

      //AFPChainXMLParser.rebuildAFPChain(flipped, ca2, ca1);

      //FatCat newCat = new FatCat();

      //Group[] twistedGroups = AFPTwister.twistOptimized(flipped,ca2,ca1);

      // FatCatAligner aligner =  newCat.getFatCatAligner();
      //aligner.setTwistedGroups(twistedGroups);			
      //newCat.display(flipped, ca2, ca1,  new ArrayList<Group>(),new ArrayList<Group>(),new ArrayList<Group>(),new ArrayList<Group>());

      String xmlNew = AFPChainXMLConverter.toXML(flipped, ca2, ca1);

      AFPChain backChain = AFPChainXMLParser.fromXML(xmlNew, ca2, ca1);
      
      AFPChain origFlip  = AFPChainFlipper.flipChain(backChain);
      //AFPChainXMLParser.rebuildAFPChain(origFlip, ca1, ca2);

      String xmlBack = AFPChainXMLConverter.toXML(origFlip);
      if ( ! xmlBack.equals(xml)){
         printFirstMismatch(xmlBack, xml);
      }
      assertEquals(xmlBack, xml);



   }






   static final String newline = System.getProperty("line.separator");
   public void printFirstMismatch(String s1, String s2){
      String[] spl1 = s1.split(newline);
      String[] spl2 = s2.split(newline);

      for (int i = 0 ; i < spl1.length ; i++){

         String line1 = spl1[i];

         if ( i >= spl2.length){
            System.err.println("s2 does not contain line " + (i+1));
            return;
         }
         String line2 = spl2[i];

         if ( line1.equals(line2)){
            continue;
         }

         System.err.println("mismatch in line: " + (i+1));

         for ( int j = 0 ; j < line1.length();j++){
            char c1 = line1.charAt(j);

            if ( j >= line2.length()){
               System.err.println("s2 is shorter than s1. length s1:" + line1.length() + " length2:" + line2.length() );
               return;
            }

            char c2 = line2.charAt(j);
            if ( c1 != c2){

               System.err.println("line1: " + line1.substring(0,j+1));
               System.err.println("line2: " + line2.substring(0,j+1));

               System.err.println("mismatch at position " + (j+1) + " c1: "+ c1 + " " + c2);

               return;
            }
         }


      }

   }
}
