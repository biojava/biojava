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

package org.biojava.bio.structure.align.fatcat;

import java.io.IOException;


import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Chain;

import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureImpl;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.xml.AFPChainFlipper;
import org.biojava.bio.structure.align.xml.AFPChainXMLConverter;
import org.biojava.bio.structure.align.xml.AFPChainXMLParser;
import org.biojava.bio.structure.io.PDBFileReader;


import junit.framework.TestCase;

public class FlipAFPChainTest extends TestCase {

	public void testFlipping(){
		try {
			Structure s1 = getStructure("1cdg", "A");		
			Structure s2 = getStructure("1tim","A");

			String name1 = "1cdg.A";
			String name2 = "1tim.A";
			
			Atom[] ca1 = StructureTools.getAtomCAArray(s1);
			Atom[] ca2 = StructureTools.getAtomCAArray(s2);

			FatCat fatCat = new FatCat();
			AFPChain afpChain = fatCat.alignRigid(ca1,ca2);
			afpChain.setName1(name1);
			afpChain.setName2(name2);
			
			

			String xml = AFPChainXMLConverter.toXML(afpChain, ca1, ca2);
			
			AFPChain newC    = AFPChainXMLParser.fromXML(xml, ca1, ca2);			
			AFPChain flipped = AFPChainFlipper.flipChain(newC);

			assertEquals(afpChain.getName1(), flipped.getName2());
			assertEquals(afpChain.getName2(),flipped.getName1());
			assertEquals(afpChain.getCa1Length(),flipped.getCa2Length());
			assertEquals(afpChain.getCa2Length(),flipped.getCa1Length());
			
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
			

		} catch (Exception e){
			e.printStackTrace();
			fail(e.getMessage());
		}
	}



	private Structure getStructure(String pdbId, String chainId) throws IOException, StructureException{
		PDBFileReader pdbpars = new PDBFileReader();
		pdbpars.setPath(AFPChainSerialisationTest.PDB_FILE_PATH);
		pdbpars.setAutoFetch(true);
		Structure structure1 = pdbpars.getStructureById(pdbId);

		Chain c = structure1.getChainByPDB(chainId);

		Structure s = new StructureImpl();
		s.addChain(c);

		return s;

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
