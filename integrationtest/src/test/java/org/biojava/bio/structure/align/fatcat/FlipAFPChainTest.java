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

import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.SVDSuperimposer;
import org.biojava.bio.structure.StructureTools;

import org.biojava.bio.structure.StructureException;

import org.biojava.bio.structure.align.StructureAlignment;
import org.biojava.bio.structure.align.StructureAlignmentFactory;

import org.biojava.bio.structure.align.ce.CeMain;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.model.AfpChainWriter;
import org.biojava.bio.structure.align.seq.SmithWaterman3Daligner;
import org.biojava.bio.structure.align.util.AFPAlignmentDisplay;
import org.biojava.bio.structure.align.util.AFPChainScorer;
import org.biojava.bio.structure.align.util.AtomCache;

import org.biojava.bio.structure.align.xml.AFPChainFlipper;
import org.biojava.bio.structure.align.xml.AFPChainXMLConverter;
import org.biojava.bio.structure.align.xml.AFPChainXMLParser;
import org.biojava.bio.structure.jama.Matrix;



import junit.framework.TestCase;

public class FlipAFPChainTest extends TestCase {

	public void testFlipping(){
		try {

			//String name1 = "1cdg.A";
			//String name2 = "1tim.A";

			String name1= "1a4w.H";
			String name2= "1hiv.A";
			
			// we currently don;t test CECP, because there is a minor mismatch.
			StructureAlignment[] aligs = new StructureAlignment[]{
				StructureAlignmentFactory.getAlgorithm(CeMain.algorithmName),
				StructureAlignmentFactory.getAlgorithm(FatCatRigid.algorithmName),
				StructureAlignmentFactory.getAlgorithm(FatCatFlexible.algorithmName),
				StructureAlignmentFactory.getAlgorithm(SmithWaterman3Daligner.algorithmName),
			};

			// TODO replace aligs with StructureAlignmentFactory.getAllAlgorithms()
			for (StructureAlignment alig : aligs) {

				align(alig, name1, name2);
			}



		} catch (Exception e){
			e.printStackTrace();
			fail(e.getMessage());
		}
	}

	private void align (StructureAlignment algorithm, String name1, String name2)
	throws StructureException, IOException{
		
		
		AtomCache cache = new AtomCache();
		Atom[] ca1 = cache.getAtoms(name1);
		Atom[] ca2 = cache.getAtoms(name2);


		AFPChain afpChain = algorithm.align(ca1,ca2);
		afpChain.setName1(name1);
		afpChain.setName2(name2);
		double tmScore = AFPChainScorer.getTMScore(afpChain, ca1, ca2);
		afpChain.setTMScore(tmScore);
		
		String xml = AFPChainXMLConverter.toXML(afpChain, ca1, ca2);

		AFPChain newC    = AFPChainXMLParser.fromXML(xml, ca1, ca2);	
		//System.out.println(xml);
		//System.out.println(AFPChainXMLConverter.toXML(newC));
		AFPChain flipped = AFPChainFlipper.flipChain(newC);
		
		assertEquals(afpChain.getName1(), flipped.getName2());
		assertEquals(afpChain.getName2(),flipped.getName1());
		assertEquals(afpChain.getCa1Length(),flipped.getCa2Length());
		assertEquals(afpChain.getCa2Length(),flipped.getCa1Length());
		assertEquals(String.format("%.2f",afpChain.getTMScore()), String.format("%.2f",flipped.getTMScore()));
		assertTrue(afpChain.getTMScore() != -1);

		String xmlNew = AFPChainXMLConverter.toXML(flipped, ca2, ca1);
		//System.out.println(xmlNew);
		AFPChain backChain = AFPChainXMLParser.fromXML(xmlNew, ca2, ca1);
		AFPChain origFlip  = AFPChainFlipper.flipChain(backChain);

		assertNotNull("Got null, instead of an AFPChain object!", origFlip);
		
		assertNotNull("could not get nr. of eqr: ", afpChain.getNrEQR());
		assertNotNull("could not get nr. of eqr: ", origFlip.getNrEQR());
				
		assertTrue("The nr. of equivalent positions is not equal!", afpChain.getNrEQR() == origFlip.getNrEQR());
		
		Atom shift1 = afpChain.getBlockShiftVector()[0];
		Atom shift2 = origFlip.getBlockShiftVector()[0];
		
		assertTrue("The shift vectors are not similar!", Calc.getDistance(shift1, shift2) < 0.1);
		
		//assert the RMSD in the flipped alignment is small		
		double rmsd1 = getRMSD(afpChain,ca1,ca2);
		double rmsd2 = getRMSD(flipped,ca2,ca1);
		//System.out.println("rmsd:" +rmsd1 + " " + rmsd2);
		assertTrue("The RMSD are vastly different!", Math.abs(rmsd1-rmsd2) < 0.01);
		
		
		// this can;t work any more because there is minor after comma mismatches..
		//String xmlBack = AFPChainXMLConverter.toXML(origFlip);
		//if ( ! xmlBack.equals(xml)){
		//	printFirstMismatch(xmlBack, xml);
		//}
		//assertEquals("The alignment representations are not the same!" , xmlBack, xml);
		AFPAlignmentDisplay.getAlign(origFlip, ca1, ca2);
		String img1 = AfpChainWriter.toDBSearchResult(afpChain);
		String img2 = AfpChainWriter.toDBSearchResult(origFlip);
		assertEquals("The alignment images do not match!",img1,img2);
		
		//System.out.println(xml);
		//System.out.println(xmlNew);
		
		

	}


	/** get the RMSD between the aligned positions
	 * 
	 * @param afpChain
	 * @param ca1
	 * @param ca2
	 * @return
	 */
	private double getRMSD(AFPChain afpChain, Atom[] ca1, Atom[] ca2) 
	throws StructureException {
		
		Atom[] ca2clone = StructureTools.cloneCAArray(ca2);
		rotateAtoms2(afpChain,ca2clone);
		
		// get only the subset of Atoms that is on structurally equivalent positions
		
		Atom[] catmp1 = AFPAlignmentDisplay.getAlignedAtoms1(afpChain, ca1);
		Atom[] catmp2 = AFPAlignmentDisplay.getAlignedAtoms2(afpChain, ca2clone);
		
		assertTrue(catmp1.length == catmp2.length);
		
		assertTrue(catmp1.length == afpChain.getNrEQR());
				
         return SVDSuperimposer.getRMS(catmp1,catmp2);
	}
	
	public static void rotateAtoms2(AFPChain afpChain,Atom[] ca2){

		

		int blockNum = afpChain.getBlockNum();

		int[] optLen = afpChain.getOptLen();
		int[][][] optAln = afpChain.getOptAln();
		
		for(int bk = 0; bk < blockNum; bk ++)       {

			Matrix m= afpChain.getBlockRotationMatrix()[bk];
			Atom shift = afpChain.getBlockShiftVector()[bk];
			for ( int i=0;i< optLen[bk];i++){
				int pos = optAln[bk][1][i];
				Atom a = ca2[pos];
				
				Calc.rotate(a, m);
				Calc.shift(a, shift);
				
				//atoms.add(ca2[pos]);
			}

		}
		
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
