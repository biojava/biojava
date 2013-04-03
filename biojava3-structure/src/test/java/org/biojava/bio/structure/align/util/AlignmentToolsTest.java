/*
 *                  BioJava development code
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
 * Created on Jun 8, 2007
 *
 */
package org.biojava.bio.structure.align.util;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.StructureAlignment;
import org.biojava.bio.structure.align.StructureAlignmentFactory;
import org.biojava.bio.structure.align.ce.CeCPMain;
import org.biojava.bio.structure.align.ce.CeMain;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AtomCache;

import junit.framework.TestCase;

public class AlignmentToolsTest extends TestCase {
	
	public void testIsSequential() throws StructureException, IOException {
		AtomCache cache = new AtomCache();
		
		String name1, name2;
		Atom[] ca1, ca2;
		AFPChain afpChain;
		StructureAlignment ce;
		
		
		// CP case
		name1="1QDM.A"; // swaposin
		name2="1NKL"; // saposin
		
		ca1=cache.getAtoms(name1);
		ca2=cache.getAtoms(name2);
		
		ce = StructureAlignmentFactory.getAlgorithm(CeCPMain.algorithmName);
		afpChain = ce.align(ca1,ca2);
		
		assertFalse("CeCPMain should give non-sequential alignments (between blocks).",AlignmentTools.isSequentialAlignment(afpChain,false));
		assertFalse("CeCPMain should give non-sequential alignments (within blocks).",AlignmentTools.isSequentialAlignment(afpChain,true));

		// linear case		
		ce = StructureAlignmentFactory.getAlgorithm(CeMain.algorithmName);
		afpChain = ce.align(ca1,ca2);
		
		assertTrue("CeMain should give sequential alignments (within blocks).",AlignmentTools.isSequentialAlignment(afpChain,true));
		assertTrue("CeMain should give sequential alignments (between blocks).",AlignmentTools.isSequentialAlignment(afpChain,false));

		// now change the block interior a bit
		
		int[][][] optAln = afpChain.getOptAln();
		int tmp;
		tmp = optAln[0][0][2];
		optAln[0][0][2] = optAln[0][0][1];
		optAln[0][0][1] = tmp;
		tmp = optAln[0][1][2];
		optAln[0][1][2] = optAln[0][1][1];
		optAln[0][1][1] = tmp;
		
		assertTrue("Modifying block interior shouldn't effect block sequence.",AlignmentTools.isSequentialAlignment(afpChain,false));
		assertFalse("Modifying block interior should be not sequential.",AlignmentTools.isSequentialAlignment(afpChain,true));

	}
	
	public void testGetSymmetryOrderForMaps() {
		int order;
		final int maxSymmetry = 8;
		final float minimumMetricChange = .5f;// be liberal, since we have small alignments
		
		// noisy C3 alignment
		Map<Integer,Integer> alignment1 = new HashMap<Integer,Integer>();
		alignment1.put(1, 5);
		alignment1.put(2, 6);
		alignment1.put(4, 7);
		alignment1.put(6, 9);
		alignment1.put(7, 11);
		alignment1.put(9, 2);
		alignment1.put(10, 3);
		alignment1.put(11, 4);

		Map<Integer,Integer> identity = new AlignmentTools.IdentityMap<Integer>();
		
		order = AlignmentTools.getSymmetryOrder(alignment1, identity, maxSymmetry, minimumMetricChange);
		assertEquals("Wrong order for alignment 1",3,order);
		
		// sequential alignment. Should be order 1, but we report this as "no symmetry"
		//TODO Change default return value in getSymmetry
		Map<Integer,Integer> alignment2 = new HashMap<Integer,Integer>();
		for(int i=1;i<10;i++) {
			alignment2.put(i, i+1);
		}
		
		order = AlignmentTools.getSymmetryOrder(alignment2, identity, maxSymmetry, minimumMetricChange);
		assertEquals("Wrong order for alignment 2",1,order);
		
		// now try to get symmetry order with an imperfect identity
		Map<Integer,Integer> alignment3 = new HashMap<Integer,Integer>();
		alignment3.put(1, 15);
		alignment3.put(2, 16);
		alignment3.put(4, 17);
		alignment3.put(6, 19);
		alignment3.put(7, 21);
		alignment3.put(9, 12);
		alignment3.put(10, 13);
		alignment3.put(11, 14);
		
		Map<Integer,Integer> identityMinus10 = new HashMap<Integer,Integer>();
		for(int i=1;i<=11;i++) {
			identityMinus10.put(i+10, i);
		}
		
		
		order = AlignmentTools.getSymmetryOrder(alignment3, identityMinus10, maxSymmetry, minimumMetricChange);
		assertEquals("Wrong order for alignment 3 with I(x)=x-10",3,order);
		
		/* These tests don't work because there are no paths longer than maxSymmetry, so they hit 0 error (NaN metric change)
		//Stringent minimumMetric values cause it to miss the alignment
		order = AlignmentTools.getSymmetryOrder(alignment3, identityMinus10, maxSymmetry, .001f);
		assertEquals("Wrong order for alignment 1 with I(x)=x+1 & minMetric=.01",1,order);

		order = AlignmentTools.getSymmetryOrder(alignment1, identity, maxSymmetry, .001f);
		assertEquals("Wrong order for alignment 1 & minMetric=.01",1,order);
		*/
	}
	
	public void testGuessSequentialAlignment() {
		// noisy C3 alignment
		Map<Integer,Integer> alignment1 = new HashMap<Integer,Integer>();
		alignment1.put(1, 5);
		alignment1.put(2, 6);
		alignment1.put(4, 7);
		alignment1.put(6, 9);
		alignment1.put(7, 11);
		alignment1.put(9, 2);
		alignment1.put(10, 3);
		alignment1.put(11, 4);
		
		// Sequential version of the alignment
		Map<Integer,Integer> sequentialForward = new HashMap<Integer,Integer>();
		sequentialForward.put(1, 2);
		sequentialForward.put(2, 3);
		sequentialForward.put(4, 4);
		sequentialForward.put(6, 5);
		sequentialForward.put(7, 6);
		sequentialForward.put(9, 7);
		sequentialForward.put(10, 9);
		sequentialForward.put(11, 11);
		
		// inverse of sequentialForward
		Map<Integer,Integer> sequentialBackward = new HashMap<Integer,Integer>();
		sequentialBackward.put(2, 1);
		sequentialBackward.put(3, 2);
		sequentialBackward.put(4, 4);
		sequentialBackward.put(5, 6);
		sequentialBackward.put(6, 7);
		sequentialBackward.put(7, 9);
		sequentialBackward.put(9, 10);
		sequentialBackward.put(11, 11);
		
	
		Map<Integer,Integer> result;
		
		result = AlignmentTools.guessSequentialAlignment(alignment1, false);
		assertEquals("Wrong forward alignment",sequentialForward,result);
	
		result = AlignmentTools.guessSequentialAlignment(alignment1, true);
		assertEquals("Wrong backward alignment",sequentialBackward,result);
	}
	
	 
	public void testGetSymmetryOrderWithCECP() throws IOException, StructureException {
		
		String name1,name2;
		int trueOrder;
		
		// Two highly-symmetric circularly permuted proteins (swaposin)
		name1 = "1QDM.A";
		name2 = "1NKL";
		trueOrder = 2;
		
//		// Non-symmetric
//		name1 = "1NLS.A";
//		name2 = "1RIN.A";
//		trueOrder = 1;
//		
//		// non-symmetric
//		name1 = "1ATG.A";
//		name2 = "2B4L.A";
//		trueOrder = 1;
//		
//		name1 = "1TIM.A";
//		name2 = "1CDG";
//		trueOrder = 1;
//		
//		name1 = "1a22.A";
//		name2 = "2ffx.J";
		
		
		AtomCache cache = new AtomCache();
		Atom[] ca1 = cache.getAtoms(name1);
		Atom[] ca2 = cache.getAtoms(name2);
		
		StructureAlignment cecp = StructureAlignmentFactory.getAlgorithm(CeCPMain.algorithmName);
		long startAlignmentTime = System.currentTimeMillis();
		AFPChain afpChain = cecp.align(ca1, ca2);
		long alignmentTime = System.currentTimeMillis() - startAlignmentTime;
		
		final int maxSymmetry = 8;
		final float minimumMetricChange = .4f;
		
		long startSymmetryOrderTime = System.currentTimeMillis();
		int order = AlignmentTools.getSymmetryOrder(afpChain, maxSymmetry, minimumMetricChange);
		long symmetryOrderTime = System.currentTimeMillis() - startSymmetryOrderTime;
		
		//System.out.println("Len1\tLen2\tAlignT\tOrderT\tOrder");
		//System.out.format("%d\t%d\t%f\t%f\t%d", ca1.length,ca2.length,
		//		alignmentTime/1000.,symmetryOrderTime/1000.,order);
		//System.out.println();

		assertEquals("Wrong order found for "+name1+" vs "+name2,trueOrder,order);
	}
	
	public void testApplyAlignment() {
		// noisy C3 alignment
		Map<Integer,Integer> alignment1 = new HashMap<Integer,Integer>();
		alignment1.put(1, 5);
		alignment1.put(2, 6);
		alignment1.put(4, 7);
		alignment1.put(6, 9);
		alignment1.put(7, 11);
		alignment1.put(9, 2);
		alignment1.put(10, 3);
		alignment1.put(11, 4);

		Map<Integer,Integer> image1 = new HashMap<Integer,Integer>();
		image1.put(1, null);
		image1.put(2, 9);
		image1.put(4, 11);
		image1.put(6, 2);
		image1.put(7, 4);
		image1.put(9, 6);
		image1.put(10, null);
		image1.put(11, 7);
		//image1.put(5, null);
		//image1.put(3, null);
		//TODO handle nulls consistently. Either include all of them, or none.
		
		Map<Integer,Integer> result1 = AlignmentTools.applyAlignment(alignment1,2);
		assertEquals("Alignment1 incorrectly applied",image1,result1);
	}
	
	public void testApplyAlignmentNonIdentical() {
		// noisy C3 alignment
		Map<Integer,Integer> alignment1 = new HashMap<Integer,Integer>();
		alignment1.put(1, 15);
		alignment1.put(2, 16);
		alignment1.put(4, 17);
		alignment1.put(6, 19);
		alignment1.put(7, 21);
		alignment1.put(9, 12);
		alignment1.put(10, 13);
		alignment1.put(11, 14);

		Map<Integer,Integer> image1 = new HashMap<Integer,Integer>();
		image1.put(1, null);
		image1.put(2, 19);
		image1.put(4, 21);
		image1.put(6, 12);
		image1.put(7, 14);
		image1.put(9, 16);
		image1.put(10, null);
		image1.put(11, 17);
		//image1.put(5, null);
		//image1.put(3, null);
		
		Map<Integer,Integer> identity1 = new HashMap<Integer, Integer>();
		for(int i=1;i<12;i++) {
			identity1.put(i+10,i);
		}
		Map<Integer,Integer> result1 = AlignmentTools.applyAlignment(alignment1,identity1,2);
		assertEquals("Alignment1 incorrectly applied with identity x->x-10",image1,result1);
	}
}
