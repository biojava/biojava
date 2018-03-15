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
package org.biojava.nbio.structure.test.align;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.StructureAlignment;
import org.biojava.nbio.structure.align.StructureAlignmentFactory;
import org.biojava.nbio.structure.align.ce.CeCPMain;
import org.biojava.nbio.structure.align.ce.CeMain;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.util.AlignmentTools;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.align.xml.AFPChainXMLParser;
import org.junit.Assert;
import org.junit.Test;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

public class AlignmentToolsTest {

	@Test
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

		Assert.assertFalse("CeCPMain should give non-sequential alignments (between blocks).", AlignmentTools.isSequentialAlignment(afpChain, false));
		Assert.assertFalse("CeCPMain should give non-sequential alignments (within blocks).", AlignmentTools.isSequentialAlignment(afpChain, true));

		// linear case
		ce = StructureAlignmentFactory.getAlgorithm(CeMain.algorithmName);
		afpChain = ce.align(ca1,ca2);

		Assert.assertTrue("CeMain should give sequential alignments (within blocks).", AlignmentTools.isSequentialAlignment(afpChain, true));
		Assert.assertTrue("CeMain should give sequential alignments (between blocks).", AlignmentTools.isSequentialAlignment(afpChain, false));

		// now change the block interior a bit

		int[][][] optAln = afpChain.getOptAln();
		int tmp;
		tmp = optAln[0][0][2];
		optAln[0][0][2] = optAln[0][0][1];
		optAln[0][0][1] = tmp;
		tmp = optAln[0][1][2];
		optAln[0][1][2] = optAln[0][1][1];
		optAln[0][1][1] = tmp;

		Assert.assertTrue("Modifying block interior shouldn't effect block sequence.", AlignmentTools.isSequentialAlignment(afpChain, false));
		Assert.assertFalse("Modifying block interior should be not sequential.", AlignmentTools.isSequentialAlignment(afpChain, true));

	}

	@Test
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
		Assert.assertEquals("Wrong order for alignment 1", 3, order);

		// sequential alignment. Should be order 1, but we report this as "no symmetry"
		//TODO Change default return value in getSymmetry
		Map<Integer,Integer> alignment2 = new HashMap<Integer,Integer>();
		for(int i=1;i<10;i++) {
			alignment2.put(i, i+1);
		}

		order = AlignmentTools.getSymmetryOrder(alignment2, identity, maxSymmetry, minimumMetricChange);
		Assert.assertEquals("Wrong order for alignment 2", 1, order);

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
		Assert.assertEquals("Wrong order for alignment 3 with I(x)=x-10", 3, order);

		/* These tests don't work because there are no paths longer than maxSymmetry, so they hit 0 error (NaN metric change)
		//Stringent minimumMetric values cause it to miss the alignment
		order = AlignmentTools.getSymmetryOrder(alignment3, identityMinus10, maxSymmetry, .001f);
		assertEquals("Wrong order for alignment 1 with I(x)=x+1 & minMetric=.01",1,order);

		order = AlignmentTools.getSymmetryOrder(alignment1, identity, maxSymmetry, .001f);
		assertEquals("Wrong order for alignment 1 & minMetric=.01",1,order);
		*/
	}

	@Test
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
		Assert.assertEquals("Wrong forward alignment", sequentialForward, result);

		result = AlignmentTools.guessSequentialAlignment(alignment1, true);
		Assert.assertEquals("Wrong backward alignment", sequentialBackward, result);
	}


	@Test
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
		//long startAlignmentTime = System.currentTimeMillis();
		AFPChain afpChain = cecp.align(ca1, ca2);
		//long alignmentTime = System.currentTimeMillis() - startAlignmentTime;

		final int maxSymmetry = 8;
		final float minimumMetricChange = .4f;

		//long startSymmetryOrderTime = System.currentTimeMillis();
		int order = AlignmentTools.getSymmetryOrder(afpChain, maxSymmetry, minimumMetricChange);
		//long symmetryOrderTime = System.currentTimeMillis() - startSymmetryOrderTime;

		//System.out.println("Len1\tLen2\tAlignT\tOrderT\tOrder");
		//System.out.format("%d\t%d\t%f\t%f\t%d", ca1.length,ca2.length,
		//		alignmentTime/1000.,symmetryOrderTime/1000.,order);
		//System.out.println();

		Assert.assertEquals("Wrong order found for " + name1 + " vs " + name2, trueOrder, order);
	}

	@Test
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
		Assert.assertEquals("Alignment1 incorrectly applied", image1, result1);
	}

	@Test
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
		Assert.assertEquals("Alignment1 incorrectly applied with identity x->x-10", image1, result1);
	}

	@Test
	public void testToConciseAlignmentString() {
		Map<Integer,Integer> test;
		String result,expected;
		int i=0;

		test = new HashMap<Integer, Integer>();
		test.put(1, 2);
		test.put(2, 3);
		test.put(3, 4);
		test.put(7, 8);
		expected = "1>2>3>4 7>8";

		result = AlignmentTools.toConciseAlignmentString(test);
		Assert.assertEquals((i++) + ". Linear strings.", expected, result);


		test = new HashMap<Integer, Integer>();
		test.put(1, 2);
		test.put(2, 3);
		test.put(3, 1);
		test.put(7, 7);
		expected = "1>2>3>1 7>7";

		result = AlignmentTools.toConciseAlignmentString(test);
		Assert.assertEquals((i++) + ". Cycles.", expected, result);

		test = new HashMap<Integer, Integer>();
		test.put(1, 2);
		test.put(2, 3);
		test.put(3, 1);
		test.put(7, 7);
		expected = "1>2>3>1 7>7";

		result = AlignmentTools.toConciseAlignmentString(test);
		Assert.assertEquals((i++) + ". Complex.", expected, result);

		test = new HashMap<Integer, Integer>();
		test.put(1, 2);
		test.put(2, 3);
		test.put(3, 4);
		test.put(4, 5);
		test.put(5, 6);
		test.put(6, 7);
		test.put(7, 3);
		test.put(8, 4);
		test.put(9, 11);
		test.put(11, 10);
		test.put(10, 9);
		expected = "1>2>3>4>5>6>7>3 8>4 9>11>10>9";

		result = AlignmentTools.toConciseAlignmentString(test);
		Assert.assertEquals((i++) + ". Complex.", expected, result);

		//This test highlights a suboptimal case, where more paths are used than necessary.
		test.remove(2);
		//expected = "1>2 8>4>5>6>7>3 9>11>10>9"; //more optimal, but would require depth first search
		expected = "1>2 3>4>5>6>7>3 8>4 9>11>10>9";

		result = AlignmentTools.toConciseAlignmentString(test);
		Assert.assertEquals((i++) + ". Sub-optimal arrangement", expected, result);

		Map <Integer,Double> test2 = new HashMap<Integer, Double>();
		test2.put(1, 12.);
		test2.put(2, 13.);
		test2.put(3, 14.);
		test2.put(4, 15.);
		test2.put(5, 16.);
		test2.put(6, 17.);
		test2.put(7, 13.);
		test2.put(8, 14.);
		test2.put(9, 21.);
		test2.put(11, 20.);
		test2.put(10, 19.);
		expected = "1>2>3>4>5>6>7>3 8>4 9>11>10>9";

		Map <Double, Integer> inverse = new OffsetMap(-10);

		result = AlignmentTools.toConciseAlignmentString(test2,inverse);
		Assert.assertEquals((i++) + ". Inverse.", expected, result);


	}

	/**
	 * Tests that {@link AlignmentTools#updateSuperposition(AFPChain, Atom[], Atom[])} calculates the correct RMSD and TM-score for an AFPChain of 1 block.
	 * TODO: Write a test with 2 blocks
	 */
	@Test
	public void testUpdateSuperposition() throws IOException, StructureException {
		Structure s = StructureTools.getStructure("31BI");
		Atom[] ca1 = StructureTools.getRepresentativeAtomArray(s);
		Atom[] ca2 = StructureTools.getRepresentativeAtomArray(s);
		StringBuilder sb = new StringBuilder();
		BufferedReader br = new BufferedReader(new FileReader("src/test/resources/align/31BI_symm_align.xml"));
		String line = "";
		while ((line = br.readLine()) != null) {
			sb.append(line);
		}
		br.close();
		AFPChain afpChain = AFPChainXMLParser.fromXML(sb.toString(), ca1, ca2);
		afpChain.setTMScore(-1);
		afpChain.setTotalRmsdOpt(-1);
		AlignmentTools.updateSuperposition(afpChain, ca1, ca2);
		Assert.assertEquals("TM-score is wrong", 0.62779, afpChain.getTMScore(), 0.001);
		Assert.assertEquals("RMSD is wrong", 2.50569, afpChain.getTotalRmsdOpt(), 0.001);
	}

	/**
	 * Maps (Double d)->((int)(d+offset))
	 * @author blivens
	 *
	 */
	public static class OffsetMap extends AbstractMap<Double,Integer>
	{
		private int offset;
		public OffsetMap(int offset) {
			this.offset = offset;
		}
		@Override
		public Integer get(Object key) {
			if(key==null) return null;
			return (int)((Double)key+offset);
		}

		/**
		 * Always returns the empty set
		 */
		@Override
		public Set<java.util.Map.Entry<Double,Integer>> entrySet() {
			return Collections.emptySet();
		}

		@Override
		public boolean containsKey(Object key) {
			return true;
		}
	}
	
}
