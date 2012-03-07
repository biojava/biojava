/**
 * 
 */
package org.biojava.bio.structure.align.ce;

//import static org.junit.Assert.*;

import java.lang.reflect.InvocationTargetException;
import java.util.Arrays;

import junit.framework.TestCase;

import org.biojava.bio.structure.AminoAcidImpl;
import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.AtomImpl;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.ChainImpl;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.ResidueNumber;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureTools;

import org.biojava.bio.structure.align.ce.CECalculator;
import org.biojava.bio.structure.align.ce.CeCPMain;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.io.PDBParseException;
import org.junit.*;


/**
 * @author Spencer Bliven
 *
 */
public class CeCPMainTest extends TestCase {

	@Test
	public void testFilterDuplicateAFPs() throws Exception {
		int[][][] dupAlign = new int[1][2][];

		int ca2len = 12;
		dupAlign[0][0] = new int[] { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13 };
		dupAlign[0][1] = new int[] { 3, 5, 6, 7, 8, 9,10,11, 0+ca2len, 1+ca2len, 2+ca2len, 3+ca2len, 4+ca2len, 7+ca2len };
		Atom[] ca1,ca2;
		
		ca1 = makeDummyCA(dupAlign[0][0].length);
		ca2 = makeDummyCA(ca2len);
		ca2 = StructureTools.duplicateCA2(ca2);
		AFPChain afp = makeDummyAFPChain(dupAlign, ca1, ca2);

//		AFPChain newAFP = (AFPChain) filterDuplicateAFPs.invoke(null, afp, new CECalculator(null), ca1, ca2);
		AFPChain newAFP = CeCPMain.filterDuplicateAFPs(afp, new CECalculator(null), ca1, ca2, 0);

		int[][][] align = newAFP.getOptAln();
		int[] blkLen = newAFP.getOptLen();
		// optimal alignment should be
		//  1  2  3  4  5  6  7 | 8  9 10 11 12
		//  5  6  7  8  9 10 11 | 0  1  2  3  4
		int[][][] expected = new int[][][] {
				new int[][] {
						new int[] { 1,  2,  3,  4,  5,  6,  7, },
						new int[] { 5,  6,  7,  8,  9, 10, 11, },
				},
				new int[][] {
						new int[] { 8,  9, 10, 11, 12, },
						new int[] { 0,  1,  2,  3,  4, },
				},
		};
		
		int[] expectedLen = new int[] { expected[0][0].length, expected[1][0].length };
		
		assertTrue(Arrays.deepEquals(expected, align));
		assertTrue(Arrays.equals(expectedLen, blkLen));
		
	}

	@Test
	public void testFilterDuplicateAFPsMinLenCTerm() throws PDBParseException, StructureException {
		int[][][] startAln, filteredAln;
		int[] filteredLen;
		int ca2len;
		
		ca2len = 10;
		startAln = new int[][][] {
				new int[][] {
						new int[] {  0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,},
						new int[] {  0, 3, 5, 6, 7, 8, 9, 0+ca2len, 1+ca2len, 2+ca2len, 3+ca2len,},
				},
		};
		
		Atom[] ca1, ca2;
		AFPChain afpChain,result;
		ca1 = makeDummyCA(startAln[0][0].length);
		ca2 = makeDummyCA(ca2len);
		ca2 = StructureTools.duplicateCA2(ca2);
		afpChain = makeDummyAFPChain(startAln, ca1, ca2);

		
		
		// Best block with minCPlen 0-3
		filteredAln = new int[][][] {
				new int[][] {
						new int[] {  1, 2, 3, 4, 5, 6,},
						new int[] {  3, 5, 6, 7, 8, 9,},
				},
				new int[][] {
						new int[] { 7, 8, 9,},
						new int[] { 0, 1, 2,},
				},
		};
		filteredLen = new int[] { filteredAln[0][0].length, filteredAln[1][0].length };
		
		
		
		int minCPlength;
		for(minCPlength=0;minCPlength<4;minCPlength++) {
			result = CeCPMain.filterDuplicateAFPs(afpChain, new CECalculator(null), ca1, ca2,minCPlength);
		
			assertTrue("Wrong optAln for minCPlength="+minCPlength,Arrays.deepEquals(filteredAln, result.getOptAln()));
			assertTrue("Wrong optLen for minCPlength="+minCPlength,Arrays.equals(filteredLen, result.getOptLen()));
	}

		// For minCPlength=4, filtering changes
		minCPlength=4;
		filteredAln = new int[][][] {
				new int[][] {
						new int[] {  2, 3, 4, 5, 6,},
						new int[] {  5, 6, 7, 8, 9,},
				},
				new int[][] {
						new int[] { 7, 8, 9,10,},
						new int[] { 0, 1, 2, 3,},
				},
		};
		filteredLen = new int[] { filteredAln[0][0].length, filteredAln[1][0].length };

		result = CeCPMain.filterDuplicateAFPs(afpChain, new CECalculator(null), ca1, ca2,minCPlength);

		assertTrue("Wrong optAln for minCPlength="+minCPlength,Arrays.deepEquals(filteredAln, result.getOptAln()));
		assertTrue("Wrong optLen for minCPlength="+minCPlength,Arrays.equals(filteredLen, result.getOptLen()));

		// For minCPlength=5, filtering changes
		minCPlength=5;
		filteredAln = new int[][][] {
				new int[][] {
						new int[] {  0, 1, 2, 3, 4, 5, 6,},
						new int[] {  0, 3, 5, 6, 7, 8, 9,},
				},
		};
		filteredLen = new int[] { filteredAln[0][0].length };

		result = CeCPMain.filterDuplicateAFPs(afpChain, new CECalculator(null), ca1, ca2,minCPlength);

		assertTrue("Wrong optAln for minCPlength="+minCPlength,Arrays.deepEquals(filteredAln, result.getOptAln()));
		assertTrue("Wrong optLen for minCPlength="+minCPlength,Arrays.equals(filteredLen, result.getOptLen()));

		minCPlength=7;
		result = CeCPMain.filterDuplicateAFPs(afpChain, new CECalculator(null), ca1, ca2,minCPlength);

		assertTrue("Wrong optAln for minCPlength="+minCPlength,Arrays.deepEquals(filteredAln, result.getOptAln()));
		assertTrue("Wrong optLen for minCPlength="+minCPlength,Arrays.equals(filteredLen, result.getOptLen()));

		// Eventually, no alignment!
		minCPlength=8;
		filteredAln = new int[0][][];
		filteredLen = new int[0];


		result = CeCPMain.filterDuplicateAFPs(afpChain, new CECalculator(null), ca1, ca2,minCPlength);

		assertTrue("Wrong optAln for minCPlength="+minCPlength,Arrays.deepEquals(filteredAln, result.getOptAln()));
		assertTrue("Wrong optLen for minCPlength="+minCPlength,Arrays.equals(filteredLen, result.getOptLen()));

	}

	@Test
	public void testFilterDuplicateAFPsMinLenNTerm() throws PDBParseException, StructureException {
		int[][][] startAln, filteredAln;
		int ca2len;

		// minCPLen == 5
		ca2len = 10;
		startAln = new int[][][] {
				new int[][] {
						new int[] {  0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,},
						new int[] {  6, 8, 9,10,11,12,13,14,16,17,19,},
				},
		};

		// The longest alignment would include the second 0-3
		// However, this leads to a short block
		filteredAln = new int[][][] {
				new int[][] {
						new int[] {  1, 2,},
						new int[] {  8, 9,},
				},
				new int[][] {
						new int[] {  3, 4, 5, 6, 7, 8, 9,},
						new int[] {  0, 1, 2, 3, 4, 6, 7,},
				},
		};

		int[] filteredLen = new int[] { filteredAln[0][0].length, filteredAln[1][0].length };

		Atom[] ca1, ca2;
		AFPChain afpChain,result;
		ca1 = makeDummyCA(startAln[0][0].length);
		ca2 = makeDummyCA(ca2len);
		ca2 = StructureTools.duplicateCA2(ca2);
		afpChain = makeDummyAFPChain(startAln, ca1, ca2);

		int minCPlength;
		for(minCPlength=0;minCPlength<3;minCPlength++) {
			result = CeCPMain.filterDuplicateAFPs(afpChain, new CECalculator(null), ca1, ca2,minCPlength);

			assertTrue("Wrong optAln for minCPlength="+minCPlength,Arrays.deepEquals(filteredAln, result.getOptAln()));
			assertTrue("Wrong optLen for minCPlength="+minCPlength,Arrays.equals(filteredLen, result.getOptLen()));
		}

		// For minCPlength=3, filtering changes
		minCPlength=3;
		filteredAln = new int[][][] {
				new int[][] {
						new int[] {  0, 1, 2,},
						new int[] {  6, 8, 9,},
				},
				new int[][] {
						new int[] {  3, 4, 5, 6, 7,},
						new int[] {  0, 1, 2, 3, 4,},
				},
		};
		filteredLen = new int[] { filteredAln[0][0].length, filteredAln[1][0].length };

		result = CeCPMain.filterDuplicateAFPs(afpChain, new CECalculator(null), ca1, ca2,minCPlength);

		assertTrue("Wrong optAln for minCPlength="+minCPlength,Arrays.deepEquals(filteredAln, result.getOptAln()));
		assertTrue("Wrong optLen for minCPlength="+minCPlength,Arrays.equals(filteredLen, result.getOptLen()));

		// For minCPlength=4, filtering changes
		minCPlength=5;
		filteredAln = new int[][][] {
				new int[][] {
						new int[] {  3, 4, 5, 6, 7, 8, 9,10,},
						new int[] {  0, 1, 2, 3, 4, 6, 7, 9,},
				},
		};
		filteredLen = new int[] { filteredAln[0][0].length };

		result = CeCPMain.filterDuplicateAFPs(afpChain, new CECalculator(null), ca1, ca2,minCPlength);

		assertTrue("Wrong optAln for minCPlength="+minCPlength,Arrays.deepEquals(filteredAln, result.getOptAln()));
		assertTrue("Wrong optLen for minCPlength="+minCPlength,Arrays.equals(filteredLen, result.getOptLen()));

		minCPlength=8;
		result = CeCPMain.filterDuplicateAFPs(afpChain, new CECalculator(null), ca1, ca2,minCPlength);

		assertTrue("Wrong optAln for minCPlength="+minCPlength,Arrays.deepEquals(filteredAln, result.getOptAln()));
		assertTrue("Wrong optLen for minCPlength="+minCPlength,Arrays.equals(filteredLen, result.getOptLen()));

		// Eventually, no alignment!
		minCPlength=9;
		filteredAln = new int[0][][];
		filteredLen = new int[0];


		result = CeCPMain.filterDuplicateAFPs(afpChain, new CECalculator(null), ca1, ca2,minCPlength);

		assertTrue("Wrong optAln for minCPlength="+minCPlength,Arrays.deepEquals(filteredAln, result.getOptAln()));
		assertTrue("Wrong optLen for minCPlength="+minCPlength,Arrays.equals(filteredLen, result.getOptLen()));

	}


	/**
	 * Creates a minimal AFPChain from the specified alignment and proteins
	 * @param dupAlign
	 * @param ca1
	 * @param ca2
	 * @return
	 */
	private AFPChain makeDummyAFPChain(int[][][] dupAlign, Atom[] ca1,Atom[] ca2) {
		AFPChain afp = new AFPChain();
		afp.setOptAln(dupAlign);
		afp.setOptLength(dupAlign[0][1].length);
		afp.setCa1Length(ca1.length);
		afp.setCa2Length(ca2.length);
		afp.setBlockNum(1);
		afp.setOptLen(new int[] {dupAlign[0][1].length});
		return afp;
	}
	
	
	@Test
	public void testCalculateMinCP() throws SecurityException, NoSuchMethodException, IllegalArgumentException, IllegalAccessException, InvocationTargetException {
		int[] block;
		int ca2len;

		block = new int[] { 4,5,6,8,11,12,14,15, };
		ca2len = 10;

		int minCPlength;
		CeCPMain.CPRange cpRange;

		minCPlength = 0;
		cpRange = CeCPMain.calculateMinCP(block, block.length, ca2len, minCPlength);
		assertEquals("Wrong minCPnterm for minCPlen="+minCPlength,11,cpRange.n);
		assertEquals("Wrong minCPcterm for minCPlen="+minCPlength,8,cpRange.c);

		minCPlength = 1;
		cpRange = CeCPMain.calculateMinCP(block, block.length, ca2len, minCPlength);
		assertEquals("Wrong minCPnterm for minCPlen="+minCPlength,8,cpRange.n);
		assertEquals("Wrong minCPcterm for minCPlen="+minCPlength,11,cpRange.c);

		minCPlength = 2;
		cpRange = CeCPMain.calculateMinCP(block, block.length, ca2len, minCPlength);
		assertEquals("Wrong minCPnterm for minCPlen="+minCPlength,6,cpRange.n);
		assertEquals("Wrong minCPcterm for minCPlen="+minCPlength,12,cpRange.c);

		minCPlength = 3;
		cpRange = CeCPMain.calculateMinCP(block, block.length, ca2len, minCPlength);
		assertEquals("Wrong minCPnterm for minCPlen="+minCPlength,5,cpRange.n);
		assertEquals("Wrong minCPcterm for minCPlen="+minCPlength,14,cpRange.c);

		minCPlength = 4;
		cpRange = CeCPMain.calculateMinCP(block, block.length, ca2len, minCPlength);
		assertEquals("Wrong minCPnterm for minCPlen="+minCPlength,4,cpRange.n);
		assertEquals("Wrong minCPcterm for minCPlen="+minCPlength,15,cpRange.c);

		minCPlength = 5;
		cpRange = CeCPMain.calculateMinCP(block, block.length, ca2len, minCPlength);
		assertEquals("Wrong minCPnterm for minCPlen="+minCPlength,-1,cpRange.n);
		assertEquals("Wrong minCPcterm for minCPlen="+minCPlength,20,cpRange.c);

		block = new int[] {0,9,10,19};
		ca2len = 10;

		minCPlength = 0;
		cpRange = CeCPMain.calculateMinCP(block, block.length, ca2len, minCPlength);
		assertEquals("Wrong minCPnterm for minCPlen="+minCPlength,10,cpRange.n);
		assertEquals("Wrong minCPcterm for minCPlen="+minCPlength,9,cpRange.c);

		minCPlength = 1;
		cpRange = CeCPMain.calculateMinCP(block, block.length, ca2len, minCPlength);
		assertEquals("Wrong minCPnterm for minCPlen="+minCPlength,9,cpRange.n);
		assertEquals("Wrong minCPcterm for minCPlen="+minCPlength,10,cpRange.c);

		minCPlength = 2;
		cpRange = CeCPMain.calculateMinCP(block, block.length, ca2len, minCPlength);
		assertEquals("Wrong minCPnterm for minCPlen="+minCPlength,0,cpRange.n);
		assertEquals("Wrong minCPcterm for minCPlen="+minCPlength,19,cpRange.c);

		minCPlength = 3;
		cpRange = CeCPMain.calculateMinCP(block, block.length, ca2len, minCPlength);
		assertEquals("Wrong minCPnterm for minCPlen="+minCPlength,-1,cpRange.n);
		assertEquals("Wrong minCPcterm for minCPlen="+minCPlength,20,cpRange.c);

		minCPlength = 4;
		cpRange = CeCPMain.calculateMinCP(block, block.length, ca2len, minCPlength);
		assertEquals("Wrong minCPnterm for minCPlen="+minCPlength,-1,cpRange.n);
		assertEquals("Wrong minCPcterm for minCPlen="+minCPlength,20,cpRange.c);

	}

	/**
	 * Makes dummy CA atoms at 1A intervals
	 * 
	 * @param len
	 * @return
	 * @throws PDBParseException
	 */
	private Atom[] makeDummyCA(int len) throws PDBParseException {
		Atom[] ca1;
		Chain chain1 = new ChainImpl();
		//Some dummy Atoms. Just check they're unique
		ca1 = new Atom[len];
		for(int i=0;i<len;i++) {
			ca1[i] = new AtomImpl();
			ca1[i].setFullName(" CA ");
			ca1[i].setName("CA");
			ca1[i].setCoords(new double[] { i, 0, 0 });
			Group aa = new AminoAcidImpl();
			aa.setPDBName("GLY");
			aa.setResidueNumber( ResidueNumber.fromString(i+""));
			aa.addAtom(ca1[i]);
			chain1.addGroup(aa);
		}
		return ca1;
	}
	
	public void testCECP1(){

		String name1 = "PDP:3A2KAc";
		String name2 = "d1wy5a2";

		try {
			CeCPMain algorithm = new CeCPMain();

			AtomCache cache = new AtomCache();

			Atom[] ca1 = cache.getAtoms(name1);
			Atom[] ca2 = cache.getAtoms(name2);

			AFPChain afpChain = algorithm.align(ca1, ca2);
			CECalculator calculator = algorithm.getCECalculator();

			//               System.out.println(calculator.get);
			//StructureAlignmentJmol jmol =
			//StructureAlignmentDisplay.display(afpChain, ca1, ca2);
			if ( ! (afpChain.getBlockNum() == 1)){
				System.out.println(calculator.getLcmp());
				System.out.println(afpChain.toFatcat(ca1, ca2));
			}
			assertEquals(1,afpChain.getBlockNum());




		} catch (Exception e){
			e.printStackTrace();
		}

	}

}
