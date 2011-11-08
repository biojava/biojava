/**
 * 
 */
package org.biojava.bio.structure.align.ce;

//import static org.junit.Assert.*;

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
import org.biojava.bio.structure.io.PDBParseException;
//import org.junit.*;


/**
 * @author Spencer Bliven
 *
 */
public class CeCPMainTest extends TestCase {

	//@Test
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
		AFPChain newAFP = CeCPMain.filterDuplicateAFPs(afp, new CECalculator(null), ca1, ca2);

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

	//@Test
	public void testFilterDuplicateAFPsMinLen() throws PDBParseException, StructureException {
		int[][][] startAln, filteredAln;
		int ca2len;
		
		// minCPLen == 5
		ca2len = 10;
		startAln = new int[][][] {
				new int[][] {
						new int[] {  0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,},
						new int[] {  0, 3, 5, 6, 7, 8, 9, 0+ca2len, 1+ca2len, 2+ca2len, 3+ca2len,},
				},
		};
		
		// The longest alignment would include the second 0-3
		// However, this leads to a short block
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
		
		int[] filteredLen = new int[] { filteredAln[0][0].length, filteredAln[1][0].length };
		
		Atom[] ca1, ca2;
		AFPChain afpChain;
		ca1 = makeDummyCA(startAln[0][0].length);
		ca2 = makeDummyCA(ca2len);
		ca2 = StructureTools.duplicateCA2(ca2);
		afpChain = makeDummyAFPChain(startAln, ca1, ca2);
		afpChain = CeCPMain.filterDuplicateAFPs(afpChain, new CECalculator(null), ca1, ca2);
		
		assertTrue(Arrays.deepEquals(filteredAln, afpChain.getOptAln()));
		assertTrue(Arrays.equals(filteredLen, afpChain.getOptLen()));

	}


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
}
