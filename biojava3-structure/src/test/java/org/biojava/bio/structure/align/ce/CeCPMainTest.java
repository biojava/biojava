/**
 * 
 */
package org.biojava.bio.structure.align.ce;

import static org.junit.Assert.*;

import java.lang.reflect.Method;

import org.biojava.bio.structure.AminoAcidImpl;
import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.AtomImpl;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.ChainImpl;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.ResidueNumber;
import org.biojava.bio.structure.align.ce.CECalculator;
import org.biojava.bio.structure.align.ce.CeCPMain;
import org.biojava.bio.structure.align.model.AFPChain;
import org.junit.*;


/**
 * @author Spencer Bliven
 *
 */
public class CeCPMainTest {

	@Test
	public void testFilterDuplicateAFPs() throws Exception {
		Method filterDuplicateAFPs = CeCPMain.class.getDeclaredMethod(
				"filterDuplicateAFPs", AFPChain.class, CECalculator.class, Atom[].class,Atom[].class);
		filterDuplicateAFPs.setAccessible(true);


		int[][][] dupAlign = new int[1][2][];

		int ca2len = 12;
		dupAlign[0][0] = new int[] { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13 };
		dupAlign[0][1] = new int[] { 3, 5, 6, 7, 8, 9,10,11, 0+ca2len, 1+ca2len, 2+ca2len, 3+ca2len, 4+ca2len, 7+ca2len };

		Chain chain1 = new ChainImpl();
		//Some dummy Atoms. Just check they're unique
		Atom[] ca1 = new Atom[14];
		for(int i=0;i<ca1.length;i++) {
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
		Atom[] ca2 = new Atom[ca2len*2];
		Chain chain2 = new ChainImpl();
		for(int i=0;i<ca2.length;i++) {
			ca2[i] = new AtomImpl();
			ca2[i].setFullName(" CA ");
			ca2[i].setName("CA");
			ca2[i].setCoords(new double[] { i%ca2len, 1, 1 });
			Group aa = new AminoAcidImpl();
			aa.addAtom(ca2[i]);
			aa.setPDBName("GLY");
			aa.setResidueNumber( ResidueNumber.fromString(i+""));
			chain2.addGroup(aa);
		}

		AFPChain afp = new AFPChain();
		afp.setOptAln(dupAlign);
		afp.setOptLength(dupAlign[0][1].length);
		afp.setCa2Length(ca2.length);
		afp.setBlockNum(1);
		afp.setOptLen(new int[] {dupAlign[0][1].length});

		AFPChain newAFP = (AFPChain) filterDuplicateAFPs.invoke(null, afp, new CECalculator(null), ca1, ca2);

		int[][][] align = newAFP.getOptAln();
		int[] blkLen = newAFP.getOptLen();
		// optimal alignment should be
		//  1  2  3  4  5  6  7 | 8  9 10 11 12
		//  5  6  7  8  9 10 11 | 0  1  2  3  4

		assertEquals(2,align.length);
		assertEquals(2,align[0].length);
		assertEquals(7,align[0][0].length);
		assertEquals(7,align[0][1].length);
		assertEquals(2,align[1].length);
		assertEquals(5,align[1][0].length);
		assertEquals(5,align[1][1].length);

		assertEquals(2,blkLen.length);
		assertEquals(7,blkLen[0]);
		assertEquals(5,blkLen[1]);


		for(int i=0;i<align[0][0].length; i++) {
			assertEquals(align[0][0][i],i+1);
			assertEquals(align[0][1][i],(i+5)%12);
		}
	}
}
