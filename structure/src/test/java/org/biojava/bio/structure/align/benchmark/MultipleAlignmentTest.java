/**
 * 
 */
package org.biojava.bio.structure.align.benchmark;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.align.util.AtomCache;

import junit.framework.TestCase;

/**
 * @author Spencer Bliven
 *
 */
public class MultipleAlignmentTest extends TestCase {

	public void testGetAlignmentMatrix() {
		/* Test on the following alignment:
		 #d1hcy_2-d1lnlb1
		 HIS.194._._	HIS.41._.B
		 HIS.198._._	HIS.61._.B
		 HIS.224._._	HIS.70._.B
		 HIS.344._._	HIS.181._.B
		 */
		String[] pdbIDs = new String[] {"1hcy.B", "1lnl.B"};

		PDBResidue[][] residues = new PDBResidue[][] {
				new PDBResidue[] {
						new PDBResidue("194", "B", "HIS"), 
						new PDBResidue("198", "B", "HIS"), 
						new PDBResidue("224", "B", "HIS"),
						new PDBResidue("344", "B", "HIS")
				},
				new PDBResidue[] {
						new PDBResidue("41", "B", "HIS"), 
						new PDBResidue("61", "B", "HIS"), 
						new PDBResidue("70", "B", "HIS"),
						new PDBResidue("181", "B", "HIS")
				},
		};
		//Expected output
		int[][] resIndices = new int[][] {
				new int[] { 193, 197, 223, 343, },
				new int[] { 43, 63, 72, 183 },
		};


		try {
			MultipleAlignment align = new MultipleAlignment(pdbIDs, residues);

			AtomCache cache = new AtomCache("/tmp/", true);
			List<Atom[]> structures = new ArrayList<Atom[]>();
			for(String pdb :pdbIDs) {
				structures.add(cache.getAtoms(pdb));
			}

			int[][] alignMat = align.getAlignmentMatrix(pdbIDs, structures);
			
			assertTrue(Arrays.deepEquals(alignMat,resIndices));
						
//			if(! Arrays.deepEquals(alignMat,resIndices)) {
//				String msg = "Error in getAlignmentMatrix(). Expected:\n";
//				msg += Arrays.deepToString(resIndices);
//				msg += "\nFound:;";
//				msg += Arrays.deepToString(alignMat);
//
//				System.err.println(msg);
//			}
		} catch (Exception e) {
			e.printStackTrace();
			fail(e.getMessage());
		}

	}

}
