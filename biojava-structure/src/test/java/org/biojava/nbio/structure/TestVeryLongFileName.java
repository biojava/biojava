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
 */
package org.biojava.nbio.structure;

import org.junit.Assert;
import org.junit.Test;


public class TestVeryLongFileName {


	@Test
	public void testVeryLongFilename() {

		PDBHeader header = new PDBHeader();

		header.setTitle("IAMAVERYLONGTITLEWITHOUTANYSPACECHARACTERSJUSTTOMAKESUREWECANTESTWHATISGOINGTOHAPPENIFWETRYTOWRITETHISTOAPDBFILE.");

		header.toPDB();

		String title2 = "jCE V.1.1 : file:/Users/ap3/tmp/mareike/IAMAVERYLONGTITLEWITHOUTANYSPACECHARACTERSJUSTTOMAKESUREWECANTESTWHATISGOINGTOHAPPENIFWETRYTOWRITETHISTOAPDBFILE/pdb1_chainG.pdb vs. file:/Users/ap3/tmp/mareike/IAMAVERYLONGTITLEWITHOUTANYSPACECHARACTERSJUSTTOMAKESUREWECANTESTWHATISGOINGTOHAPPENIFWETRYTOWRITETHISTOAPDBFILE/pdb2_chainB.pdb ";
		header.setTitle(title2);
		header.toPDB();
	}
}
