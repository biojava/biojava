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
package org.biojava.nbio.structure.asa;

import junit.framework.TestCase;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;

import java.io.IOException;

/**
 * Testing of Accessible Surface Area calculations
 *
 *
 * @author duarte_j
 *
 */
public class TestAsaCalc extends TestCase {


	public void testAsa3PIU() throws StructureException, IOException {

		Structure structure = StructureIO.getStructure("3PIU");


		AsaCalculator asaCalc = new AsaCalculator(structure,
				AsaCalculator.DEFAULT_PROBE_SIZE,
				1000, 1, false);

		double totResidues = 0;
		double totAtoms = 0;

		GroupAsa[] groupAsas = asaCalc.getGroupAsas();

		double[] asas = asaCalc.calculateAsas();

		for (double asa:asas) {
			totAtoms += asa;
		}

		for (GroupAsa groupAsa: groupAsas) {
			totResidues+=groupAsa.getAsaU();

			assertTrue(groupAsa.getRelativeAsaU()<=1.0);
		}

		assertEquals(totAtoms, totResidues,0.000001);

		assertEquals(17638.0, totAtoms, 1.0);

	}
}
