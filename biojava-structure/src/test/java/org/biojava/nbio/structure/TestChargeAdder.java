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

import static org.junit.Assert.assertEquals;

import java.io.IOException;

import org.biojava.nbio.structure.io.ChargeAdder;
import org.junit.Test;

/**
 * Class of functions to test the charge adder.
 * @author Anthony Bradley
 *
 */
public class TestChargeAdder {

	
	/**
	 * Test that it works on a very basic level.
	 * @throws StructureException 
	 * @throws IOException 
	 */
	@Test
	public void testBasic() throws IOException, StructureException {
		
		// Get the structure
		Structure structure = StructureIO.getStructure("3AAE");
		ChargeAdder.addCharges(structure);
		// Now count the charges
		int chargeCount = 0;
		for(Chain chain : structure.getChains()) {
			for (Group group : chain.getAtomGroups()) {
				for (Atom atom : group.getAtoms()) {
					if (atom.getCharge()!=0) {
						chargeCount++;
					}
				}
			}
		}
		// Check that the count is as excpected
		assertEquals(425, chargeCount);	
	}
	
	
	/**
	 * Test that it can parse '?' values in the CCD.
	 * @throws StructureException 
	 * @throws IOException 
	 */
	@Test
	public void testQuestionMark() throws IOException, StructureException {
		// Get the structure
		Structure structure = StructureIO.getStructure("3PE6");
		ChargeAdder.addCharges(structure);
		// Now count the charges
		int chargeCount = 0;
		for(Chain chain : structure.getChains()) {
			for (Group group : chain.getAtomGroups()) {
				for (Atom atom : group.getAtoms()) {
					if (atom.getCharge()!=0) {
						chargeCount++;
					}
				}
			}
		}
		assertEquals(39, chargeCount);
	}
	
	
	
}
