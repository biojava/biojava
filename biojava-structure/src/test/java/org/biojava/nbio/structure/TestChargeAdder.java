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
