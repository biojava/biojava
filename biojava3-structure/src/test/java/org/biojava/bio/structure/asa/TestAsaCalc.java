package org.biojava.bio.structure.asa;

import java.io.IOException;

import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava3.structure.StructureIO;

import junit.framework.TestCase;

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
