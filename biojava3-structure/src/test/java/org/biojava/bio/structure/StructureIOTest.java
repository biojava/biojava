package org.biojava.bio.structure;

import junit.framework.TestCase;

import org.biojava3.structure.StructureIO;

public class StructureIOTest extends TestCase {
	public void testStructureIO(){
		
		String pdbId = "1gav";
		
		int nrAssembls = StructureIO.getNrBiologicalAssemblies(pdbId);
		assertEquals(1,nrAssembls);
		
		pdbId = "1hv4";
		nrAssembls = StructureIO.getNrBiologicalAssemblies(pdbId);
		assertEquals(2,nrAssembls);
		
	}
}
