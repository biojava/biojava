package org.biojava.nbio.structure.io.mmcif;

import static org.junit.Assert.assertEquals;

import java.io.IOException;

import org.biojava.nbio.structure.PDBHeader;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.junit.Test;

/**
 * A class to test the parsing of R-work from files
 * @author Anthony Bradley
 *
 */
public class TestParseHeader {
	
	/**
	 * Test we can parse R-work and R-free effectively.
	 * @throws IOException
	 * @throws StructureException
	 */
	@Test
	public void testRfactors() throws IOException, StructureException{
		
		AtomCache atomCache = new AtomCache();
		atomCache.setUseMmCif(true);
		StructureIO.setAtomCache(atomCache);
		Structure structure = StructureIO.getStructure("4cup");
		PDBHeader pdbHeader = structure.getPDBHeader();
		// Check they are the same
		assertEquals(pdbHeader.getRfree(),0.2078f, 0.000001f);
		assertEquals(pdbHeader.getRwork(),0.1763f, 0.000001f);
		
	}

}
