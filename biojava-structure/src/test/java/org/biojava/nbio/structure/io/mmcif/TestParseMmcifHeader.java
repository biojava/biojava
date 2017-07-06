package org.biojava.nbio.structure.io.mmcif;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import org.biojava.nbio.structure.PDBHeader;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.junit.Ignore;
import org.junit.Test;

/**
 * Test parsing header information from MmCif files.
 * 
 * @author Anthony Bradley
 * @author Aleix Lafita
 *
 */
public class TestParseMmcifHeader {

	/**
	 * Test we can parse R-work and R-free effectively.
	 */
	@Test
	public void testRfactors() throws IOException, StructureException {

		AtomCache atomCache = new AtomCache();
		atomCache.setUseMmCif(true);

		Structure structure = atomCache.getStructure("4cup");

		PDBHeader pdbHeader = structure.getPDBHeader();
		// Check they are the same
		assertEquals(pdbHeader.getRfree(), 0.2078f, 0.000001f);
		assertEquals(pdbHeader.getRwork(), 0.1763f, 0.000001f);

	}

	/**
	 * Test parsing dates from MMCIF file version 4.
	 */
	@Test
	public void testDatesV4() throws IOException, StructureException {
		
		ClassLoader classLoader = this.getClass().getClassLoader();
		String file4 = classLoader.getResource("org/biojava/nbio/structure/io/mmcif/1stp_v4.cif").getPath();

		Structure s = StructureIO.getStructure(file4);

		// The latest modified year should be 2011
		assertEquals(s.getPDBHeader().getModDate().getYear() + 1900, 2011);

		// The release date should be october 1992
		assertEquals(s.getPDBHeader().getRelDate().getMonth() + 1, 10);
		assertEquals(s.getPDBHeader().getRelDate().getYear() + 1900, 1992);
		
		// The deposition date should be march 1992
		assertEquals(s.getPDBHeader().getDepDate().getMonth() + 1, 3);
		assertEquals(s.getPDBHeader().getDepDate().getYear() + 1900, 1992);
		
	}
	
	/**
	 * Test parsing dates from MMCIF file version 5.
	 */
	@Test @Ignore
	public void testDatesV5() throws IOException, StructureException {
		
		ClassLoader classLoader = this.getClass().getClassLoader();
		String file4 = classLoader.getResource("org/biojava/nbio/structure/io/mmcif/1stp_v5.cif").getPath();

		Structure s = StructureIO.getStructure(file4);

		// The latest modified year should be 2011
		assertEquals(s.getPDBHeader().getModDate().getYear() + 1900, 2011);

		// The release date should be october 1992
		assertEquals(s.getPDBHeader().getRelDate().getMonth() + 1, 10);
		assertEquals(s.getPDBHeader().getRelDate().getYear() + 1900, 1992);
		
		// The deposition date should be march 1992
		assertEquals(s.getPDBHeader().getDepDate().getMonth() + 1, 3);
		assertEquals(s.getPDBHeader().getDepDate().getYear() + 1900, 1992);
		

	}

}
