package org.biojava.bio.structure.io;

import static org.junit.Assert.*;

import java.io.IOException;

import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava3.structure.StructureIO;
import org.junit.Test;

/**
 * Testing for title parsing in PDB and mmCIF files
 * especially mult-line titles in PDB
 * 
 * @author duarte_j
 *
 */
public class TestTitleParsing {

	@Test
	public void test2W6E() throws IOException, StructureException {
		
		// an entry with a title in multiple lines in PDB file
		
		AtomCache cache = new AtomCache();
		
		StructureIO.setAtomCache(cache); 

		cache.setUseMmCif(false);
		Structure sPdb = StructureIO.getStructure("2W6E");
		cache.setUseMmCif(true);
		Structure sCif = StructureIO.getStructure("2W6E");
		
 
		// we can only compare titles by first forcing lower case, since the cases don't coincide cif vs pdb
		assertEquals(sPdb.getPDBHeader().getTitle().toLowerCase(),sCif.getPDBHeader().getTitle().toLowerCase());
		
		assertEquals(1, sPdb.getPDBHeader().getExperimentalTechniques().size());
		assertEquals(1, sCif.getPDBHeader().getExperimentalTechniques().size());
		
		assertEquals(sPdb.getPDBHeader().getResolution(),sCif.getPDBHeader().getResolution(),0.001);
		
		// an x-ray entry
		assertTrue(sPdb.isCrystallographic());
		assertTrue(sCif.isCrystallographic());

		assertFalse(sPdb.isNmr());
		assertFalse(sCif.isNmr());

		
	}

}
