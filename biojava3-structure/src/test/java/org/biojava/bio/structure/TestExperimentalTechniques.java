package org.biojava.bio.structure;

import static org.junit.Assert.*;

import java.io.IOException;
import java.util.Iterator;

import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava3.structure.StructureIO;
import org.junit.Test;

public class TestExperimentalTechniques {

	@Test
	public void test4LNC() throws IOException, StructureException {
		
		// a multiple experimental techniques PDB entry (X-RAY + NEUTRON DIFFRACTION)
		
		AtomCache cache = new AtomCache("/tmp",false);
		
		StructureIO.setAtomCache(cache); 

		cache.setUseMmCif(false);
		Structure sPdb = StructureIO.getStructure("4LNC");
		cache.setUseMmCif(true);
		Structure sCif = StructureIO.getStructure("4LNC");
		
		comparePdbToCif(sPdb, sCif);
		
		assertEquals(2, sPdb.getPDBHeader().getExperimentalTechniques().size());
		assertEquals(2, sCif.getPDBHeader().getExperimentalTechniques().size());
		
		
		
		assertTrue(sPdb.isCrystallographic());
		assertTrue(sCif.isCrystallographic());

		assertFalse(sPdb.isNmr());
		assertFalse(sCif.isNmr());

		
	}
	
	@Test
	public void test3ZPK() throws IOException, StructureException {
		
		// a multiple experimental techniques PDB entry (EM + SOLUTION NMR)
		
		AtomCache cache = new AtomCache("/tmp",false);
		
		StructureIO.setAtomCache(cache); 

		cache.setUseMmCif(false);
		Structure sPdb = StructureIO.getStructure("3ZPK");
		cache.setUseMmCif(true);
		Structure sCif = StructureIO.getStructure("3ZPK");
		
		comparePdbToCif(sPdb, sCif);
		
		assertEquals(2, sPdb.getPDBHeader().getExperimentalTechniques().size());
		assertEquals(2, sCif.getPDBHeader().getExperimentalTechniques().size());
		
		
		
		assertFalse(sPdb.isCrystallographic());
		assertFalse(sCif.isCrystallographic());

		assertTrue(sPdb.isNmr());
		assertTrue(sCif.isNmr());

		
	}
	
	@Test
	public void test2B6O() throws IOException, StructureException {
		
		// a single experimental technique ELECTRON CRYSTALLOGRAPHY entry
		
		AtomCache cache = new AtomCache("/tmp",false);
		
		StructureIO.setAtomCache(cache); 

		cache.setUseMmCif(false);
		Structure sPdb = StructureIO.getStructure("2B6O");
		cache.setUseMmCif(true);
		Structure sCif = StructureIO.getStructure("2B6O");
		
		comparePdbToCif(sPdb, sCif);
		
		assertEquals(1, sPdb.getPDBHeader().getExperimentalTechniques().size());
		assertEquals(1, sCif.getPDBHeader().getExperimentalTechniques().size());
		
		
		
		assertTrue(sPdb.isCrystallographic());
		assertTrue(sCif.isCrystallographic());

		assertFalse(sPdb.isNmr());
		assertFalse(sCif.isNmr());

		
	}
	
	@Test
	public void test4CSO() throws IOException, StructureException {
		
		// a single experimental technique (X-RAY) entry
		
		AtomCache cache = new AtomCache("/tmp",false);
		
		StructureIO.setAtomCache(cache); 

		cache.setUseMmCif(false);
		Structure sPdb = StructureIO.getStructure("4CSO");
		cache.setUseMmCif(true);
		Structure sCif = StructureIO.getStructure("4CSO");
		
		comparePdbToCif(sPdb, sCif);
		
		assertEquals(1, sPdb.getPDBHeader().getExperimentalTechniques().size());
		assertEquals(1, sCif.getPDBHeader().getExperimentalTechniques().size());
		
		
		
		assertTrue(sPdb.isCrystallographic());
		assertTrue(sCif.isCrystallographic());

		assertFalse(sPdb.isNmr());
		assertFalse(sCif.isNmr());

		
	}
	
	private void comparePdbToCif(Structure sPdb, Structure sCif) {
		assertNotNull(sPdb.getPDBHeader().getExperimentalTechniques());
		assertNotNull(sCif.getPDBHeader().getExperimentalTechniques());
		
		Iterator<ExperimentalTechnique> itCif = sCif.getPDBHeader().getExperimentalTechniques().iterator();
		for (ExperimentalTechnique et:sPdb.getPDBHeader().getExperimentalTechniques()) {
			assertEquals(et, itCif.next());
		}
	}

}
