package org.biojava.bio.structure.xtal;

import java.io.IOException;

import org.biojava.bio.structure.ExperimentalTechnique;
import org.biojava.bio.structure.PDBCrystallographicInfo;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava3.structure.StructureIO;
import org.junit.Test;

import static org.junit.Assert.*;

/**
 * Testing of crystallographic info parsing in both pdb and mmCIF files
 * 
 * @author duarte_j
 *
 */
public class TestCrystalInfo {

	private static final float DELTA = 0.000001f;

	@Test
	public void test1NMR() throws IOException, StructureException { 

		AtomCache cache = new AtomCache("/tmp",false);
		
		StructureIO.setAtomCache(cache); 

		cache.setUseMmCif(false);
		Structure s1 = StructureIO.getStructure("1NMR");
		assertFalse(s1.isCrystallographic());
		assertTrue(s1.isNmr());
		assertEquals(s1.getPDBHeader().getExperimentalTechniques().iterator().next(),ExperimentalTechnique.SOLUTION_NMR);
		
		cache.setUseMmCif(true);
		Structure s2 = StructureIO.getStructure("1NMR");
		assertFalse(s2.isCrystallographic());
		assertTrue(s2.isNmr());
		assertEquals(s2.getPDBHeader().getExperimentalTechniques().iterator().next(),ExperimentalTechnique.SOLUTION_NMR);


		testCrystallographicInfo(s1, s2);

	}
	
	@Test
	public void test1B8G() throws IOException, StructureException { 

		AtomCache cache = new AtomCache("/tmp",false);
		
		StructureIO.setAtomCache(cache); 

		cache.setUseMmCif(false);
		Structure s1 = StructureIO.getStructure("1B8G");
		assertTrue(s1.isCrystallographic());
		assertFalse(s1.isNmr());
		assertEquals(s1.getPDBHeader().getExperimentalTechniques().iterator().next(),ExperimentalTechnique.XRAY_DIFFRACTION);

		
		cache.setUseMmCif(true);
		Structure s2 = StructureIO.getStructure("1B8G");
		assertTrue(s2.isCrystallographic());
		assertFalse(s2.isNmr());
		assertEquals(s2.getPDBHeader().getExperimentalTechniques().iterator().next(),ExperimentalTechnique.XRAY_DIFFRACTION);

		testCrystallographicInfo(s1, s2);

	}
	
	@Test
	public void test4M7P() throws IOException, StructureException { 

		// multimodel x-ray structure
		
		AtomCache cache = new AtomCache("/tmp",false);
		
		StructureIO.setAtomCache(cache); 

		cache.setUseMmCif(false);
		Structure s1 = StructureIO.getStructure("4M7P");
		assertTrue(s1.isCrystallographic());
		assertFalse(s1.isNmr());
		assertTrue(s1.nrModels()>1);
		assertEquals(s1.getPDBHeader().getExperimentalTechniques().iterator().next(),ExperimentalTechnique.XRAY_DIFFRACTION);

		
		cache.setUseMmCif(true);
		Structure s2 = StructureIO.getStructure("4M7P");
		assertTrue(s2.isCrystallographic());
		assertFalse(s2.isNmr());
		assertTrue(s2.nrModels()>1);
		assertEquals(s2.getPDBHeader().getExperimentalTechniques().iterator().next(),ExperimentalTechnique.XRAY_DIFFRACTION);

		testCrystallographicInfo(s1, s2);

	}
	
	@Test
	public void test2MBQ() throws IOException, StructureException { 

		// single model NMR structure
		
		AtomCache cache = new AtomCache("/tmp",false);
		
		StructureIO.setAtomCache(cache); 

		cache.setUseMmCif(false);
		Structure s1 = StructureIO.getStructure("2MBQ");
		assertFalse(s1.isCrystallographic());
		assertTrue(s1.isNmr());
		assertFalse(s1.nrModels()>1);
		assertEquals(s1.getPDBHeader().getExperimentalTechniques().iterator().next(),ExperimentalTechnique.SOLUTION_NMR);

		
		cache.setUseMmCif(true);
		Structure s2 = StructureIO.getStructure("2MBQ");
		assertFalse(s2.isCrystallographic());
		assertTrue(s2.isNmr());
		assertFalse(s2.nrModels()>1);
		assertEquals(s2.getPDBHeader().getExperimentalTechniques().iterator().next(),ExperimentalTechnique.SOLUTION_NMR);


		testCrystallographicInfo(s1, s2);

	}
	
	private void testCrystallographicInfo(Structure s1, Structure s2) {
		PDBCrystallographicInfo xtalInfo = s1.getPDBHeader().getCrystallographicInfo();
		PDBCrystallographicInfo xtalInfo2 = s2.getPDBHeader().getCrystallographicInfo();
		
		
		if (xtalInfo==null && xtalInfo2==null) return; // both null: NMR or something similar, nothing to test

		CrystalCell cell1 = xtalInfo.getCrystalCell();
		CrystalCell cell2 = xtalInfo2.getCrystalCell();

		if (s1.isNmr()) assertTrue(cell1==null); 
		if (s2.isNmr()) assertTrue(cell2==null);
		
		if (cell1==null && cell2==null) return;
		
		assertEquals(xtalInfo.getA(), xtalInfo2.getA(),DELTA);
		assertEquals(xtalInfo.getB(), xtalInfo2.getB(),DELTA);
		assertEquals(xtalInfo.getC(), xtalInfo2.getC(),DELTA);
		assertEquals(xtalInfo.getAlpha(), xtalInfo2.getAlpha(),DELTA);
		assertEquals(xtalInfo.getBeta(), xtalInfo2.getBeta(),DELTA);
		assertEquals(xtalInfo.getGamma(), xtalInfo2.getGamma(),DELTA);
		assertEquals(xtalInfo.getZ(), xtalInfo2.getZ());

		assertEquals(xtalInfo.getSpaceGroup(),xtalInfo2.getSpaceGroup());
	}
}
