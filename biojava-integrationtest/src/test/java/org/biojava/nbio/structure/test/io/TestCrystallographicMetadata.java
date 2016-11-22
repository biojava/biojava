package org.biojava.nbio.structure.test.io;

import org.junit.Test;
import static org.junit.Assert.*;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.align.util.AtomCache;

/**
 * A test for the parsing of some crystallographic metadata: non standard space group and non standard coordinate frame convention.
 *  
 * 
 * For more info see https://github.com/eppic-team/owl/issues/4 and https://github.com/eppic-team/eppic/issues/37
 * 
 * 
 * 
 * @author Jose Duarte
 * @since 4.2.5
 */
public class TestCrystallographicMetadata {

	
	@Test
	public void test4hhb() throws Exception {
		
		AtomCache cache = new AtomCache();
		// at the moment implemented only in mmcif
		cache.setUseMmCif(true); 
		StructureIO.setAtomCache(cache);

		Structure s = StructureIO.getStructure("4hhb");

		// 4hhb is one of the few entries that aren't in the standard coordinate frame convention
		assertTrue(s.getCrystallographicInfo().isNonStandardCoordFrameConvention());
		
		// 4hhn has a standard SG
		assertFalse(s.getCrystallographicInfo().isNonStandardSg());
		assertNotNull(s.getCrystallographicInfo().getSpaceGroup());
	}

	@Test
	public void test1smt() throws Exception {
		
		AtomCache cache = new AtomCache();
		// at the moment implemented only in mmcif
		cache.setUseMmCif(true); 
		StructureIO.setAtomCache(cache);

		Structure s = StructureIO.getStructure("1smt");

		// 1smt is a normal entry, should be standard
		assertFalse(s.getCrystallographicInfo().isNonStandardCoordFrameConvention());
	
		// 1smt has a standard SG
		assertFalse(s.getCrystallographicInfo().isNonStandardSg());
		assertNotNull(s.getCrystallographicInfo().getSpaceGroup());
		
	}
	
	@Test
	public void test1zna() throws Exception {
		AtomCache cache = new AtomCache();
		// at the moment implemented only in mmcif
		cache.setUseMmCif(true); 
		StructureIO.setAtomCache(cache);

		Structure s = StructureIO.getStructure("1zna");

		// 1zna is one of the few entries that has a non-standard SG
		assertTrue(s.getCrystallographicInfo().isNonStandardSg());
		assertNull(s.getCrystallographicInfo().getSpaceGroup());
	}
		

}
