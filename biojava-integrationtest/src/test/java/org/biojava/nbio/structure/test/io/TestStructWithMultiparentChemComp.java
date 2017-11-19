package org.biojava.nbio.structure.test.io;

import org.junit.Test;
import static org.junit.Assert.*;


import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.align.util.AtomCache;

public class TestStructWithMultiparentChemComp {
	
	@Test
	public void test4Q7U() throws Exception {
		
		AtomCache cache = new AtomCache();
		cache.setUseMmCif(true);
		StructureIO.setAtomCache(cache);
		
 		Structure s = StructureIO.getStructure("4q7u");
 		assertNotNull(s);
 		Chain c = s.getPolyChain("A");
 		
 		String seq = c.getSeqResSequence();
 		assertNotNull(seq);
 		assertEquals(245, seq.length());
 		assertFalse(seq.contains("?"));
 		assertTrue(seq.contains("X"));
		
	}

}
