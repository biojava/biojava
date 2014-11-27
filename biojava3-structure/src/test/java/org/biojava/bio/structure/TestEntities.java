package org.biojava.bio.structure;

import static org.junit.Assert.*;

import java.io.IOException;

import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.io.FileParsingParameters;
import org.junit.Test;

public class TestEntities {

	@Test
	public void testEntities() throws IOException, StructureException {
		
		AtomCache cache = new AtomCache();		
		FileParsingParameters params = new FileParsingParameters();
		params.setAlignSeqRes(true);
		cache.setFileParsingParams(params);
		
		// a structure with 6 identical chains in AU
		Structure s = cache.getStructure("3hbx");
		
		assertEquals(1,s.getEntities().size());
		
		Chain firstChain = s.getChainByPDB("A");
		
		
		for (Chain chain: s.getChains()) {
			assertTrue(s.getEntity(chain.getChainID()).getRepresentative() == firstChain);
		}
	}
	

}
