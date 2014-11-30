package org.biojava.bio.structure;

import static org.junit.Assert.*;

import java.io.IOException;

import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.io.FileParsingParameters;
import org.junit.Test;

public class TestEntities {

	@Test
	public void testEntitiesMmCif() throws IOException, StructureException {
		
		AtomCache cache = new AtomCache();
		cache.setUseMmCif(true);
		FileParsingParameters params = new FileParsingParameters();
		params.setAlignSeqRes(true);
		cache.setFileParsingParams(params);
		
		// a structure with 6 identical chains in AU
		Structure s = cache.getStructure("3hbx");
		
		System.out.println("Entities for 3hbx: ");
		for (Entity ent:s.getEntities()) {
			System.out.print(ent.getRepresentative().getChainID()+":");
			for (Chain c:ent.getMembers()) {
				System.out.print(" "+c.getChainID());
			}
			System.out.println();
		}
		
		assertEquals(1,s.getEntities().size());		
		
		Chain firstChain = s.getChainByPDB("A");
		
		
		for (Chain chain: s.getChains()) {
			assertTrue(s.getEntity(chain.getChainID()).getRepresentative() == firstChain);
		}
	}
	
	@Test
	public void testEntitiesPdb() throws IOException, StructureException {
		
		AtomCache cache = new AtomCache();
		cache.setUseMmCif(false);
		FileParsingParameters params = new FileParsingParameters();
		params.setAlignSeqRes(true);
		cache.setFileParsingParams(params);
		
		// a structure with 6 identical chains in AU
		Structure s = cache.getStructure("3hbx");
		
		System.out.println("Entities for 3hbx: ");
		for (Entity ent:s.getEntities()) {
			System.out.print(ent.getRepresentative().getChainID()+":");
			for (Chain c:ent.getMembers()) {
				System.out.print(" "+c.getChainID());
			}
			System.out.println();
		}
		
		assertEquals(1,s.getEntities().size());
		
		Chain firstChain = s.getChainByPDB("A");
		
		
		for (Chain chain: s.getChains()) {
			assertTrue(s.getEntity(chain.getChainID()).getRepresentative() == firstChain);
		}
	}

	@Test
	public void testEntitiesPdbNoAlignSeqres() throws IOException, StructureException {
		
		AtomCache cache = new AtomCache();
		cache.setUseMmCif(false);
		FileParsingParameters params = new FileParsingParameters();
		params.setAlignSeqRes(false);
		cache.setFileParsingParams(params);
		
		// a structure with 6 identical chains in AU
		Structure s = cache.getStructure("3hbx");
		
		System.out.println("Entities for 3hbx: ");
		for (Entity ent:s.getEntities()) {
			System.out.print(ent.getRepresentative().getChainID()+":");
			for (Chain c:ent.getMembers()) {
				System.out.print(" "+c.getChainID());
			}
			System.out.println();
		}
		
		assertEquals(1,s.getEntities().size());
		
		Chain firstChain = s.getChainByPDB("A");
		
		
		for (Chain chain: s.getChains()) {
			assertTrue(s.getEntity(chain.getChainID()).getRepresentative() == firstChain);
		}
	}


}
