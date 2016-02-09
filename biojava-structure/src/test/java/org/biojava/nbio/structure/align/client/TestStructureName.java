package org.biojava.nbio.structure.align.client;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;

import org.biojava.nbio.structure.StructureException;
import org.junit.Test;

public class TestStructureName {

	@Test
	public void testMultiCharChainIds() throws StructureException {
		
		String str = "4V4F.AL";
		
		StructureName sn = new StructureName(str);
		
		assertEquals("AL", sn.getChainId());
		assertEquals("4V4F", sn.getPdbId());
		
		str = "4v4f.AL";
		sn = new StructureName(str);
		
		assertEquals("AL", sn.getChainId());
		assertEquals("4V4F", sn.getPdbId());

		str = "4v4f.al";
		sn = new StructureName(str);
		
		assertEquals("al", sn.getChainId());
		assertEquals("4V4F", sn.getPdbId());

		
		str = "4v4f.ABCD";
		sn = new StructureName(str);
		
		assertEquals("ABCD", sn.getChainId());
		assertEquals("4V4F", sn.getPdbId());

		
		// More than 4 characters should work too. In principle there's no limit in mmCIF, though the PDB is 
		// restricting chain ids to 4 chars 
		str = "4v4f.ABCDEFGHIJ";
		sn = new StructureName(str);
		
		assertEquals("ABCDEFGHIJ", sn.getChainId());
		assertEquals("4V4F", sn.getPdbId());


	}
	
	@Test
	public void testSingleCharChainIds() throws StructureException {
		
		String str = "1SMT.A";
		
		StructureName sn = new StructureName(str);
		
		assertEquals("A", sn.getChainId());
		assertEquals("1SMT", sn.getPdbId());
		
		str = "1SMT.a";
		sn = new StructureName(str);
		
		assertEquals("a", sn.getChainId());
		assertEquals("1SMT", sn.getPdbId());

		
	}

	@Test
	public void testFiles() throws IOException {
		File f = new File("hopefully_this_doesnt_exist_nalkjas3");
		assertFalse(f.exists());
		assertNull(f.getParentFile());
		StructureName sn = new StructureName("hopefully_this_doesnt_exist_nalkjas3");
		assertFalse(sn.isFile());
		assertTrue(sn.isPdbId());
	}

}
