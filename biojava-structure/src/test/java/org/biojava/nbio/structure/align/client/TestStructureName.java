package org.biojava.nbio.structure.align.client;

import static org.junit.Assert.*;
import static org.biojava.nbio.structure.align.client.StructureName.Source.*;

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
	/**
	 * Test explicit prefixes
	 * @throws StructureException
	 */
	@Test
	public void testPrefixes() throws StructureException {
		StructureName sn;
		
		// Basic case
		sn = new StructureName("PDB:4hhb");
		assertTrue(sn.isPdbId());
		assertTrue(sn.getSource() == PDB);
		assertEquals("4HHB",sn.getPdbId());
		sn = new StructureName("PDB:4hhb.A:1-50");
		assertTrue(sn.isPdbId());
		assertTrue(sn.getSource() == PDB);
		assertEquals("4HHB",sn.getPdbId());
		// Invalid strings work too, they just don't load
		sn = new StructureName("PDB:x");
		assertTrue(sn.isPdbId());
		assertTrue(sn.getSource() == PDB);
		assertEquals("X",sn.getPdbId());
		// SCOP
		sn = new StructureName("SCOP:d2gs2a_");
		assertTrue(sn.isScopName());
		assertTrue(sn.getSource() == SCOP);
		assertEquals("2GS2",sn.getPdbId());
		// CATH
		sn = new StructureName("CATH:1qvrC03");
		assertTrue(sn.isCathID());
		assertTrue(sn.getSource() == CATH);
		assertEquals("1QVR",sn.getPdbId());
		// PDP
		sn = new StructureName("PDP:4HHBAa");
		assertTrue(sn.isPDPDomain());
		assertTrue(sn.getSource() == PDP);
		assertEquals("4HHB",sn.getPdbId());
		// URL
		sn = new StructureName("URL:http://www.rcsb.org/pdb/files/1B8G.pdb.gz");
		assertTrue(sn.isURL());
		assertTrue(sn.getSource() == URL);
		assertEquals("1B8G",sn.getPdbId());
		sn = new StructureName("URL:file:///4hhb.pdb");
		assertTrue(sn.isURL());
		assertTrue(sn.getSource() == URL);
		assertEquals("4HHB",sn.getPdbId());
//		// File
//		sn = new StructureName("FILE:~/4hhb.pdb");
//		assertTrue(sn.isFile());
//		assertTrue(sn.getSource() == FILE);
//		assertEquals("4HHB",sn.getPdbId());
//		// files are slightly different from URLs
//		sn = new StructureName("file:/4hhb.pdb");
//		assertTrue(sn.isFile());
//		assertTrue(sn.getSource() == FILE);
//		assertEquals("4HHB",sn.getPdbId());
//		// files are slightly different from URLs
//		sn = new StructureName("file:/4hhb_other.pdb");
//		assertTrue(sn.isFile());
//		assertTrue(sn.getSource() == FILE);
//		assertEquals("4hhb_other",sn.getPdbId());

		// ECOD
		sn = new StructureName("e1lyw.1");
		assertTrue(sn.isEcodDomain());
		assertTrue(sn.getSource() == ECOD);
		assertEquals("1LYW",sn.getPdbId());
		// BIO
		sn = new StructureName("BIO:2ehz:1");
		assertTrue(sn.isBioAssembly());
		assertTrue(sn.getSource() == BIO);
		assertEquals("2EHZ",sn.getPdbId());
		
		// Invalid prefix
		sn = new StructureName("XXX:2ehz");
		assertTrue(sn.isPdbId());
		assertTrue(sn.getSource() == PDB);
		assertEquals("XXX:2EHZ",sn.getPdbId());

	}

}
