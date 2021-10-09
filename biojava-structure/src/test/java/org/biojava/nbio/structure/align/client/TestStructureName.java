/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 */
package org.biojava.nbio.structure.align.client;

import static org.biojava.nbio.structure.align.client.StructureName.Source.BIO;
import static org.biojava.nbio.structure.align.client.StructureName.Source.CATH;
import static org.biojava.nbio.structure.align.client.StructureName.Source.ECOD;
import static org.biojava.nbio.structure.align.client.StructureName.Source.FILE;
import static org.biojava.nbio.structure.align.client.StructureName.Source.PDB;
import static org.biojava.nbio.structure.align.client.StructureName.Source.SCOP;
import static org.biojava.nbio.structure.align.client.StructureName.Source.URL;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URL;

import org.biojava.nbio.structure.PDBId;
import org.biojava.nbio.structure.StructureException;
import org.junit.Ignore;
import org.junit.Test;

public class TestStructureName {

	@Test
	public void testMultiCharChainIds() throws StructureException {

		String str = "4V4F.AL";

		StructureName sn = new StructureName(str);

		assertEquals("AL", sn.getChainId());
		assertEquals(new PDBId("4V4F"), sn.getPDBId());

		str = "4v4f.AL";
		sn = new StructureName(str);

		assertEquals("AL", sn.getChainId());
		assertEquals(new PDBId("4V4F"), sn.getPDBId());

		str = "4v4f.al";
		sn = new StructureName(str);

		assertEquals("al", sn.getChainId());
		assertEquals(new PDBId("4V4F"), sn.getPDBId());


		str = "4v4f.ABCD";
		sn = new StructureName(str);

		assertEquals("ABCD", sn.getChainId());
		assertEquals(new PDBId("4V4F"), sn.getPDBId());


		// More than 4 characters should work too. In principle there's no limit in mmCIF, though the PDB is
		// restricting chain ids to 4 chars
		str = "4v4f.ABCDEFGHIJ";
		sn = new StructureName(str);

		assertEquals("ABCDEFGHIJ", sn.getChainId());
		assertEquals(new PDBId("4V4F"), sn.getPDBId());


	}

	@Test
	public void testSingleCharChainIds() throws StructureException {

		String str = "1SMT.A";

		StructureName sn = new StructureName(str);

		assertEquals("A", sn.getChainId());
		assertEquals(new PDBId("1SMT"), sn.getPDBId());

		str = "1SMT.a";
		sn = new StructureName(str);

		assertEquals("a", sn.getChainId());
		assertEquals(new PDBId("1SMT"), sn.getPDBId());


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
		assertEquals(new PDBId("4HHB"),sn.getPDBId());
		sn = new StructureName("PDB:4hhb.A:1-50");
		assertTrue(sn.isPdbId());
		assertTrue(sn.getSource() == PDB);
		assertEquals(new PDBId("4HHB"),sn.getPDBId());
		// Invalid strings work too, they just don't load
		sn = new StructureName("PDB:x");
		assertTrue(sn.isPdbId());
		assertTrue(sn.getSource() == PDB);
		assertEquals(null, sn.getPDBId());
		// SCOP
		sn = new StructureName("SCOP:d2gs2a_");
		assertTrue(sn.isScopName());
		assertTrue(sn.getSource() == SCOP);
		assertEquals(new PDBId("2GS2"),sn.getPDBId());
		// CATH
		sn = new StructureName("CATH:1qvrC03");
		assertTrue(sn.isCathID());
		assertTrue(sn.getSource() == CATH);
		assertEquals(new PDBId("1QVR"),sn.getPDBId());
		// URL
		sn = new StructureName("URL:http://www.rcsb.org/pdb/files/1B8G.pdb.gz");
		assertTrue(sn.isURL());
		assertTrue(sn.getSource() == URL);
		assertEquals(new PDBId("1B8G"),sn.getPDBId());
		sn = new StructureName("URL:file:///4hhb.pdb");
		assertTrue(sn.isURL());
		assertTrue(sn.getSource() == URL);
		assertEquals(new PDBId("4HHB"),sn.getPDBId());
		// File: expand home directory (invalid URL)
		sn = new StructureName("FILE:~/4hhb.pdb");
		assertTrue(sn.isFile());
		assertTrue(sn.getSource() == FILE);
		assertEquals(new PDBId("4HHB"),sn.getPDBId());
		// Relative file (invalid URL)
		sn = new StructureName("file:4hhb.pdb");
		assertTrue(sn.isFile());
		assertTrue(sn.getSource() == FILE);
		assertEquals(new PDBId("4HHB"),sn.getPDBId());
		// Absolute paths are valid URLs
		sn = new StructureName("file:/4hhb_other.pdb");
		assertTrue(sn.isURL());
		assertTrue(sn.getSource() == URL);
		assertEquals(new PDBId("4HHB"),sn.getPDBId());

		// ECOD
		sn = new StructureName("e1lyw.1");
		assertTrue(sn.isEcodDomain());
		assertTrue(sn.getSource() == ECOD);
		assertEquals(new PDBId("1LYW"),sn.getPDBId());
		// BIO
		sn = new StructureName("BIO:2ehz:1");
		assertTrue(sn.isBioAssembly());
		assertTrue(sn.getSource() == BIO);
		assertEquals(new PDBId("2EHZ"),sn.getPDBId());

		// Invalid prefix
		sn = new StructureName("XXX:2ehz");
		assertTrue(sn.isPdbId());
		assertTrue(sn.getSource() == PDB);
//		assertEquals("XXX:2ehz",sn.getPDBId());
		assertEquals(null,sn.getPDBId()); // As long as it is EXPLICITLY defined to wrong type, it should NOT be recognized.

	}
	/**
	 * Test without prefixes
	 * @throws StructureException
	 */
	@Test
	public void testGuesses() throws StructureException {
		StructureName sn;

		// Basic case
		sn = new StructureName("4hhb");
		assertTrue(sn.isPdbId());
		assertTrue(sn.getSource() == PDB);
		assertEquals(new PDBId("4HHB"),sn.getPDBId());
		sn = new StructureName("4hhb.A:1-50");
		assertTrue(sn.isPdbId());
		assertTrue(sn.getSource() == PDB);
		assertEquals(new PDBId("4HHB"),sn.getPDBId());
		// Invalid strings should NOT work.
		sn = new StructureName("x");
		assertTrue(sn.isPdbId());
		assertTrue(sn.getSource() == PDB);
		assertEquals(null,sn.getPDBId());
		// SCOP
		sn = new StructureName("d2gs2a_");
		assertTrue(sn.isScopName());
		assertTrue(sn.getSource() == SCOP);
		assertEquals(new PDBId("2GS2"),sn.getPDBId());
		// CATH
		sn = new StructureName("1qvrC03");
		assertTrue(sn.isCathID());
		assertTrue(sn.getSource() == CATH);
		assertEquals(new PDBId("1QVR"),sn.getPDBId());
		// URL
		sn = new StructureName("http://www.rcsb.org/pdb/files/1B8G.pdb.gz");
		assertTrue(sn.isURL());
		assertTrue(sn.getSource() == URL);
		assertEquals(new PDBId("1B8G"),sn.getPDBId());
		sn = new StructureName("file:///4hhb.pdb");
		assertTrue(sn.isURL());
		assertTrue(sn.getSource() == URL);
		assertEquals(new PDBId("4HHB"),sn.getPDBId());


		// Files are hard to test, since they rely on existing files
		// You can run these tests locally after updating the hard-coded paths
		//sn = new StructureName("~/pdb/4hhb.pdb");
		//assertTrue(sn.isFile());
		//assertTrue(sn.getSource() == FILE);
		//assertEquals("4HHB",sn.getPdbId());
		//sn = new StructureName("/Users/blivens/pdb/4hhb.pdb");
		//assertTrue(sn.isFile());
		//assertTrue(sn.getSource() == FILE);
		//assertEquals("4HHB",sn.getPdbId());
		//sn = new StructureName("~/pdb/1f9m-1.pdb");
		//assertTrue(sn.isFile());
		//assertTrue(sn.getSource() == FILE);
		//assertEquals("1F9M",sn.getPdbId());

		// ECOD
		sn = new StructureName("e1lyw.1");
		assertTrue(sn.isEcodDomain());
		assertTrue(sn.getSource() == ECOD);
		assertEquals(new PDBId("1LYW"),sn.getPDBId());
		// BIO is not guessed
		sn = new StructureName("2ehz:1");
		assertFalse(sn.isBioAssembly());
		assertTrue(sn.getSource() == PDB);
		assertEquals(null, sn.getPDBId()); //Not Guessed
	}

	// Not really a test, but rather documenting Java's URL behavior
	@Ignore
	@Test
	public void testURLs() throws MalformedURLException {
		URL url;

		// Tilde doesn't get expanded
		url = new URL("file://~/1abc.pdb");
		assertEquals("/1abc.pdb", url.getPath());
		assertEquals("~",url.getHost());
		url = new URL("file:///~/1abc.pdb");
		assertEquals("/~/1abc.pdb", url.getPath());
		assertEquals("",url.getHost());

		// Supports omitting the initial slashes
		url = new URL("file:~/1abc.pdb");
		assertEquals("~/1abc.pdb", url.getPath());
		assertEquals("",url.getHost());

		// proper case. Three slashes gives empty host
		url = new URL("file:///1abc.pdb");
		assertEquals("/1abc.pdb", url.getPath());
		assertEquals("",url.getHost());

		// Two slashes triggers host
		url = new URL("file://1abc.pdb");
		assertEquals("", url.getPath());
		assertEquals("1abc.pdb",url.getHost());

		// One slash treated like zero slashes
		url = new URL("file:/1abc.pdb");
		assertEquals("/1abc.pdb", url.getPath());
		assertEquals("",url.getHost());
		assertEquals("file",url.getProtocol());

		// Surprise! url: prefix already works
		url = new URL("url:file://localhost/1abc.pdb");
		assertEquals("/1abc.pdb", url.getPath());
		assertEquals("localhost",url.getHost());
		assertEquals("file",url.getProtocol());
		url = new URL("URL:file://localhost/1abc.pdb");
		assertEquals("/1abc.pdb", url.getPath());
		assertEquals("localhost",url.getHost());
		assertEquals("file",url.getProtocol());

		// But doubling the file prefix doesn't. Is that OK?
		url = new URL("file:file://localhost/1abc.pdb");
		assertEquals("file://localhost/1abc.pdb", url.getPath());
		assertEquals("",url.getHost());
		assertEquals("file",url.getProtocol());

	}

}
