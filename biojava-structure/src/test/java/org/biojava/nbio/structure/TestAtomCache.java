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
 * Created on Mar 1, 2010
 * Author: Andreas Prlic
 *
 */

package org.biojava.nbio.structure;

import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.io.LocalPDBDirectory;
import org.biojava.nbio.structure.io.LocalPDBDirectory.FetchBehavior;
import org.biojava.nbio.structure.io.LocalPDBDirectory.ObsoleteBehavior;
import org.biojava.nbio.structure.io.MMCIFFileReader;
import org.biojava.nbio.structure.io.PDBFileReader;
import org.junit.Before;
import org.junit.Test;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.*;

// TODO dmyersturnbull: we should merge TestAtomCache and AtomCacheTest
public class TestAtomCache {

	public static final String lineSplit = System.getProperty("file.separator");
	private AtomCache cache;

	@Before
	public void setUp() {
		cache = new AtomCache();

		// Delete files which were cached in previous tests
		String[] uncacheIDs = new String[] {
				"1cmw", "1hhb","4hhb"
		};

		List<LocalPDBDirectory> readers = new ArrayList<LocalPDBDirectory>();
		readers.add(new MMCIFFileReader(cache.getPath()) );
		readers.add(new PDBFileReader(cache.getPath()) );
		for(LocalPDBDirectory reader : readers) {
			reader.setFetchBehavior(cache.getFetchBehavior());
			reader.setObsoleteBehavior(cache.getObsoleteBehavior());

			for(String pdbId : uncacheIDs) {
				reader.deleteStructure(pdbId);
			}
		}
	}

	// TODO dmyersturnbull: Which of these syntaxes do we support? We should re-enable after
	@Test
	public void testAtomCacheNameParsing() throws IOException, StructureException {


		String name1= "4hhb";
		Structure s = cache.getStructure(name1);
		assertNotNull(s);
		assertTrue(s.getPolyChains().size() == 4);

		String name2 = "4hhb.C";
		String chainId2 = "C";
		s = cache.getStructure(name2);

		assertTrue(s.getChains().size() == 1);
		Chain c = s.getChainByPDB(chainId2);
		assertEquals(c.getName(),chainId2);


		// Colon separators removed in BioJava 4.1.0
		String name2b = "4hhb:A";
		try {
			s = cache.getStructure(name2b);
			fail("Invalid structure format");
		} catch(IOException e) {
		} catch(StructureException e) {
		}


		// Numeric chain IDs are allowed but deprecated.
		String name3 = "4hhb.1";
		String chainId3 = "B";
		s = cache.getStructure(name3);
		assertNotNull(s);
		assertTrue(s.getChains().size() == 1);

		c = s.getChainByPDB(chainId3);
		assertEquals(c.getChainID(),chainId3);


		String name4 = "4hhb.A:10-20,B:10-20,C:10-20";
		s = cache.getStructure(name4);
		assertNotNull(s);

		assertEquals(s.getChains().size(), 3);

		c = s.getChainByPDB("B");
		assertEquals(c.getAtomLength(),11);

		String name5 = "4hhb.(A:10-20,A:30-40)";
		s =cache.getStructure(name5);
		assertNotNull(s);

		assertEquals(s.getChains().size(),1 );
		c = s.getChainByPDB("A");
		assertEquals(c.getAtomLength(),22);

		try {
			// This syntax used to work, since the first paren is treated as a separator
			String name6 = "4hhb(A:10-20,A:30-40)";
			s =cache.getStructure(name6);
			fail("A chain separator is required after the ID since 4.2.0");
		} catch(StructureException e) {}

		// Works since we detect a separator
		String name8 = "4hhb.(C)";
		s = cache.getStructure(name8);

		assertTrue(s.getChains().size() == 1);
		c = s.getChainByPDB(chainId2);
		assertEquals(c.getChainID(),chainId2);

	}

	@Test(expected=IOException.class)
	public void testObsoleteId() throws StructureException, IOException {
		cache.setFetchBehavior(FetchBehavior.FETCH_FILES);
		cache.setObsoleteBehavior(ObsoleteBehavior.THROW_EXCEPTION);

		// OBSOLETE PDB; should throw an exception
		cache.setUseMmCif(false);
		cache.getStructure("1HHB");
	}

	// note: we expect an IOException because 1CMW is obsolete and hasn't got a replacement
	@Test
	public void testFetchCurrent1CMW() throws IOException, StructureException {

		cache.setFetchBehavior(FetchBehavior.FETCH_FILES);
		cache.setObsoleteBehavior(ObsoleteBehavior.FETCH_CURRENT);

		// OBSOLETE PDB; should throw an exception
		cache.setUseMmCif(false);
		try {
			cache.getStructure("1CMW");
			fail("Obsolete structure should throw exception");
		} catch(IOException e) {}

		cache.setUseMmCif(true);
		try {
			cache.getStructure("1CMW");
			fail("Obsolete structure should throw exception");
		} catch(IOException e) {}
	}

	// 1HHB is obsolete with a replacement
	@Test
	public void testFetchCurrent1HHB() throws IOException, StructureException {

		cache.setFetchBehavior(FetchBehavior.FETCH_FILES);
		cache.setObsoleteBehavior(ObsoleteBehavior.FETCH_CURRENT);

		cache.setUseMmCif(false);
		Structure s = cache.getStructure("1HHB");
		assertEquals("Failed to get the current ID for 1HHB.","4HHB",s.getPDBCode());

		cache.setUseMmCif(true);
		s = cache.getStructure("1HHB");
		assertEquals("Failed to get the current ID for 1HHB.","4HHB",s.getPDBCode());
	}

	// Fetching obsolete directly
	@Test
	public void testFetchObsolete() throws IOException, StructureException {
		cache.setFetchBehavior(FetchBehavior.FETCH_FILES);
		cache.setObsoleteBehavior(ObsoleteBehavior.FETCH_OBSOLETE);

		Structure s;

		cache.setUseMmCif(false);
		s = cache.getStructure("1CMW");
		assertEquals("Failed to get OBSOLETE file 1CMW.","1CMW", s.getPDBCode());

		s = cache.getStructure("1HHB");
		assertEquals("Failed to get OBSOLETE file 1HHB.","1HHB", s.getPDBCode());

		cache.setUseMmCif(true);
		s = cache.getStructure("1CMW");
		assertEquals("Failed to get OBSOLETE file 1CMW.","1CMW", s.getPDBCode());

		s = cache.getStructure("1HHB");
		assertEquals("Failed to get OBSOLETE file 1HHB.","1HHB", s.getPDBCode());

	}

}
