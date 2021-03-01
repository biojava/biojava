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
import org.biojava.nbio.structure.io.CifFileReader;
import org.biojava.nbio.structure.io.LocalPDBDirectory;
import org.biojava.nbio.structure.io.LocalPDBDirectory.FetchBehavior;
import org.biojava.nbio.structure.io.LocalPDBDirectory.ObsoleteBehavior;
import org.biojava.nbio.structure.io.PDBFileReader;
import org.biojava.nbio.structure.io.StructureFiletype;
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
	public void setUp() throws IOException {
		cache = new AtomCache();

		// Delete files which were cached in previous tests
		String[] uncacheIDs = new String[] {
				"1cmw", "1hhb","4hhb"
		};

		List<LocalPDBDirectory> readers = new ArrayList<LocalPDBDirectory>();
		readers.add(new CifFileReader(cache.getPath()) );
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
		assertEquals(4,s.getPolyChains().size());

		String name2 = "4hhb.C";
		String chainId2 = "C";
		s = cache.getStructure(name2);

		assertEquals(1,s.getPolyChains().size());
		// Chain name 'C' corresponds to three IDs: C (141 res), I (1 HEM), and M (59 water)
		assertEquals(3,s.getChains().size());
		Chain c = s.getPolyChainByPDB(chainId2);
		assertEquals(chainId2,c.getName());

		// Number of groups: Polymer + water + ligand
		assertEquals(141,c.getAtomLength());
		assertEquals(141, s.getChainByIndex(0).getAtomLength());
		assertEquals(1, s.getChainByIndex(1).getAtomLength());
		assertEquals(59, s.getChainByIndex(2).getAtomLength());

		// Colon separators removed in BioJava 4.1.0
		String name2b = "4hhb:A";
		try {
			s = cache.getStructure(name2b);
			fail("Invalid structure format");
		} catch(StructureException e) {
		}


		// Numeric chain IDs are allowed but deprecated.
		String name3 = "4hhb.1";
		String chainId3 = "B";
		s = cache.getStructure(name3);
		assertNotNull(s);
		assertEquals(1,s.getPolyChains().size());

		c = s.getPolyChainByPDB(chainId3);
		assertEquals(chainId3,c.getName());


		String name4 = "4hhb.A:10-20,B:10-20,C:10-20";
		s = cache.getStructure(name4);
		assertNotNull(s);

		assertEquals(3,s.getPolyChains().size());
		assertEquals(3,s.getChains().size());

		c = s.getPolyChainByPDB("B");
		assertEquals(11,c.getAtomLength());

		String name5 = "4hhb.(A:10-20,A:30-40)";
		s =cache.getStructure(name5);
		assertNotNull(s);

		assertEquals(1,s.getPolyChains().size() );
		// two chains: A (22 res), and G (1 HEM)
		assertEquals(2,s.getChains().size() );
		c = s.getPolyChainByPDB("A");
		assertEquals(22,c.getAtomLength());

		try {
			// This syntax used to work, since the first paren is treated as a separator
			String name6 = "4hhb(A:10-20,A:30-40)";
			s =cache.getStructure(name6);
			fail("A chain separator is required after the ID since 4.2.0");
		} catch(StructureException e) {}

		// Works since we detect a separator
		String name8 = "4hhb.(C)";
		s = cache.getStructure(name8);

		assertEquals(1,s.getPolyChains().size());
		c = s.getPolyChainByPDB(chainId2);
		assertEquals(chainId2,c.getName());

	}

	@Test(expected=IOException.class)
	public void testObsoleteId() throws StructureException, IOException {
		cache.setFetchBehavior(FetchBehavior.FETCH_FILES);
		cache.setObsoleteBehavior(ObsoleteBehavior.THROW_EXCEPTION);

		// OBSOLETE PDB; should throw an exception
		cache.setFiletype(StructureFiletype.PDB);
		cache.getStructure("1HHB");
	}

	// note: we expect an IOException because 1CMW is obsolete and hasn't got a replacement
	@Test
	public void testFetchCurrent1CMW() throws IOException, StructureException {

		cache.setFetchBehavior(FetchBehavior.FETCH_FILES);
		cache.setObsoleteBehavior(ObsoleteBehavior.FETCH_CURRENT);

		// OBSOLETE PDB; should throw an exception
		cache.setFiletype(StructureFiletype.PDB);
		try {
			cache.getStructure("1CMW");
			fail("Obsolete structure should throw exception");
		} catch(IOException e) {}

		cache.setFiletype(StructureFiletype.CIF);
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

		cache.setFiletype(StructureFiletype.PDB);
		Structure s = cache.getStructure("1HHB");
		assertEquals("Failed to get the current ID for 1HHB.","4HHB",s.getPDBCode());

		cache.setFiletype(StructureFiletype.CIF);
		s = cache.getStructure("1HHB");
		assertEquals("Failed to get the current ID for 1HHB.","4HHB",s.getPDBCode());
	}

	// Fetching obsolete directly
	@Test
	public void testFetchObsolete() throws IOException, StructureException {
		cache.setFetchBehavior(FetchBehavior.FETCH_FILES);
		cache.setObsoleteBehavior(ObsoleteBehavior.FETCH_OBSOLETE);

		Structure s;
		cache.setFiletype(StructureFiletype.PDB);
		s = cache.getStructure("1CMW");
		assertEquals("Failed to get OBSOLETE file 1CMW.","1CMW", s.getPDBCode());

		s = cache.getStructure("1HHB");
		assertEquals("Failed to get OBSOLETE file 1HHB.","1HHB", s.getPDBCode());

		cache.setFiletype(StructureFiletype.CIF);
		s = cache.getStructure("1CMW");
		assertEquals("Failed to get OBSOLETE file 1CMW.","1CMW", s.getPDBCode());

		s = cache.getStructure("1HHB");
		assertEquals("Failed to get OBSOLETE file 1HHB.","1HHB", s.getPDBCode());

	}

	@Test
	public void testGetScopDomain() throws IOException, StructureException {
		String name = "d2gs2a_";

		Structure s = cache.getStructure(name);
		assertNotNull("Failed to fetch structure from SCOP ID", s);
		assertEquals("2gs2.A", s.getName());
	}

	@Test
	public void testSettingFileParsingType(){
		AtomCache cache = new AtomCache();

		//test defaults
		// first is mmtf, second is mmcif
		testFlags(cache, false, false, true);

		// now change the values
		cache.setFiletype(StructureFiletype.CIF);
		testFlags(cache, false, true, false);

		cache.setFiletype(StructureFiletype.MMTF);
		testFlags(cache, true, false, false);

		// this sets to use PDB!
		cache.setFiletype(StructureFiletype.PDB);
		testFlags(cache, false, false, false);

		// back to MMTF
		cache.setFiletype(StructureFiletype.MMTF);
		testFlags(cache, true, false, false);

		// back to parsing PDB
		cache.setFiletype(StructureFiletype.PDB);
		testFlags(cache, false, false, false);
	}


	/** test the flags for parsing in the atom cache
	 *
	 * @param cache
	 * @param useMmTf
	 * @param useMmCif
	 */
	private void testFlags(AtomCache cache ,boolean useMmTf, boolean useMmCif, boolean useBcif) {
		assertEquals("flag for parsing mmtf is set to " + cache.getFiletype() + " but should be " + useMmTf,
				cache.getFiletype() == StructureFiletype.MMTF, useMmTf);
		assertEquals("flag for parsing mmcif is set to " + cache.getFiletype() + " but should be set to " + useMmCif,
				cache.getFiletype() == StructureFiletype.CIF, useMmCif);
		assertEquals("flag for parsing bcif is set to " + cache.getFiletype() + " but should be set to " + useBcif,
				cache.getFiletype() == StructureFiletype.BCIF, useBcif);
	}
}
