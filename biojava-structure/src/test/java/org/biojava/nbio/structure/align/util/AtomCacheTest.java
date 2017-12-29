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
package org.biojava.nbio.structure.align.util;

import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.io.LocalPDBDirectory;
import org.biojava.nbio.structure.io.LocalPDBDirectory.FetchBehavior;
import org.biojava.nbio.structure.io.LocalPDBDirectory.ObsoleteBehavior;
import org.biojava.nbio.structure.io.MMCIFFileReader;
import org.biojava.nbio.structure.scop.ScopDatabase;
import org.biojava.nbio.structure.scop.ScopFactory;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.List;
import java.util.Locale;

import static org.junit.Assert.*;


/**
 * A test for {@link AtomCache}.
 * @author dmyerstu
 * @since 3.0.6
 */
public class AtomCacheTest {

	private AtomCache cache;
	private String previousPDB_DIR;

	@Before
	public void setUp() {
		previousPDB_DIR = System.getProperty(UserConfiguration.PDB_DIR, null);
		cache = new AtomCache();
		cache.setObsoleteBehavior(ObsoleteBehavior.FETCH_OBSOLETE);
		// Use a fixed SCOP version for stability
		ScopFactory.setScopDatabase(ScopFactory.VERSION_1_75B);
	}

	@After
	public void tearDown() {
		if (previousPDB_DIR != null)
			System.setProperty(UserConfiguration.PDB_DIR, previousPDB_DIR);
	}

	/**
	 * Tests {@link AtomCache#getStructureForDomain(String)} on a multi-chain domain with no ligands but an explicit range (not whole-chain).
	 */
	@Test
	public void testGetStructureForDomain1() throws IOException, StructureException {
		String ranges = "A:328-396,B:518-527";
		Structure whole = cache.getStructure("1h6w");
		AtomPositionMap map = new AtomPositionMap(StructureTools.getAllAtomArray(whole), AtomPositionMap.ANYTHING_MATCHER);
		List<ResidueRangeAndLength> rrs = ResidueRangeAndLength.parseMultiple(ranges, map);
		int expectedLengthA = rrs.get(0).getLength();
		int expectedLengthB = rrs.get(1).getLength();
		Structure structure = cache.getStructureForDomain("d1h6w.2");
		assertEquals(2, structure.getPolyChains().size());
		Chain a = structure.getPolyChainByPDB("A");
		Chain b = structure.getPolyChainByPDB("B");
		assertEquals(expectedLengthA, a.getAtomGroups().size());
		assertEquals(expectedLengthB, b.getAtomGroups().size());
	}

	/**
	 * Tests {@link AtomCache#getStructureForDomain(String)} on a multi-chain domain with two zinc ligands that occurs after the TER. The ligands are in chains E and F, so they should not be included in the domain.
	 */
	@Test
	public void testGetStructureForDomain2() throws IOException, StructureException {
		String ranges = "A:,B:";
		Structure whole = cache.getStructure("1I3O");
		AtomPositionMap map = new AtomPositionMap(StructureTools.getAllAtomArray(whole), AtomPositionMap.ANYTHING_MATCHER);
		List<ResidueRangeAndLength> rrs = ResidueRangeAndLength.parseMultiple(ranges, map);
		int expectedLengthA = rrs.get(0).getLength();
		int expectedLengthB = rrs.get(1).getLength();
		Structure structure = cache.getStructureForDomain("d1i3o.1");
		assertEquals(2, structure.getPolyChains().size());
		Chain a = structure.getPolyChainByPDB("A");
		Chain b = structure.getPolyChainByPDB("B");
		// since biojava 5.0 we have no ligand or water molecules in the polymer chains, we have to subtract the 3 water molecules
		assertEquals(expectedLengthA - 3, a.getAtomGroups().size());
		// since biojava 5.0 we have no ligand or water molecules in the polymer chains, we have to subtract the 4 water molecules
		assertEquals(expectedLengthB - 4, b.getAtomGroups().size());
		List<Group> ligandsA = StructureTools.filterLigands(b.getAtomGroups());
		assertEquals(0, ligandsA.size());
		List<Group> ligandsB = StructureTools.filterLigands(b.getAtomGroups());
		assertEquals(0, ligandsB.size());
	}

	/**
	 * Tests {@link AtomCache#getStructureForDomain(String)} on a single-chain domain with two zinc ligands that occurs after the TER.
	 */
	@Test
	public void testGetStructureForDomain3() throws IOException, StructureException {
		String ranges = "E:";
		Structure whole = cache.getStructure("1I3O");
		AtomPositionMap map = new AtomPositionMap(StructureTools.getAllAtomArray(whole), AtomPositionMap.ANYTHING_MATCHER);
		List<ResidueRangeAndLength> rrs = ResidueRangeAndLength.parseMultiple(ranges, map);
		int expectedLengthE = rrs.get(0).getLength();
		Structure structure = cache.getStructureForDomain("d1i3oe_");
		assertEquals(1, structure.getPolyChains().size());
		Chain e = structure.getPolyChainByPDB("E");
		// since biojava 5.0 we have no ligand molecules in the polymer chains, we have to subtract the 2 zinc molecules
		assertEquals(expectedLengthE - 2, e.getAtomGroups().size());

		Chain eligands = structure.getNonPolyChainsByPDB("E").get(0);
		List<Group> ligandsE = StructureTools.filterLigands(eligands.getAtomGroups());
		assertEquals(1, ligandsE.size());
	}

	/**
	 * Test parsing of chain-less ranges (present in SCOP < 1.73)
	 * @throws IOException
	 * @throws StructureException
	 */
	@Test
	public void testGetStructureForChainlessDomains() throws IOException, StructureException {
		ScopDatabase scop = ScopFactory.getSCOP(ScopFactory.VERSION_1_71); // Uses the range '1-135' without a chain
		Structure structure = cache.getStructureForDomain("d1hcy_1",scop);

		//System.out.println(cache.getStructure("1hcy"));
		//System.out.println(structure);
		assertEquals(1, structure.getPolyChains().size());
		Chain a = structure.getPolyChainByPDB("A");
		int expectedLengthA = 135;
		assertEquals(expectedLengthA, a.getAtomGroups().size());


		assertTrue(structure.hasNonPolyChain("G"));
		assertTrue(structure.hasNonPolyChain("H"));

		Chain copper  = structure.getNonPolyChain("I");
		assertEquals(1,copper.getAtomGroups().size());

	}

	@Test
	public void testSetPath_withTilde() throws Exception {
		cache.setPath("~" + File.separator);

		assertEquals(System.getProperty("user.home") + File.separator, cache.getPath());
	}

	@Test
	public void testNewInstanceWithTilder() throws Exception {
		AtomCache cache1 = new AtomCache("~" + File.separator);

		assertEquals(System.getProperty("user.home") + File.separator, cache1.getPath());
	}

	@Test
	public void testFetchBehavior() throws IOException, ParseException {
		// really more of a LocalPDBDirectory test, but throw it in with AtomCache
		String pdbId = "1hh0"; // A small structure, since we download it multiple times
		LocalPDBDirectory reader = new MMCIFFileReader(cache.getPath());

		// delete
		reader.deleteStructure(pdbId);
		assertNull("Failed to delete previous version",reader.getLocalFile(pdbId));

		// LOCAL_ONLY fails
		reader.setFetchBehavior(FetchBehavior.LOCAL_ONLY);
		Structure s;
		try {
			s = reader.getStructureById(pdbId);
			fail("LOCAL_ONLY shouldn't download files");
		} catch(IOException e) {
			assertTrue("Wrong IOException reason", e.getMessage().contains("configured not to download"));
		}

		// delete
		reader.deleteStructure(pdbId);
		assertNull("Failed to delete previous version",reader.getLocalFile(pdbId));

		// fetch from server
		reader.setFetchBehavior(FetchBehavior.FETCH_FILES);
		s = reader.getStructureById(pdbId);
		assertNotNull("Failed to fetch structure",s);
		File location = reader.getLocalFile(pdbId);

		long prerem = LocalPDBDirectory.LAST_REMEDIATION_DATE-1000*60*60*25; // 25 hours before the remediation
		location.setLastModified(prerem);
		assertEquals(prerem,location.lastModified()); //sanity check

		// force refetching
		reader.setFetchBehavior(FetchBehavior.FORCE_DOWNLOAD);
		s = reader.getStructureById(pdbId);
		assertNotNull("Failed to fetch structure",s);
		location = reader.getLocalFile(pdbId);
		assertTrue(location.exists());
		long currMod = location.lastModified();
		assertTrue("Not re-downloaded", currMod > prerem);

		// Now LOCAL_ONLY should work
		reader.setFetchBehavior(FetchBehavior.LOCAL_ONLY);
		s = reader.getStructureById(pdbId);
		assertNotNull("Failed to fetch structure",s);

		// Check remediation
		location.setLastModified(prerem);

		// Shouldn't re-fetch
		reader.setFetchBehavior(FetchBehavior.FETCH_FILES);
		s = reader.getStructureById(pdbId);
		location = reader.getLocalFile(pdbId);
		assertTrue(location.exists());
		assertEquals("Falsely re-downloaded", prerem,location.lastModified());

		// Now should re-fetch
		reader.setFetchBehavior(FetchBehavior.FETCH_REMEDIATED);
		s = reader.getStructureById(pdbId);
		assertNotNull("Failed to fetch structure",s);
		location = reader.getLocalFile(pdbId);
		assertTrue(location.exists());
		currMod = location.lastModified();
		assertTrue("Not re-downloaded", currMod > prerem);

		// test FETCH_IF_OUTDATED: change existing file timestamp to 2000 and try refetching (the file is from March 2009)
		SimpleDateFormat formatter = new SimpleDateFormat("yyyy/MM/dd", Locale.US);
		Date d = formatter.parse("2000/01/01");
		location.setLastModified(d.getTime());
		reader.setFetchBehavior(FetchBehavior.FETCH_IF_OUTDATED);
		s = reader.getStructureById(pdbId);
		assertNotNull("Failed to fetch structure",s);
		currMod = location.lastModified();
		assertTrue("Not re-downloaded", currMod>d.getTime());

		// try again: should not download
		reader.setFetchBehavior(FetchBehavior.FETCH_IF_OUTDATED);
		location = reader.getLocalFile(pdbId);
		currMod = location.lastModified();
		s = reader.getStructureById(pdbId);
		assertEquals("Falsely re-downloaded", currMod, location.lastModified());

	}

	@Test
	public void testSeqRes() throws StructureException, IOException {
		String name;
		StructureIdentifier id;
		Structure full, reduced;
		Chain chain;
		List<Group> seqres;

		// normal structure
		name = "1hh0";
		id = new SubstructureIdentifier(name);
		
		full = id.loadStructure(cache);
		assertEquals("Wrong number of models in full "+name,1,full.nrModels());
		assertEquals("Wrong number of chains in full "+name,1,full.getChains().size());
		chain = full.getChainByIndex(0);
		seqres = chain.getSeqResGroups();
		assertEquals("Wrong seqres length in full "+name,46,seqres.size());
		
		reduced = id.reduce(full);
		assertEquals("Wrong number of models in reduced "+name,1,reduced.nrModels());
		assertEquals("Wrong number of chains in reduced "+name,1,reduced.getChains().size());
		chain = reduced.getChainByIndex(0);
		seqres = chain.getSeqResGroups();
		assertEquals("Wrong seqres length in reduced "+name,46,seqres.size());

		// single chain
		name = "1hh0.A";
		id = new SubstructureIdentifier(name);
		
		full = id.loadStructure(cache);
		assertEquals("Wrong number of models in full "+name,1,full.nrModels());
		assertEquals("Wrong number of chains in full "+name,1,full.getChains().size());
		chain = full.getChainByIndex(0);
		seqres = chain.getSeqResGroups();
		assertEquals("Wrong seqres length in full "+name,46,seqres.size());
		
		reduced = id.reduce(full);
		assertEquals("Wrong number of models in reduced "+name,1,reduced.nrModels());
		assertEquals("Wrong number of chains in reduced "+name,1,reduced.getChains().size());
		chain = reduced.getChainByIndex(0);
		seqres = chain.getSeqResGroups();
		assertEquals("Wrong seqres length in reduced "+name,46,seqres.size());

		// subrange
		name = "1hh0.A:10-20";
		id = new SubstructureIdentifier(name);
		
		full = id.loadStructure(cache);
		assertEquals("Wrong number of models in full "+name,1,full.nrModels());
		assertEquals("Wrong number of chains in full "+name,1,full.getChains().size());
		chain = full.getChainByIndex(0);
		seqres = chain.getSeqResGroups();
		assertEquals("Wrong seqres length in full "+name,46,seqres.size());
		assertEquals("Wrong SeqNum at first group in full",1,(int)chain.getAtomGroup(0).getResidueNumber().getSeqNum());

		reduced = id.reduce(full);
		assertEquals("Wrong number of models in reduced "+name,1,reduced.nrModels());
		assertEquals("Wrong number of chains in reduced "+name,1,reduced.getChains().size());
		chain = reduced.getChainByIndex(0);
		seqres = chain.getSeqResGroups();
		assertEquals("Wrong seqres length in reduced "+name,46,seqres.size());
		
		assertEquals("Wrong SeqNum at first group in reduced",10,(int)chain.getAtomGroup(0).getResidueNumber().getSeqNum());

	}
	
}
