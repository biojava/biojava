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
package org.biojava.nbio.structure.test.io;

import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.io.LocalPDBDirectory.ObsoleteBehavior;
import org.biojava.nbio.structure.quaternary.BioAssemblyInfo;
import org.biojava.nbio.structure.xtal.CrystalCell;
import org.junit.After;
import org.junit.BeforeClass;
import org.junit.ComparisonFailure;
import org.junit.Ignore;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import javax.vecmath.Matrix4d;
import java.io.*;
import java.util.*;

import static org.junit.Assert.*;

/**
 * A test to make sure both PDB and mmCIF parsers can parse
 * properly large samples of the PDB.
 *
 * Will take very long to run, thus they are ignored by default.
 * To run them use, for the 1000 entries one:
 * <pre>
 * mvn -Dtest=TestLongPdbVsMmCifParsing#testLongPdbVsMmCif test
 * </pre>
 * or for the 10000 entries:
 * <pre>
 * mvn -Dtest=TestLongPdbVsMmCifParsing#testVeryLongPdbVsMmCif test
 * </pre>
 *
 *
 * @author duarte_j
 *
 */
public class TestLongPdbVsMmCifParsing {

	private static final Logger logger = LoggerFactory.getLogger(TestLongPdbVsMmCifParsing.class);

	private static final String TEST_LARGE_SET_FILE = "/random_1000_set.list";
	private static final String TEST_VERY_LARGE_SET_FILE = "/random_10000_set.list";

	private static final int DOTS_PER_LINE = 100;

	private static final float DELTA = 0.01f;
	private static final float DELTA_RESOLUTION = 0.01f;
	private static final float DELTA_RFREE = 0.01f;

	/**
	 * The maximum number of PDBs for which we allow a mismatch of mol_ids (entity_ids) between PDB and mmCIF files
	 * If more mismatches than this, the test will fail.
	 * As of 2014.12.04 there are 7 mismatches
	 */
	private static final int   MAX_ALLOWED_MOL_ID_MISMATCHES = 10;

	private static AtomCache cache;
	private static FileParsingParameters params;

	private String pdbId;
	
	private int countTested = 0;

	private HashSet<String> pdbIdsWithMismatchingMolIds;

	@BeforeClass
	public static void setUpBeforeClass() {
		cache = new AtomCache();

		System.out.println("##### Starting long test. THIS CAN TAKE UP TO 1 HOUR TO COMPLETE!");
		System.out.println("##### Using PDB/mmCIF cache dir: "+cache.getPath());
		System.out.println("##### Each dot is a PDB entry being tested. "+DOTS_PER_LINE+" dots per line");

		// disallow the use of the default /tmp dir, to make sure PDB_DIR is set
		if (cache.getPath().equals(System.getProperty("java.io.tmpdir")) ||
			(cache.getPath().equals(System.getProperty("java.io.tmpdir")+File.separator))    ) {

			throw new IllegalArgumentException("PDB_DIR has not been set or it is set to the default temp directory. Please set PDB_DIR to run this test");
		};

		params = new FileParsingParameters();
		cache.setFileParsingParams(params);
		cache.setObsoleteBehavior(ObsoleteBehavior.THROW_EXCEPTION);
	}

	@Ignore
	@Test
	public void testLongPdbVsMmCif() throws IOException, StructureException {

		List<String> pdbIds = readTestSetFile(TEST_LARGE_SET_FILE);

		testAll(pdbIds);

	}

	@Ignore
	@Test
	public void testVeryLongPdbVsMmCif() throws IOException, StructureException {

		List<String> pdbIds = readTestSetFile(TEST_VERY_LARGE_SET_FILE);

		testAll(pdbIds);

	}

	@Ignore
	@Test
	public void testSingle() throws IOException, StructureException {
		testAll(Arrays.asList("4kro"));
	}

	@After
	public void printInfo() {
		if (pdbId!=null)
			System.out.println("\n##### ----> Last tested PDB entry was: "+pdbId + " ("+ countTested + " done so far)");
	}

	private void testAll(List<String> pdbIds) throws IOException, StructureException {

		pdbIdsWithMismatchingMolIds = new HashSet<String>();

		long start = System.currentTimeMillis();

		System.out.println("##### Total of "+pdbIds.size()+" PDB entries to test");

		for (int i = 0; i<pdbIds.size(); i++) {
			pdbId = pdbIds.get(i);
			
			countTested = i + 1;
			
			System.out.print(".");

			testSingleEntry(pdbId);

			if ( ( (i+1)%DOTS_PER_LINE )==0 ) System.out.println();
			
		}

		pdbId = null; // to avoid printing the message if tests pass for all PDB entries

		long end = System.currentTimeMillis();

		checkWarnings();

		System.out.printf("\nDone in %5.1f minutes\n", (end-start)/60000.0);
	}

	private void checkWarnings() {
		if (pdbIdsWithMismatchingMolIds.size()>0)
			System.out.println("A total of "+pdbIdsWithMismatchingMolIds.size()+" PDB entries have mismatches in their Compound mol_ids (entity_ids)");

		assertTrue("Mismatching mol_id (entity_id) between pdb and cif above the maximum allowed ("+MAX_ALLOWED_MOL_ID_MISMATCHES+")",
				pdbIdsWithMismatchingMolIds.size()<MAX_ALLOWED_MOL_ID_MISMATCHES);
	}

	private void testSingleEntry(String pdbId) throws IOException, StructureException {

		Structure sCif = getCifStructure(pdbId);
		Structure sPdb = getPdbStructure(pdbId);

		assertNotNull(sCif);
		assertNotNull(sPdb);

		try {

			testStructureMethods(sPdb, sCif);

			testHeader(sPdb, sCif);

			testChains(sPdb, sCif);

		} catch (ComparisonFailure e) {
			System.out.println("\nComparison failure! Values follow:");
			System.out.println("Actual  : "+e.getActual());
			System.out.println("Expected: "+e.getExpected());
			throw e;
		}

	}

	private void testStructureMethods(Structure sPdb, Structure sCif) {

		assertEquals("failed isNmr:",sPdb.isNmr(), sCif.isNmr());
		assertEquals("failed isCrystallographic:",sPdb.isCrystallographic(), sCif.isCrystallographic());
		assertEquals("failed nrModels:",sPdb.nrModels(), sCif.nrModels());

		assertEquals("failed for getPdbCode:",sPdb.getPDBCode(),sCif.getPDBCode());

		assertFalse(sPdb.isBiologicalAssembly());
		assertFalse(sCif.isBiologicalAssembly());

		// TODO journal article not parsed in mmCIF parser
		//assertEquals("failed hasJournalArticle",sPdb.hasJournalArticle(),sCif.hasJournalArticle());

		// entity type should always be present
		for (EntityInfo e: sPdb.getEntityInfos()) {
			assertNotNull(e.getType());
		}

		for (EntityInfo e: sCif.getEntityInfos()) {
			assertNotNull(e.getType());
		}
		
		// entities: there's quite some inconsistencies here between pdb and cif:
		// sugar polymers are not in pdb at all: we avoid them
		boolean canCompareEntityCounts = true;
		for (EntityInfo e:sCif.getEntityInfos()) {
			if (e.getDescription().contains("SUGAR")) canCompareEntityCounts = false;
		}
		if (canCompareEntityCounts) {
			int entCountCif = 0;
			for (EntityInfo e: sCif.getEntityInfos()) {
				if (e.getType() == EntityType.POLYMER) 
					entCountCif++; 

			}
			int entCountPdb = 0;
			for (EntityInfo e:sPdb.getEntityInfos()) {
				if (e.getType() == EntityType.POLYMER) 
					entCountPdb++;
			}

			assertEquals("failed number of non-sugar polymeric Entities pdb vs cif", entCountPdb, entCountCif);
		}

		// ss bonds
		// 4ab9 contains an error in ssbond in pdb file (misses 1 ssbond)
		// 2bdi contains also errors, the counts in both differ a lot 80 vs 92
		if (!sPdb.getPDBCode().equals("4AB9") && !sPdb.getPDBCode().equals("2BDI"))
			assertEquals("number of ss bonds should coincide pdb vs cif", sPdb.getSSBonds().size(), sCif.getSSBonds().size());

	}

	private void testHeader(Structure sPdb, Structure sCif) {

		PDBHeader hPdb = sPdb.getPDBHeader();
		PDBHeader hCif = sCif.getPDBHeader();

		boolean isNmr = sPdb.isNmr();
		boolean isCrystallographic = sPdb.isCrystallographic();

		assertNotNull(hPdb);
		assertNotNull(hCif);

		assertEquals("failed for PDB id (getIdCode)",hPdb.getIdCode(),hCif.getIdCode());

		assertNotNull("pdb authors null",hPdb.getAuthors());
		assertNotNull("cif authors null",hCif.getAuthors());
		// I suppose 2 is a safe bet for authors length...
		assertTrue("authors length should be at least 2",hPdb.getAuthors().length()>=2);
		// for authors we strip spaces in case of ambiguities with names
		// there's too much variability in authors, commenting out, e.g. for 1zjo they don't coincide
		//assertEquals("failed getAuthors:",
		//		hPdb.getAuthors().toLowerCase().replaceAll(" ", ""),
		//		hCif.getAuthors().toLowerCase().replaceAll(" ", ""));

		assertNotNull("pdb classification null in pdb",hPdb.getClassification());
		assertNotNull("cif classification null in cif",hCif.getClassification());
		// there's too much variability in classification between pdb and mmcif, e.g. in 3ofb they don't coincide
		//assertEquals("failed getClassification:",hPdb.getClassification().toLowerCase(), hCif.getClassification().toLowerCase());

		// description is set in CIF parser to same as classification (_struct_keywords.pdbx_keywords field)
		// while in PDB parser it is simply not set
		//assertNotNull("pdb description null",hPdb.getDescription());
		assertNotNull("cif description null",hCif.getDescription());
		//assertEquals("failed getDescription:",hPdb.getDescription().toLowerCase(), hCif.getDescription().toLowerCase());

		assertEquals("failed getDepDate:",hPdb.getDepDate(), hCif.getDepDate());
		assertEquals("failed getModDate:",hPdb.getModDate(), hCif.getModDate());

		assertNotNull(hPdb.getExperimentalTechniques());
		assertNotNull(hCif.getExperimentalTechniques());
		assertTrue(hPdb.getExperimentalTechniques().size()>0);
		assertEquals("failed for getExperimentalTechniques",hPdb.getExperimentalTechniques(),hCif.getExperimentalTechniques());

		// for some Electron Microscopy/Crystallography entries (e.g. 3iz2) the resolution in mmCIF is not present in the usual place
		if (!hPdb.getExperimentalTechniques().contains(ExperimentalTechnique.ELECTRON_CRYSTALLOGRAPHY) &&
				!hPdb.getExperimentalTechniques().contains(ExperimentalTechnique.ELECTRON_MICROSCOPY)) {
			assertEquals("failed getResolution:",hPdb.getResolution(), hCif.getResolution(), DELTA_RESOLUTION);
		}

		// JRNL record is sometimes missing (e.g. 21bi) and thus is null, we can't test for nulls here in the general case
		//assertNotNull("journal article null",hPdb.getJournalArticle());
		// TODO journal article not parsed in mmCIF parser
		// TODO when fixed in mmCIF parser, compare PDB to mmCIF values if not null
		//assertNotNull("journal article null",hCif.getJournalArticle());

		assertNotNull("title null in pdb",hPdb.getTitle());
		assertNotNull("title null in cif",hCif.getTitle());
		// for titles we strip spaces in case of ambiguities with spacing
		assertEquals("failed for getTitle",
					hPdb.getTitle().toLowerCase().replaceAll(" ", ""),
					hCif.getTitle().toLowerCase().replaceAll(" ", ""));

		// tests specific to experimental techniques
		if (isNmr) {
			assertEquals("resolution is not the default value in NMR structure",
					PDBHeader.DEFAULT_RESOLUTION, hPdb.getResolution(), DELTA_RESOLUTION);
		}

		if (!isCrystallographic) {
			assertEquals("rfree is not the default value in non-crystallographic structure in pdb",
					PDBHeader.DEFAULT_RFREE, DELTA_RFREE, hPdb.getRfree());
			assertEquals("rfree is not the default value in non-crystallographic structure in cif",
					PDBHeader.DEFAULT_RFREE, DELTA_RFREE, hCif.getRfree());
		}

		if (isCrystallographic) {

			assertEquals("failed for Rfree:",hPdb.getRfree(), hCif.getRfree(), DELTA_RFREE);

			assertNotNull("getCrystallographicInfo is null in pdb",hPdb.getCrystallographicInfo());
			assertNotNull("getCrystallographicInfo is null in cif",hCif.getCrystallographicInfo());

			PDBCrystallographicInfo ciPdb = hPdb.getCrystallographicInfo();
			PDBCrystallographicInfo ciCif = hCif.getCrystallographicInfo();

			assertNotNull("space group null in pdb", ciPdb.getSpaceGroup());
			assertNotNull("space group null in cif", ciCif.getSpaceGroup());
			assertNotNull("crystal cell null in pdb",ciPdb.getCrystalCell());
			assertNotNull("crystal cell null in cif",ciCif.getCrystalCell());
			assertEquals("failed for space group short symbol pdb vs cif",
					ciPdb.getSpaceGroup().getShortSymbol(), ciCif.getSpaceGroup().getShortSymbol());

			CrystalCell ccPdb = ciPdb.getCrystalCell();
			CrystalCell ccCif = ciCif.getCrystalCell();

			assertEquals("failed for cell A:",ccPdb.getA(),ccCif.getA(),DELTA);
			assertEquals("failed for cell B:",ccPdb.getB(),ccCif.getB(),DELTA);
			assertEquals("failed for cell C:",ccPdb.getC(),ccCif.getC(),DELTA);
			assertEquals("failed for cell Alpha:",ccPdb.getAlpha(),ccCif.getAlpha(),DELTA);
			assertEquals("failed for cell Beta:",ccPdb.getBeta(),ccCif.getBeta(),DELTA);
			assertEquals("failed for cell Gamma:",ccPdb.getGamma(),ccCif.getGamma(),DELTA);


			if (ciPdb.getNcsOperators()==null) {
				assertTrue(ciCif.getNcsOperators()==null);
			} else {

				Matrix4d[] ncsOpersPdb = ciPdb.getNcsOperators();
				Matrix4d[] ncsOpersCif = ciCif.getNcsOperators();

				assertEquals("Number of NCS operators don't coincide", ncsOpersPdb.length, ncsOpersCif.length);

				for (int i=0;i<ncsOpersPdb.length;i++) {
					assertTrue("NCS operators "+i+" don't coincide",ncsOpersPdb[i].epsilonEquals(ncsOpersCif[i], 0.0001));
				}
			}
		}

		// biological assemblies
		// a) we don't test in non-crystallographic case because annotation is inconsistent between PDB and mmCIF,
		//    e.g. 2kli (NMR) has bioassembly annotation in mmCIF but not in PDB
		// b) we don't test virus entries (we check via looking at ncs operators==null):
		//    they are inconsistent PDB vs mmCIF (e.g. 1pgw has no oligomeric size in PDB, and 120 in mmCIF)
		if (isCrystallographic && hPdb.getCrystallographicInfo().getNcsOperators()==null
				// 1ruh, 2ms2, 2r06: virus proteins with data consistency issue: it's missing the MTRXn record (so it appears as ncs operators==null)
				&& (!sPdb.getPDBCode().equalsIgnoreCase("1ruh"))
				&& (!sPdb.getPDBCode().equalsIgnoreCase("2ms2"))
				&& (!sPdb.getPDBCode().equalsIgnoreCase("2r06"))) {

			assertEquals("Number of bioassemblies doesn't coincide",
					hPdb.getNrBioAssemblies(), hCif.getNrBioAssemblies());

			Map<Integer,BioAssemblyInfo> batPdb = hPdb.getBioAssemblies();
			Map<Integer,BioAssemblyInfo> batCif = hCif.getBioAssemblies();

			assertEquals("Size of bioassemblies map doesn't coincide with nr of bioassemblies",
					hPdb.getNrBioAssemblies(),batPdb.size());
			assertEquals("Size of bioassemblies maps don't coincide",batPdb.size(), batCif.size());

			for (int id:batPdb.keySet()) {
				assertTrue("Bioassembly id is not contained in mmCIF",batCif.containsKey(id));
				// there's an inconsistency in 4amh pdb vs mmCIF in mmSize
				if (sPdb.getPDBCode().equalsIgnoreCase("4amh")) continue;

				assertEquals("Macromolecular size of assembly "+id+" doesn't coincide",
						batPdb.get(id).getMacromolecularSize(), batCif.get(id).getMacromolecularSize());
			}
		}
	}

	private void testChains(Structure sPdb, Structure sCif) throws StructureException {
		assertNotNull(sPdb.getChains());
		assertNotNull(sCif.getChains());
		
		// sugar chains are badly annotated and inconsistent between pdb/mmcif
		// let's skip this test if we have sugar entities 

		if (!containsSugar(sCif)) {

			assertEquals(sPdb.getPolyChains().size(), sCif.getPolyChains().size());

			// some entries like 3c5e are inconsistent in residue numbering for UNL (unknown) residues between pdb and mmcif
			// skipping this test for them
			if (!containsUNL(sCif)) {
				assertEquals(sPdb.getNonPolyChains().size(), sCif.getNonPolyChains().size());
			}

			assertEquals(sPdb.getWaterChains().size(), sCif.getWaterChains().size());
			
			if (!containsUNL(sCif)) {
				assertEquals(sPdb.getChains().size(),sCif.getChains().size());
			}

		}

		

		
		Set<String> chainIds = new TreeSet<String>();
		for (Chain chain:sPdb.getPolyChains()){
			chainIds.add(chain.getName());
		}

		for (String chainId:chainIds) {
			testSingleChain(sPdb.getPolyChainByPDB(chainId), sCif.getPolyChainByPDB(chainId));
		}
	}

	private void testSingleChain(Chain cPdb, Chain cCif) {
		assertNotNull(cPdb);
		assertNotNull(cCif);

		String chainId = cPdb.getName();

		assertEquals("failed for getName():",cPdb.getName(),cCif.getName());
		// TODO no internalChainID if parsed from PDB, should an ID be assigned following the same rules as in mmCIF?
		//assertEquals("failed for getInternalChainID():",cPdb.getInternalChainID(),cCif.getInternalChainID());
		assertNotNull("getId is null",cCif.getId());
		assertTrue("id used in mmCIF files must be at most 4 characters",cCif.getId().length()<=4);
		assertEquals("chainID must be 1 character only, failed for pdb", 1, cPdb.getName().length());
		assertEquals("chainID must be 1 character only, failed for cif", 1, cCif.getName().length());

		// getCompound() is some times null for badly formatted PDB files (e.g. 4a10, all waters are in a separate chain F)
		if (isPolymer(cPdb)) {
			assertNotNull("getCompound is null in pdb (chain "+chainId+")",cPdb.getEntityInfo());
			assertNotNull("getCompound is null in cif (chain "+chainId+")",cCif.getEntityInfo());

			// for some badly formatted entries there are mismatches of mol_ids on pdb cs mmcif, e.g. 2efw
			// we thus count them and only warn at the end
			int molIdPdb = cPdb.getEntityInfo().getMolId();
			int molIdCif = cCif.getEntityInfo().getMolId();
			if (molIdPdb!=molIdCif) {
				logger.warn("Mismatching mol_id (entity_id) for {}. pdb: {}, mmCIF: {}",pdbId,molIdPdb,molIdCif);
				pdbIdsWithMismatchingMolIds.add(pdbId);
			}
		}


		assertNotNull("getParent is null in pdb (chain "+chainId+")",cPdb.getStructure());
		assertNotNull("getParent is null in cif (chain "+chainId+")",cCif.getStructure());


		assertEquals("failed for getAtomLength (chain "+chainId+"):",cPdb.getAtomLength(),cCif.getAtomLength());

		// entries with polymers composed of all unknowns (giving only-X sequences) can't be aligned seqres-to-atom (for PDB files)
		// we've got to skip them because they won't have seqres groups
		// e.g. is 1jnv chain A

		if (cPdb.getAtomSequence().matches("^X+$")) return;

		// note for getSeqResLength to work one needs the setAlignSeqRes option in the parsers

		assertEquals("failed for getSeqResLength pdb vs cif (chain "+chainId+"):",
						cPdb.getSeqResLength(),cCif.getSeqResLength());
		assertEquals("failed for getSeqResGroups().size pdb vs cif",
				cPdb.getSeqResGroups().size(), cCif.getSeqResGroups().size());
		assertEquals("getSeqResLength and getSeqResGroups.size should coincide in pdb:",
				cPdb.getSeqResLength(),cPdb.getSeqResGroups().size());
		assertEquals("getSeqResLength and getSeqResGroups.size should coincide in cif:",
				cCif.getSeqResLength(),cCif.getSeqResGroups().size());


		assertEquals("failed for getAtomLength:",cPdb.getAtomLength(),cCif.getAtomLength());
		assertEquals("failed for getAtomGroups().size pdb vs cif",
				cPdb.getAtomGroups().size(), cCif.getAtomGroups().size());
		assertEquals("getAtomLength and getAtomGroups.size should coincide in pdb:",
				cPdb.getAtomLength(),cPdb.getAtomGroups().size());
		assertEquals("getAtomLength and getAtomGroups.size should coincide in cif:",
				cCif.getAtomLength(),cCif.getAtomGroups().size());

		assertEquals("failed for getAtomGroups(GroupType.AMINOACID) pdb vs cif:",
				cPdb.getAtomGroups(GroupType.AMINOACID).size(),cCif.getAtomGroups(GroupType.AMINOACID).size());
		assertEquals("failed for getAtomGroups(GroupType.HETATM) pdb vs cif:",
				cPdb.getAtomGroups(GroupType.HETATM).size(),cCif.getAtomGroups(GroupType.HETATM).size());
		assertEquals("failed for getAtomGroups(GroupType.NUCLEOTIDE) pdb vs cif:",
				cPdb.getAtomGroups(GroupType.NUCLEOTIDE).size(),cCif.getAtomGroups(GroupType.NUCLEOTIDE).size());

		// In 4imj, chain F there's an  alignment ambiguity because of a repeat, so the seqres to atom alignment
		// doesn't work properly for it, we skip the rest of the test for this chain
		if (cPdb.getStructure().getPDBCode().equals("4IMJ") && cPdb.getName().equals("F")) return;

		assertEquals("failed for getSeqResGroups(GroupType.AMINOACID) pdb vs cif:",
				cPdb.getSeqResGroups(GroupType.AMINOACID).size(),cCif.getSeqResGroups(GroupType.AMINOACID).size());
		assertEquals("failed for getAtomGroups(GroupType.HETATM) pdb vs cif:",
				cPdb.getSeqResGroups(GroupType.HETATM).size(),cCif.getSeqResGroups(GroupType.HETATM).size());
		assertEquals("failed for getAtomGroups(GroupType.NUCLEOTIDE) pdb vs cif:",
				cPdb.getSeqResGroups(GroupType.NUCLEOTIDE).size(),cCif.getSeqResGroups(GroupType.NUCLEOTIDE).size());



		assertTrue("getAtomLength must be at least 1 in length (chain "+chainId+")",cPdb.getAtomLength()>=1);

		if (isPolymer(cPdb)) {
			// some badly formatted PDB files (e.g. 4a10, all waters are in a separate chain F) have 0 seqres length for some chains
			assertTrue("getSeqResLength must be at least 1 in length (chain "+chainId+")",cPdb.getSeqResLength()>=1);
		}

		// in the current implementation this is not a valid test, entries that have aminoacid residues in
		// ligands, e.g. 3o6g won't pass this test
		//assertTrue("getSeqResLength ("+cPdb.getSeqResLength()+") must be >= than getAtomGroups(GroupType.AMINOACID).size() ("+
		//		cPdb.getAtomGroups(GroupType.AMINOACID).size()+") (chain "+chainName+")",
		//		cPdb.getSeqResLength()>=cPdb.getAtomGroups(GroupType.AMINOACID).size());

		int allAtomGroupsSizePdb =
				cPdb.getAtomGroups(GroupType.AMINOACID).size()+
				cPdb.getAtomGroups(GroupType.HETATM).size()+
				cPdb.getAtomGroups(GroupType.NUCLEOTIDE).size();
		int allAtomGroupsSizeCif =
				cCif.getAtomGroups(GroupType.AMINOACID).size()+
				cCif.getAtomGroups(GroupType.HETATM).size()+
				cCif.getAtomGroups(GroupType.NUCLEOTIDE).size();

		assertEquals("failed for sum of all atom group sizes (hetatm+nucleotide+aminoacid) pdb vs mmcif",allAtomGroupsSizePdb,allAtomGroupsSizeCif);

		assertEquals("failed for getAtomLength==hetatm+aminos+nucleotide",cPdb.getAtomLength(), allAtomGroupsSizePdb);

		int allSeqResGroupsSizePdb =
				cPdb.getSeqResGroups(GroupType.AMINOACID).size()+
				cPdb.getSeqResGroups(GroupType.HETATM).size()+
				cPdb.getSeqResGroups(GroupType.NUCLEOTIDE).size();
		int allSeqResGroupsSizeCif =
				cCif.getSeqResGroups(GroupType.AMINOACID).size()+
				cCif.getSeqResGroups(GroupType.HETATM).size()+
				cCif.getSeqResGroups(GroupType.NUCLEOTIDE).size();

		assertEquals("failed for sum of all seqres group sizes (hetatm+nucleotide+aminoacid) pdb vs mmcif",allSeqResGroupsSizePdb,allSeqResGroupsSizeCif);

		assertEquals("failed for getSeqResLength==hetatm+aminos+nucleotide",cPdb.getSeqResLength(), allSeqResGroupsSizePdb);

	}


	private Structure getPdbStructure(String pdbId) throws IOException, StructureException {
		cache.setUseMmCif(false);
		// set parsing params here:
		params.setAlignSeqRes(true);
		//params.setLoadChemCompInfo(true);
		params.setParseBioAssembly(true);

		return cache.getStructure(pdbId);

	}

	private Structure getCifStructure(String pdbId) throws IOException, StructureException {
		cache.setUseMmCif(true);
		// set parsing params here:
		params.setAlignSeqRes(true);
		//params.setLoadChemCompInfo(true);
		params.setParseBioAssembly(true);

		return cache.getStructure(pdbId);

	}

	/**
	 * Reads a file containing a list of PDB codes.
	 * Lines starting with "#" will be treated as comments
	 * Will stop reading after finding an empty line, this is useful to quickly test a modified list.
	 * @param testSetFile
	 * @return
	 * @throws IOException
	 */
	private List<String> readTestSetFile(String testSetFile) throws IOException {

		InputStream inStream = this.getClass().getResourceAsStream(testSetFile);
		BufferedReader br = new BufferedReader(new InputStreamReader(inStream));

		List<String> list = new ArrayList<String>();

		String line;
		while ((line=br.readLine())!=null) {
			if (line.startsWith("#")) continue;
			if (line.isEmpty()) break;

			if (!line.matches("\\d\\w\\w\\w"))
				throw new IllegalArgumentException("The input test set "+testSetFile+" contains an invalid PDB code: "+line);

			list.add(line);
		}
		br.close();

		return list;
	}

	private boolean isPolymer(Chain chain) {

		for (Group group : chain.getSeqResGroups()) {
			if ((group instanceof AminoAcid) || (group instanceof NucleotideImpl)) {
				return true;
			}
		}

		// not a single amino-acid or nucleotide, must be something not polymeric
		return false;
	}
	
	private boolean containsSugar(Structure s) {
		for (EntityInfo e:s.getEntityInfos()) {
			if (e.getDescription().contains("SUGAR")) return true;
		}
		return false;
	}
	
	private boolean containsUNL(Structure s) {
		for (Chain c:s.getNonPolyChains()) {
			for (Group g:c.getAtomGroups()) {
				if (g.getPDBName().equals("UNL")) return true;
			}
		}
		return false;
	}
}
