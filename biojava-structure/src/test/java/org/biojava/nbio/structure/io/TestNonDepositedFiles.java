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
package org.biojava.nbio.structure.io;

import static org.junit.Assert.*;

import java.io.BufferedReader;
import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.List;
import java.util.zip.GZIPInputStream;

import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.EntityInfo;
import org.biojava.nbio.structure.EntityType;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.io.mmcif.MMcifParser;
import org.biojava.nbio.structure.io.mmcif.SimpleMMcifConsumer;
import org.biojava.nbio.structure.io.mmcif.SimpleMMcifParser;
import org.biojava.nbio.structure.xtal.CrystalCell;
import org.junit.Test;

/**
 * Tests for non-deposited PDB/mmCIF files, i.e. any kind of "raw" file
 * lacking significant parts of the headers.
 *
 * Some things tested:
 * - heuristics to guess isNMR, isCrystallographic
 *
 * @author Jose Duarte
 *
 */
public class TestNonDepositedFiles {

	@Test
	public void test1B8GnoSeqresPdb() throws IOException, StructureException {

		InputStream inStream = new GZIPInputStream(this.getClass().getResourceAsStream("/org/biojava/nbio/structure/io/1b8g_raw.pdb.gz"));
		assertNotNull(inStream);

		PDBFileParser pdbpars = new PDBFileParser();
		FileParsingParameters params = new FileParsingParameters();
		params.setAlignSeqRes(true);
		pdbpars.setFileParsingParameters(params);

		Structure s = pdbpars.parsePDBFile(inStream) ;

		assertNotNull(s);

		assertTrue(s.isCrystallographic());
		assertFalse(s.isNmr());
		assertTrue(s.nrModels()==1);
		assertNull(s.getPDBHeader().getExperimentalTechniques());

		assertNotNull(s.getCrystallographicInfo().getCrystalCell());
		assertNotNull(s.getCrystallographicInfo().getSpaceGroup());

		assertEquals(s.getCrystallographicInfo().getSpaceGroup().getShortSymbol(),"P 1 21 1");

		CrystalCell cell = s.getCrystallographicInfo().getCrystalCell();
		assertTrue(cell.isCellReasonable());

		// TODO get the scale matrix from the PDB file and check it against the calculated one:
		//cell.checkScaleMatrixConsistency(scaleMatrix);
		//cell.checkScaleMatrix(scaleMatrix);

		// 2 protein chanis, 2 nonpoly PLP chains, 2 water chains
		assertEquals(6,s.getChains().size());

		// checking that heuristics in CompoundFinder work. We should have 1 polymer entity (protein) + 1 nonpoly entity (PLP) + 1 water entity
		assertEquals(3, s.getEntityInfos().size());
		assertEquals(EntityType.POLYMER, s.getEntityById(1).getType());

		//System.out.println("Chains from incomplete header file: ");
		//checkChains(s);



		// trying without seqAlignSeqRes
		params.setAlignSeqRes(false);
		inStream = new GZIPInputStream(this.getClass().getResourceAsStream("/org/biojava/nbio/structure/io/1b8g_raw.pdb.gz"));
		s = pdbpars.parsePDBFile(inStream);
		assertNotNull(s);

		assertEquals(6,s.getChains().size());

		assertEquals(3, s.getEntityInfos().size());
		assertEquals(EntityType.POLYMER, s.getEntityById(1).getType());
	}

	//@Test
	public void test1B8G() throws IOException, StructureException {

		AtomCache cache = new AtomCache();

		StructureIO.setAtomCache(cache);

		cache.setUseMmCif(true);
		Structure s = StructureIO.getStructure("1B8G");

		System.out.println("Chains from full deposited file: ");
		checkChains(s);
	}

	@Test
	public void test3C5F() throws IOException, StructureException {
		InputStream inStream = new GZIPInputStream(this.getClass().getResourceAsStream("/org/biojava/nbio/structure/io/3c5f_raw.pdb.gz"));
		assertNotNull(inStream);

		PDBFileParser pdbpars = new PDBFileParser();
		FileParsingParameters params = new FileParsingParameters();
		params.setAlignSeqRes(true);
		pdbpars.setFileParsingParameters(params);

		Structure s = pdbpars.parsePDBFile(inStream) ;

		// multi-model X-ray diffraction entry, thus:
		assertFalse(s.isNmr());
		assertTrue(s.isCrystallographic());
		assertTrue(s.nrModels()>1);
		assertNull(s.getPDBHeader().getExperimentalTechniques());

	}

	@Test
	public void test4B19() throws IOException, StructureException {
		InputStream inStream = new GZIPInputStream(this.getClass().getResourceAsStream("/org/biojava/nbio/structure/io/4b19_raw.pdb.gz"));
		assertNotNull(inStream);

		PDBFileParser pdbpars = new PDBFileParser();
		FileParsingParameters params = new FileParsingParameters();
		params.setAlignSeqRes(true);
		pdbpars.setFileParsingParameters(params);

		Structure s = pdbpars.parsePDBFile(inStream) ;

		// multi-model NMR entry, thus:
		assertTrue(s.isNmr());
		assertFalse(s.isCrystallographic());
		assertTrue(s.nrModels()>1);
		assertNull(s.getPDBHeader().getExperimentalTechniques());

	}

	@Test
	public void test2M7Y() throws IOException {
		InputStream inStream = new GZIPInputStream(this.getClass().getResourceAsStream("/org/biojava/nbio/structure/io/2m7y_raw.pdb.gz"));
		assertNotNull(inStream);

		PDBFileParser pdbpars = new PDBFileParser();
		FileParsingParameters params = new FileParsingParameters();
		params.setAlignSeqRes(true);
		pdbpars.setFileParsingParameters(params);

		Structure s = pdbpars.parsePDBFile(inStream) ;

		// single-model NMR entry, thus:
		//assertTrue(s.isNmr());  // we can't detect it properly, because it's single model!
		assertFalse(s.isCrystallographic()); // at least this we can detect from the unreasonable crystal cell
		assertEquals(1, s.nrModels());
		assertNull(s.getPDBHeader().getExperimentalTechniques());

		// testing that on single chain pdb files we assign an entity type, issue #767
		assertEquals(EntityType.POLYMER, s.getEntityById(1).getType());
	}

	private void checkChains(Structure s) {
		for (Chain chain:s.getChains()) {
			int seqResLength = chain.getSeqResLength();
			int atomLength = chain.getAtomLength();
			System.out.println("chain "+chain.getId()+", atomLength: "+atomLength+", seqResLength: "+seqResLength);
			//assertTrue("atom length ("+atomLength+") should be smaller than seqResLength ("+seqResLength+")",atomLength<=seqResLength);
			System.out.println("seq res groups size: "+chain.getSeqResGroups().size());
		}
	}

	/**
	 * A test for reading a phenix-produced (ver 1.9_1692) mmCIF file.
	 * This is the file submitted to the PDB for deposition of entry 4lup
	 * See github issue #234
	 * @throws IOException
	 */
	@Test
	public void testPhenixCifFile() throws IOException {
		InputStream inStream = new GZIPInputStream(this.getClass().getResourceAsStream("/org/biojava/nbio/structure/io/4lup_phenix_output.cif.gz"));
		MMcifParser parser = new SimpleMMcifParser();

		SimpleMMcifConsumer consumer = new SimpleMMcifConsumer();

		FileParsingParameters fileParsingParams = new FileParsingParameters();
		fileParsingParams.setAlignSeqRes(true);

		consumer.setFileParsingParameters(fileParsingParams);

		parser.addMMcifConsumer(consumer);

		parser.parse(new BufferedReader(new InputStreamReader(inStream)));

		Structure s = consumer.getStructure();

		assertNotNull(s);

		assertTrue(s.isCrystallographic());

		// all ligands are into their own chains, so we have 2 proteins, 2 nucleotide chains, 1 ligand chain and 1 purely water chain
		assertEquals(6, s.getChains().size());

		// 4 entities: 1 protein, 1 nucleotide, 1 water, 1 ligand (EDO)
		assertEquals(4, s.getEntityInfos().size());
		int[] counts = countEntityTypes(s.getEntityInfos());
		assertEquals(2, counts[0]);
		assertEquals(1, counts[1]);
		assertEquals(1, counts[2]);


	}

	@Test
	public void testPhenixPdbFile() throws IOException {
		InputStream inStream = new GZIPInputStream(this.getClass().getResourceAsStream("/org/biojava/nbio/structure/io/4lup_phenix_output.pdb.gz"));

		PDBFileParser pdbpars = new PDBFileParser();
		FileParsingParameters params = new FileParsingParameters();
		params.setAlignSeqRes(true);
		pdbpars.setFileParsingParameters(params);

		Structure s = pdbpars.parsePDBFile(inStream) ;

		assertNotNull(s);

		assertTrue(s.isCrystallographic());

		// all ligands are into their own chains, so we have 2 proteins, 2 nucleotide chains, 1 ligand chain and 1 purely water chain
		assertEquals(6, s.getChains().size());

		// 4 entities: 1 protein, 1 nucleotide, 1 water, 1 ligand (EDO)
		assertEquals(4, s.getEntityInfos().size());
		int[] counts = countEntityTypes(s.getEntityInfos());
		assertEquals(2, counts[0]);
		assertEquals(1, counts[1]);
		assertEquals(1, counts[2]);

	}

	@Test
	public void testPhaserPdbFile() throws IOException {
		InputStream inStream = new GZIPInputStream(this.getClass().getResourceAsStream("/org/biojava/nbio/structure/io/4lup_phaser_output.pdb.gz"));

		PDBFileParser pdbpars = new PDBFileParser();
		FileParsingParameters params = new FileParsingParameters();
		params.setAlignSeqRes(true);
		pdbpars.setFileParsingParameters(params);

		Structure s = pdbpars.parsePDBFile(inStream) ;

		assertNotNull(s);

		assertTrue(s.isCrystallographic());

		assertEquals(2, s.getChains().size());

		assertEquals(1, s.getEntityInfos().size());
	}


	@Test
	public void testRefmacPdbFile() throws IOException {
		InputStream inStream = new GZIPInputStream(this.getClass().getResourceAsStream("/org/biojava/nbio/structure/io/rnase_refmac_output.pdb.gz"));

		PDBFileParser pdbpars = new PDBFileParser();
		FileParsingParameters params = new FileParsingParameters();
		params.setAlignSeqRes(true);
		pdbpars.setFileParsingParameters(params);

		Structure s = pdbpars.parsePDBFile(inStream) ;

		assertNotNull(s);

		assertTrue(s.isCrystallographic());

		// 2 polymer chains with 1 ligand per chain, 1 purely water chain = 5 chains
		assertEquals(5, s.getChains().size());

		// 1 polymer entity, 1 nonpoly entity, 1 water entity
		assertEquals(3, s.getEntityInfos().size());
		int[] counts = countEntityTypes(s.getEntityInfos());
		assertEquals(1, counts[0]);
		assertEquals(1, counts[1]);
		assertEquals(1, counts[2]);

	}

	/**
	 * This test represents a common situation for a non-deposited structure.
	 * When building with common crystallography software, the user often adds new
	 * ligands (or solvent) molecules as new chains.  Only prior to deposition
	 * then relabel them so that they belong to the same chain as the polymeric residues.
	 *
	 * In this case, the ligands represent valuable information and should not be discarded.
	 */
	@Test
	public void testNewLigandChain() throws IOException {
		// Test the file parsing speed when the files are already downloaded.

		InputStream pdbStream = new GZIPInputStream(this.getClass().getResourceAsStream("/ligandTest.pdb.gz"));
		InputStream cifStream = new GZIPInputStream(this.getClass().getResourceAsStream("/ligandTest.cif.gz"));

		assertNotNull(cifStream);
		assertNotNull(pdbStream);

		FileParsingParameters params = new FileParsingParameters();
		PDBFileParser pdbpars = new PDBFileParser();
		pdbpars.setFileParsingParameters(params);
		Structure s1 = pdbpars.parsePDBFile(pdbStream) ;

		// The chain B should be present with 1 ligand HEM
		Chain c1 = s1.getNonPolyChainsByPDB("B").get(0);
		assertNotNull(c1);

		int expectedNumLigands = 1;
		assertEquals(expectedNumLigands, c1.getAtomGroups().size());

		MMcifParser mmcifpars = new SimpleMMcifParser();
		SimpleMMcifConsumer consumer = new SimpleMMcifConsumer();
		consumer.setFileParsingParameters(params);
		mmcifpars.addMMcifConsumer(consumer);
		mmcifpars.parse(cifStream) ;
		Structure s2 = consumer.getStructure();

		// The chain B should be present with 1 ligand HEM
		Chain c2 = s2.getNonPolyChainsByPDB("B").get(0);
		assertNotNull(c2);
		assertEquals(expectedNumLigands, c2.getAtomGroups().size());

		// pdb and mmcif should have same number of chains
		assertEquals(s1.getChains().size(), s2.getChains().size());
	}
	
	@Test
	public void testWaterOnlyChainPdb() throws IOException {

		// following file is cut-down version of 4a10
		InputStream pdbStream = new GZIPInputStream(this.getClass().getResourceAsStream("/org/biojava/nbio/structure/io/4a10_short.pdb.gz"));

		PDBFileParser pdbpars = new PDBFileParser();
		Structure s1 = pdbpars.parsePDBFile(pdbStream) ;

		assertEquals(2, s1.getChains().size());

		Chain c1 = s1.getWaterChainByPDB("F");

		assertNotNull("Got null when looking for water-only chain with author id F", c1);

		// checking that compounds are linked
		assertNotNull(c1.getEntityInfo());

		// checking that the water molecule was assigned an ad-hoc compound
		assertEquals(2,s1.getEntityInfos().size());

	}
	
	@Test
	public void testWaterOnlyChainCif() throws IOException {

		// following file is cut-down versions of 4a10
		InputStream cifStream = new GZIPInputStream(this.getClass().getResourceAsStream("/org/biojava/nbio/structure/io/4a10_short.cif.gz"));

		MMcifParser mmcifpars = new SimpleMMcifParser();
		SimpleMMcifConsumer consumer = new SimpleMMcifConsumer();
		mmcifpars.addMMcifConsumer(consumer);
		mmcifpars.parse(cifStream) ;
		Structure s2 = consumer.getStructure();


		assertEquals(2, s2.getChains().size());

		Chain c = s2.getWaterChainByPDB("F");

		assertNotNull("Got null when looking for water-only chain with author id F", c);

		// checking that compounds are linked
		assertNotNull(c.getEntityInfo());

		// checking that the water molecule was assigned an ad-hoc compound
		assertEquals(2,s2.getEntityInfos().size());
		
		Chain cAsymId = s2.getWaterChain("E");
		assertNotNull("Got null when looking for water-only chain with asym id E", cAsymId);
		assertSame(c, cAsymId);
		
	}
	
	/**
	 * Some PDB files coming from phenix or other software can have a CRYST1 line without z and not padded with white-spaces 
	 * for the space group column.
	 * @throws IOException
	 * @since 5.0.0
	 */
	@Test
	public void testCryst1Parsing() throws IOException {
		String cryst1Line = "CRYST1   11.111   11.111  111.111  70.00  80.00  60.00 P 1";
		Structure s;
		PDBFileParser pdbPars = new PDBFileParser();
		try(InputStream is = new ByteArrayInputStream(cryst1Line.getBytes()) ) {
			s = pdbPars.parsePDBFile(is);
		}
		assertEquals("P 1", s.getPDBHeader().getCrystallographicInfo().getSpaceGroup().getShortSymbol());
	}

	private static int[] countEntityTypes(List<EntityInfo> entities) {
		int countPoly = 0;
		int countNonPoly = 0;
		int countWater = 0;
		for (EntityInfo e:entities) {
			if (e.getType()==EntityType.POLYMER) countPoly++;
			if (e.getType()==EntityType.NONPOLYMER) countNonPoly++;
			if (e.getType()==EntityType.WATER) countWater++;
		}
		int[] counts = {countPoly, countNonPoly, countWater}; 
		return counts;
		
	}
}
