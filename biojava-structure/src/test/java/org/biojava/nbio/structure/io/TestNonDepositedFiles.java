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
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.zip.GZIPInputStream;

import org.biojava.nbio.structure.Chain;
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
 * @author duarte_j
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
		
		assertEquals(2,s.getChains().size());
		
		// checking that heuristics in CompoundFinder work. We should have a single entity (compound)
		assertEquals(1, s.getCompounds().size());
		
		//System.out.println("Chains from incomplete header file: ");
		//checkChains(s);
		
		
		
		// trying without seqAlignSeqRes
		params.setAlignSeqRes(false);
		inStream = new GZIPInputStream(this.getClass().getResourceAsStream("/org/biojava/nbio/structure/io/1b8g_raw.pdb.gz"));
		s = pdbpars.parsePDBFile(inStream);
		assertNotNull(s);
		
		assertEquals(2,s.getChains().size());
				
		assertEquals(1, s.getCompounds().size());
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
		assertTrue(s.nrModels()==1);
		assertNull(s.getPDBHeader().getExperimentalTechniques());
	}
	
	private void checkChains(Structure s) {
		for (Chain chain:s.getChains()) {
			int seqResLength = chain.getSeqResLength();
			int atomLength = chain.getAtomLength();
			System.out.println("chain "+chain.getChainID()+", atomLength: "+atomLength+", seqResLength: "+seqResLength);
			//assertTrue("atom length ("+atomLength+") should be smaller than seqResLength ("+seqResLength+")",atomLength<=seqResLength);
			System.out.println("seq res groups size: "+chain.getSeqResGroups().size());
			System.out.println("num hetatom groups: "+chain.getAtomLigands().size());
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
		assertEquals(4, s.getCompounds().size());

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
		assertEquals(4, s.getCompounds().size());
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
		
		assertEquals(1, s.getCompounds().size());
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
		
		// 2 polymer chains with ligands, 1 purely water chain
		assertEquals(3, s.getChains().size());
		
		// 1 polymer entity, 1 water entity
		assertEquals(2, s.getCompounds().size());
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
	public void testNewLigandChain() throws IOException, StructureException {
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
		Chain c1 = s1.getChainByPDB("B");
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
		Chain c2 = s2.getChainByPDB("B");
		assertNotNull(c2);
		assertEquals(expectedNumLigands, c2.getAtomGroups().size());
		
		// pdb and mmcif should have same number of chains
		assertEquals(s1.getChains().size(), s2.getChains().size());
	}
	
	@Test
	public void testWaterOnlyChain() throws IOException, StructureException {
		
		// following 2 files are cut-down versions of 4a10
		InputStream pdbStream = new GZIPInputStream(this.getClass().getResourceAsStream("/org/biojava/nbio/structure/io/4a10_short.pdb.gz"));
		InputStream cifStream = new GZIPInputStream(this.getClass().getResourceAsStream("/org/biojava/nbio/structure/io/4a10_short.cif.gz"));
		
		PDBFileParser pdbpars = new PDBFileParser();
		Structure s1 = pdbpars.parsePDBFile(pdbStream) ;
		
		assertEquals(2, s1.getChains().size());
		
		Chain c1 = null;
		try {
			c1 = s1.getChainByPDB("F");
			
		} catch (StructureException e) {
			fail("Got StructureException while looking for water-only chain F");
		}
		
		// checking that compounds are linked
		assertNotNull(c1.getCompound());
		
		// checking that the water molecule was assigned an ad-hoc compound
		assertEquals(2,s1.getCompounds().size());
		
		
		
		MMcifParser mmcifpars = new SimpleMMcifParser();
		SimpleMMcifConsumer consumer = new SimpleMMcifConsumer();
		mmcifpars.addMMcifConsumer(consumer);
		mmcifpars.parse(cifStream) ;
		Structure s2 = consumer.getStructure();
		
		
		assertEquals(2, s2.getChains().size());
		
		Chain c = null;
		try {
			c = s2.getChainByPDB("F");
			
		} catch (StructureException e) {
			fail("Got StructureException while looking for water-only chain F");
		}
		
		// checking that compounds are linked
		assertNotNull(c.getCompound());
		
		// checking that the water molecule was assigned an ad-hoc compound
		assertEquals(2,s2.getCompounds().size());
	}
}
