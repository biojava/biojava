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

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.zip.GZIPInputStream;

import static org.junit.Assert.*;

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
		
		// all ligands are into their own chains, so we have 2 proteins, 2 nucleotide chains and 1 ligand chain
		assertEquals(5, s.getChains().size());
		
		assertEquals(2, s.getCompounds().size());

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
		
		// all ligands are into their own chains, so we have 2 proteins, 2 nucleotide chains and 1 ligand chain
		assertEquals(5, s.getChains().size());
		
		assertEquals(2, s.getCompounds().size());
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
		
		assertEquals(2, s.getChains().size());
		
		assertEquals(1, s.getCompounds().size());
	}

}
