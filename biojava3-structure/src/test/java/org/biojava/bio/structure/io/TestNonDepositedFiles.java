package org.biojava.bio.structure.io;

import static org.junit.Assert.*;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.zip.GZIPInputStream;

import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.io.mmcif.MMcifParser;
import org.biojava.bio.structure.io.mmcif.SimpleMMcifConsumer;
import org.biojava.bio.structure.io.mmcif.SimpleMMcifParser;
import org.biojava.bio.structure.xtal.CrystalCell;
import org.biojava3.structure.StructureIO;
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
		
		InputStream inStream = new GZIPInputStream(this.getClass().getResourceAsStream("/org/biojava/bio/structure/io/1b8g_raw.pdb.gz"));
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
		
		//System.out.println("Chains from incomplete header file: ");
		//checkChains(s);
		
		
		
		// trying without seqAlignSeqRes
		params.setAlignSeqRes(false);
		inStream = new GZIPInputStream(this.getClass().getResourceAsStream("/org/biojava/bio/structure/io/1b8g_raw.pdb.gz"));
		s = pdbpars.parsePDBFile(inStream);
		assertNotNull(s);
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
		InputStream inStream = new GZIPInputStream(this.getClass().getResourceAsStream("/org/biojava/bio/structure/io/3c5f_raw.pdb.gz"));
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
		InputStream inStream = new GZIPInputStream(this.getClass().getResourceAsStream("/org/biojava/bio/structure/io/4b19_raw.pdb.gz"));
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
		InputStream inStream = new GZIPInputStream(this.getClass().getResourceAsStream("/org/biojava/bio/structure/io/2m7y_raw.pdb.gz"));
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
	 * @throws IOException
	 */
	//@Test
	public void testPhenixFile() throws IOException {
		InputStream inStream = new GZIPInputStream(this.getClass().getResourceAsStream("/org/biojava/bio/structure/io/4lup_phenix_output.cif.gz"));
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
	}

}
