package org.biojava.nbio.structure.io.mmtf;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.nio.file.Paths;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.junit.Ignore;
import org.junit.Test;

/**
 * Test the Biojava MMTF reader.
 * 
 * @author Anthony Bradley
 * @author Aleix Lafita
 *
 */
public class TestMmtfStructureReader {

	/**
	 * Test reading an MMTF file into a BioJava structure.
	 */
	@Test
	public void testRead() throws IOException {
		
		// Get the MMTF file from the resources folder
		ClassLoader classLoader = getClass().getClassLoader();
		String resource = "org/biojava/nbio/structure/io/mmtf/4CUP.mmtf";
		
		// Load the structure into memory
		Structure structure = MmtfActions.readFromFile((
				Paths.get(classLoader.getResource(resource).getPath())));
		
		// Check header properties of the structure
		assertEquals(structure.getPDBCode(), "4CUP");
		assertEquals(MmtfUtils.dateToIsoString(structure.getPDBHeader().getDepDate()), 
				"2014-03-21");
		
		assertEquals(structure.getChains().size(), 6);
	}
	
	/**
	 * Compare structures loaded from MMCIF and MMTF files.
	 */
	@Test @Ignore
	public void compareMmcif() throws IOException, StructureException {
		
		// Get the MMTF and MMCIF files from the resources folder
		ClassLoader classLoader = getClass().getClassLoader();
		String resource = "org/biojava/nbio/structure/io/mmtf/4CUP";
		
		// Load the structures into memory
		Structure mmtf = MmtfActions.readFromFile((
				Paths.get(classLoader.getResource(resource + ".mmtf").getPath())));
		Structure mmcif = StructureIO.getStructure(classLoader.getResource(resource + ".cif").getPath());
		
		// Compare the dates of the structure
		assertEquals(mmcif.getPDBHeader().getDepDate(), 
				mmtf.getPDBHeader().getDepDate());
		
		// Compare the experimental method
		assertEquals(mmcif.getPDBHeader().getExperimentalTechniques(), 
				mmtf.getPDBHeader().getExperimentalTechniques());
		
		// Compare the SEQRES, see issue https://github.com/biojava/biojava/issues/671
		assertEquals(mmcif.getChainByIndex(0).getSeqResSequence(), 
				mmtf.getChainByIndex(0).getSeqResSequence());
		
	}

}
