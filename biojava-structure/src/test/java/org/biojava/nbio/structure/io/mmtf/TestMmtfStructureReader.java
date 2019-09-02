package org.biojava.nbio.structure.io.mmtf;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.nio.file.Paths;
import java.util.List;

import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.io.mmcif.ChemCompGroupFactory;
import org.biojava.nbio.structure.io.mmcif.DownloadChemCompProvider;
import org.junit.Ignore;
import org.junit.Test;
import static org.junit.Assert.*;

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
	@Test
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

	/**
	 * Test for issue https://github.com/biojava/biojava/issues/792
	 */
	@Test
	@Ignore("Issue not fixed yet")
	public void checkNonStandardAminoSeqresGroupsPopulated() throws StructureException, IOException {
	    // 2X3T, see issue https://github.com/biojava/biojava/issues/792
        // Load a structure in mmtf format
        AtomCache cache = new AtomCache();
        FileParsingParameters params = new FileParsingParameters();
        cache.setFileParsingParams(params);
        cache.setUseMmCif(false);
        cache.setUseMmtf(true);

        StructureIO.setAtomCache(cache);

        ChemCompGroupFactory.setChemCompProvider(new DownloadChemCompProvider());

        Structure structure1 = StructureIO.getStructure("2X3T");
        // chain E is a glycopeptide with unobserved non-standard aminoacids. Because of mmtf limitations (representing seqres sequences as 1-letter strings) the non-standard unobserved residues are read as null
        List<Group> seqresGroups = structure1.getChain("E").getSeqResGroups();
        for (Group g : seqresGroups) {
            assertNotNull("SeqRes group should not be null", g);
        }

    }

}
