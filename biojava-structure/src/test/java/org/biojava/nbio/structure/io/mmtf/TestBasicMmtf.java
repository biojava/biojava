package org.biojava.nbio.structure.io.mmtf;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;
import java.nio.file.Paths;
import java.util.ArrayList;

import org.biojava.nbio.structure.AminoAcidImpl;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.AtomImpl;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.ChainImpl;
import org.biojava.nbio.structure.Element;
import org.biojava.nbio.structure.EntityInfo;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.PDBHeader;
import org.biojava.nbio.structure.ResidueNumber;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureImpl;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

/**
 * Test that Biojava can read and write MMTF data.
 * @author Anthony Bradley
 *
 */
public class TestBasicMmtf {

    /**
     * A test folder for testing writing files.
     */
    @Rule
    public TemporaryFolder testFolder = new TemporaryFolder();
	
	
	/**
	 * Test that Biojava can read a file from the file system.
	 * @throws IOException
	 */
	@Test
	public void testRead() throws IOException {
		ClassLoader classLoader = getClass().getClassLoader();
		Structure structure = MmtfActions.readFromFile((Paths.get(classLoader.getResource("org/biojava/nbio/structure/io/mmtf/4CUP.mmtf").getPath())));
		assertEquals(structure.getPDBCode(),"4CUP");
		assertEquals(structure.getChains().size(),6);
	}
	
	/**
	 * Test the writing of Structure objects to a file.
	 * @throws IOException 
	 */
	@Test
	public void testWrite() throws IOException {
		Structure structure = new StructureImpl();
		PDBHeader pdbHeader = new PDBHeader();
		pdbHeader.setExperimentalTechnique("X-RAY DIFFRACTION");
		structure.setEntityInfos(new ArrayList<EntityInfo>());
		structure.setPDBHeader(pdbHeader);
		Chain chain = new ChainImpl();
		Group group = new AminoAcidImpl(); 
		group.setPDBName("FKF");
		Atom atom = new AtomImpl();
		atom.setName("A");
		atom.setElement(Element.Ag);
		atom.setCoords(new double[] {1.0,2.0,3.0});
		chain.addGroup(group);
		group.addAtom(atom);
		ResidueNumber residueNumber = new ResidueNumber();
		residueNumber.setInsCode('A');
		residueNumber.setSeqNum(100);
		group.setResidueNumber(residueNumber);
		structure.addChain(chain);
		File tempFile = testFolder.newFile("tmpfile");
		MmtfActions.writeToFile(structure, tempFile.toPath());		
	}
}
