package org.biojava.nbio.structure.io.mmtf;

import java.io.File;
import java.io.IOException;
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
import org.biojava.nbio.structure.io.mmcif.model.ChemComp;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

/**
 * Test the Biojava MMTF writer.
 * 
 * @author Anthony Bradley
 * @author Aleix Lafita
 *
 */
public class TestMmtfStructureWriter {

	@Rule
	public TemporaryFolder testFolder = new TemporaryFolder();
	
	/**
	 * Test the writing of Structure objects to a file.
	 * 
	 * @throws IOException
	 */
	@Test
	public void testWrite() throws IOException {
		
		// Create a structure
		Structure structure = new StructureImpl();
		
		// Add some header information
		PDBHeader pdbHeader = new PDBHeader();
		pdbHeader.setExperimentalTechnique("X-RAY DIFFRACTION");
		structure.setPDBHeader(pdbHeader);
		
		// Create one chain
		structure.setEntityInfos(new ArrayList<EntityInfo>());
		Chain chain = new ChainImpl();
		chain.setId("A");
		chain.setName("A");
		
		// Create one group
		Group group = new AminoAcidImpl();
		group.setPDBName("FKF");
		ChemComp chemComp = new ChemComp();
		chemComp.setType("TYPfdl");
		chemComp.setOne_letter_code("A");
		group.setChemComp(chemComp);
		
		// Create one Atom
		Atom atom = new AtomImpl();
		atom.setName("A");
		atom.setElement(Element.Ag);
		atom.setCoords(new double[] { 1.0, 2.0, 3.0 });
		
		// Link together the objects
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
