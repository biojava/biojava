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
package org.biojava.nbio.structure.test.symmetry;

import static org.junit.Assert.*;

import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.cluster.SubunitClustererMethod;
import org.biojava.nbio.structure.cluster.SubunitClustererParameters;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetryDetector;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetryParameters;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetryResults;
import org.junit.Test;

/**
 * Test the {@link QuatSymmetryDetector} algorithm for symmetry detection on a
 * variety of structures with different symmetries.
 * 
 * @author Peter Rose
 * @author Aleix Lafita
 *
 */
public class TestQuatSymmetryDetectorExamples {

	/**
	 * An NMR multi-model entry: 1B4C
	 * 
	 * @throws StructureException
	 * @throws IOException
	 */
	@Test
	public void testNMR() throws IOException, StructureException {

		// as of mmCIF v5 there's no bioassemblies for NMR entries, so now we use the AU (bioassembly 0) - JD 2017-08-02
		Structure pdb = StructureIO.getStructure("BIO:1b4c:0");

		SubunitClustererParameters clusterParams = new SubunitClustererParameters();
		QuatSymmetryParameters symmParams = new QuatSymmetryParameters();
		QuatSymmetryResults symmetry = QuatSymmetryDetector.calcGlobalSymmetry(
				pdb, symmParams, clusterParams);

		// C2 symmetry non pseudosymmetric
		assertEquals("C2", symmetry.getSymmetry());
		assertEquals("A2", symmetry.getStoichiometry());
		assertFalse(symmetry.isPseudoStoichiometric());

	}
	
	/** 
	 * Test a dihedral symmetry: 2VML 
	 */
	@Test
	public void testDihedral() throws IOException, StructureException {

		Structure pdb = StructureIO.getStructure("BIO:2vml:1");

		SubunitClustererParameters clusterParams = new SubunitClustererParameters();
		clusterParams.setClustererMethod(SubunitClustererMethod.SEQUENCE);
		QuatSymmetryParameters symmParams = new QuatSymmetryParameters();
		QuatSymmetryResults symmetry = QuatSymmetryDetector.calcGlobalSymmetry(
				pdb, symmParams, clusterParams);

		// D3 symmetry non pseudosymmetric
		assertEquals("D3", symmetry.getSymmetry());
		assertEquals("A6B6", symmetry.getStoichiometry());
		assertFalse(symmetry.isPseudoStoichiometric());


	}

	/**
	 * Hemoglobin has both symmetry and pseudosymmetry: 4HHB
	 * 
	 * @throws StructureException
	 * @throws IOException
	 */
	@Test
	public void testPseudosymmetry() throws IOException, StructureException {

		Structure pdb = StructureIO.getStructure("BIO:4hhb:1");

		// Non-pseudosymmetry
		SubunitClustererParameters clusterParams = new SubunitClustererParameters();
		clusterParams.setClustererMethod(SubunitClustererMethod.SEQUENCE);
		QuatSymmetryParameters symmParams = new QuatSymmetryParameters();
		QuatSymmetryResults symmetry = QuatSymmetryDetector.calcGlobalSymmetry(
				pdb, symmParams, clusterParams);

		// C2 symmetry
		assertEquals("C2", symmetry.getSymmetry());
		assertEquals("A2B2", symmetry.getStoichiometry());

		// Use pseudosymmetry (structural clustering)
		clusterParams.setClustererMethod(SubunitClustererMethod.STRUCTURE);
		symmetry = QuatSymmetryDetector.calcGlobalSymmetry(pdb, symmParams,
				clusterParams);

		// D2 pseudo-symmetry
		assertEquals("D2", symmetry.getSymmetry());
		assertEquals("A4", symmetry.getStoichiometry());
		assertTrue(symmetry.isPseudoStoichiometric());
	}

	/**
	 * A structure with no global symmetry, but local symmetry: 4P2C
	 * 
	 * @throws IOException
	 * @throws StructureException
	 */
	@Test
	public void testLocal() throws IOException, StructureException {

		Structure pdb = StructureIO.getStructure("BIO:4p2c:1");

		// Global Symmetry
		SubunitClustererParameters clusterParams = new SubunitClustererParameters();
		QuatSymmetryParameters symmParams = new QuatSymmetryParameters();
		QuatSymmetryResults symmetry = QuatSymmetryDetector.calcGlobalSymmetry(
				pdb, symmParams, clusterParams);

		// C1 global symmetry
		assertEquals("C1", symmetry.getSymmetry());
		assertEquals("A5B5C", symmetry.getStoichiometry());

		// Local symmetry
		List<QuatSymmetryResults> localSymm = QuatSymmetryDetector
				.calcLocalSymmetries(pdb, symmParams, clusterParams);

		// C5 local symmetry excluding chain A
		assertEquals(localSymm.size(), 3);
		assertEquals("C5", localSymm.get(0).getSymmetry());
		assertEquals("C5", localSymm.get(1).getSymmetry());
		assertEquals("C5", localSymm.get(2).getSymmetry());

		// Two A5 and one A5B5 stoichiometries as local symmetry
		List<String> stoich = localSymm.stream().map(t -> t.getStoichiometry())
				.collect(Collectors.toList());

		assertTrue(stoich.contains("A5"));
		assertTrue(stoich.contains("A5B5"));
	}

	/**
	 * A structure with combined internal and quaternary symmetry: 4E3E
	 * 
	 * @throws IOException
	 * @throws StructureException
	 */
	@Test
	public void testInternalSymmetry() throws IOException, StructureException {

		Structure pdb = StructureIO.getStructure("BIO:4e3e:1");

		// Internal symmetry analysis, use structural clustering
		SubunitClustererParameters cp = new SubunitClustererParameters();
		cp.setClustererMethod(SubunitClustererMethod.STRUCTURE);
		cp.setInternalSymmetry(true);
		cp.setCoverageThreshold(0.75); // Lower coverage for internal symm

		QuatSymmetryParameters symmParams = new QuatSymmetryParameters();
		QuatSymmetryResults symmetry = QuatSymmetryDetector.calcGlobalSymmetry(
				pdb, symmParams, cp);

		// D2 combined internal and quaternary symmetry
		assertEquals("D3", symmetry.getSymmetry());
		assertEquals("A6", symmetry.getStoichiometry());

	}

	/**
	 * A structure with helical symmetry: 1B47
	 * 
	 * @throws IOException
	 * @throws StructureException
	 */
	@Test
	public void testHelical() throws IOException, StructureException {

		Structure pdb = StructureIO.getStructure("BIO:1B47:1");

		SubunitClustererParameters cp = new SubunitClustererParameters();
		QuatSymmetryParameters symmParams = new QuatSymmetryParameters();
		QuatSymmetryResults symmetry = QuatSymmetryDetector.calcGlobalSymmetry(
				pdb, symmParams, cp);

		// H symmetry A3 stoichiometry
		assertEquals("H", symmetry.getSymmetry());
		assertEquals("A3", symmetry.getStoichiometry());

	}

	/**
	 * A structure with local helical symmetry: 5JLF
	 * 
	 * @throws IOException
	 * @throws StructureException
	 */
	@Test
	public void testHelicalLocal() throws IOException, StructureException {

		Structure pdb = StructureIO.getStructure("BIO:5JLF:1");

		SubunitClustererParameters cp = new SubunitClustererParameters();
		QuatSymmetryParameters symmParams = new QuatSymmetryParameters();
		List<QuatSymmetryResults> results = QuatSymmetryDetector
				.calcLocalSymmetries(pdb, symmParams, cp);

		// H symmetry A5 stoichiometry
		assertEquals("H", results.get(0).getSymmetry());
		assertEquals("A5", results.get(0).getStoichiometry());

	}
	
	/**
	 * A structure with very similar entities (clustering at 95% seq id): 4DZ8
	 * @throws IOException
	 * @throws StructureException
	 */
	@Test
	public void testPseudoIdentity95() throws IOException, StructureException {
		Structure pdb = StructureIO.getStructure("BIO:4DZ8:1");

		SubunitClustererParameters cp = new SubunitClustererParameters();
		cp.setClustererMethod(SubunitClustererMethod.IDENTITY);
		QuatSymmetryParameters symmParams = new QuatSymmetryParameters();

		QuatSymmetryResults symmetry = QuatSymmetryDetector.calcGlobalSymmetry(
				pdb, symmParams, cp);

		assertEquals("C2", symmetry.getSymmetry());
		assertEquals("A2", symmetry.getStoichiometry());
		assertFalse(symmetry.isPseudoStoichiometric());
		assertEquals(SubunitClustererMethod.IDENTITY, symmetry.getSubunitClusters().get(0).getClustererMethod());
		
	}
}
