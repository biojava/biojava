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
import java.util.*;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.cluster.SubunitClusterer;
import org.biojava.nbio.structure.cluster.SubunitClustererMethod;
import org.biojava.nbio.structure.cluster.SubunitClustererParameters;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetryDetector;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetryParameters;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetryResults;
import org.biojava.nbio.structure.symmetry.core.Stoichiometry;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Test the {@link QuatSymmetryDetector} algorithm for symmetry detection on a
 * variety of structures with different symmetries.
 * 
 * @author Peter Rose
 * @author Aleix Lafita
 *
 */
public class TestQuatSymmetryDetectorExamples {

	private static final Logger logger = LoggerFactory.getLogger(TestQuatSymmetryDetectorExamples.class);
	
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
		assertEquals("A2", symmetry.getStoichiometry().toString());
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
		assertEquals("A6B6", symmetry.getStoichiometry().toString());
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
		assertEquals("A2B2", symmetry.getStoichiometry().toString());

		// Use pseudosymmetry (structural clustering)
		clusterParams.setClustererMethod(SubunitClustererMethod.STRUCTURE);
		symmetry = QuatSymmetryDetector.calcGlobalSymmetry(pdb, symmParams,
				clusterParams);

		// D2 pseudo-symmetry
		assertEquals("D2", symmetry.getSymmetry());
		assertEquals("A4", symmetry.getStoichiometry().toString());
		assertTrue(symmetry.isPseudoStoichiometric());
	}

	/**
	 * A selection of structures with no global symmetry, but local symmetry
	 * 
	 * @throws IOException
	 * @throws StructureException
	 */
	@Test
	public void testLocal() throws IOException, StructureException {

		AtomCache atomCache = new AtomCache();
		atomCache.setUseMmtf(true);

		List<String> testIds = new ArrayList<>();
		List<String> testStoichiometries = new ArrayList<>();
		List<Map<String,String>> testLocalSymmetries = new ArrayList<>();
		Map<String,String> localSymmetries;

		testIds.add("BIO:5NUQ:1");
			testStoichiometries.add("A3B");
			localSymmetries = new HashMap<>();
				localSymmetries.put("A3","C3");
			testLocalSymmetries.add(localSymmetries);

		testIds.add("BIO:4P2C:1");
			testStoichiometries.add("A5B5C");
			localSymmetries = new HashMap<>();
				localSymmetries.put("A5B5","C5");
			testLocalSymmetries.add(localSymmetries);

		testIds.add("BIO:3J96:1");
			testStoichiometries.add("A6B4CDE");
			localSymmetries = new HashMap<>();
				localSymmetries.put("B4","C4");
			testLocalSymmetries.add(localSymmetries);

		testIds.add("BIO:5WVK:1");
			testStoichiometries.add("A2B2C2D2E2F2G2H2I2J2K2L2M2N2OPQRSTUVWXYZabcdef");
			localSymmetries = new HashMap<>();
				localSymmetries.put("A2B2C2D2E2F2G2H2I2J2K2L2M2N2","C2");
			testLocalSymmetries.add(localSymmetries);

		testIds.add("BIO:5JXT:1");
			testStoichiometries.add("A16");
			localSymmetries = new HashMap<>();
				localSymmetries.put("A8","D2");
			testLocalSymmetries.add(localSymmetries);

		testIds.add("BIO:3R8R:1");
			testStoichiometries.add("A12");
			localSymmetries = new HashMap<>();
				localSymmetries.put("A10","D5");
			testLocalSymmetries.add(localSymmetries);

		testIds.add("BIO:1O18:1");
			testStoichiometries.add("A14B6C5D5");
			localSymmetries = new HashMap<>();
				localSymmetries.put("A14","H");
			testLocalSymmetries.add(localSymmetries);

		testIds.add("BIO:4A8A:1");
			testStoichiometries.add("A12B");
			localSymmetries = new HashMap<>();
				localSymmetries.put("A12","T");
			testLocalSymmetries.add(localSymmetries);

		testIds.add("BIO:5DN6:1");
			testStoichiometries.add("A12B3C3DEFGHIJKL");
			localSymmetries = new HashMap<>();
				localSymmetries.put("B3C3","C3");
				localSymmetries.put("A12","C12");
			testLocalSymmetries.add(localSymmetries);

		testIds.add("BIO:5FL7:1");
			testStoichiometries.add("A10B3C3DE");
			localSymmetries = new HashMap<>();
				localSymmetries.put("B3C3","C3");
				localSymmetries.put("A10","C10");
			testLocalSymmetries.add(localSymmetries);

		testIds.add("BIO:2OF5:1");
			testStoichiometries.add("A7B5");
			localSymmetries = new HashMap<>();
				localSymmetries.put("A7","H");
				localSymmetries.put("A5B5","H");
			testLocalSymmetries.add(localSymmetries);

		testIds.add("BIO:6EM9:1");
			testStoichiometries.add("A8B2");
			localSymmetries = new HashMap<>();
				localSymmetries.put("B2","C2");
				localSymmetries.put("A7","H");
			testLocalSymmetries.add(localSymmetries);

		testIds.add("BIO:4NTP:1");
			testStoichiometries.add("A16");
			localSymmetries = new HashMap<>();
				localSymmetries.put("A2","C2");
				localSymmetries.put("A6","D3");
			testLocalSymmetries.add(localSymmetries);

		testIds.add("BIO:3JC9:1");
			testStoichiometries.add("A12B12C12D12E12F12G5H2");
			localSymmetries = new HashMap<>();
				localSymmetries.put("A12C12D12E12F12H2","C2");
				localSymmetries.put("A12B12C12D12E12F12","C12");
				localSymmetries.put("G5","H");
			testLocalSymmetries.add(localSymmetries);

		QuatSymmetryParameters symmParams = new QuatSymmetryParameters();
		SubunitClustererParameters clusterParams = new SubunitClustererParameters(true);
		clusterParams.setClustererMethod(SubunitClustererMethod.SEQUENCE);
		clusterParams.setSequenceIdentityThreshold(0.75);

		for(int iTest = 0; iTest<testIds.size();iTest++) {
			logger.info("Processing "+testIds.get(iTest));
			Structure pdb = atomCache.getStructure(testIds.get(iTest));
			Stoichiometry composition = SubunitClusterer.cluster(pdb,clusterParams);

			// no global symmetry
			QuatSymmetryResults globalSymmetry = QuatSymmetryDetector.calcGlobalSymmetry(composition,symmParams);
			assertEquals("Unexpected global symmetry in "+testIds.get(iTest),
					"C1", globalSymmetry.getSymmetry());
			assertEquals("Unexpected global stoichiometry in "+testIds.get(iTest),
					testStoichiometries.get(iTest), globalSymmetry.getStoichiometry().toString());

			List<QuatSymmetryResults> foundLocal = QuatSymmetryDetector.calcLocalSymmetries(composition, symmParams);
			Map<String,String> refLocal = testLocalSymmetries.get(iTest);

			for (QuatSymmetryResults local:foundLocal) {
				logger.info("Found stoichiometry "+local.getStoichiometry().toString()+" with symmetry "+local.getSymmetry());
				assertTrue("Stoichiometry "+local.getStoichiometry().toString()+" not expected for "+testIds.get(iTest),
						refLocal.keySet().contains(local.getStoichiometry().toString()));

				assertEquals("Symmetry "+local.getSymmetry()+" with stoichiometry "+local.getStoichiometry().toString()+
								" not expected for "+testIds.get(iTest),
						refLocal.get(local.getStoichiometry().toString()),local.getSymmetry());
			}
		}
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
		cp.setClustererMethod(SubunitClustererMethod.SEQUENCE_STRUCTURE);
		cp.setInternalSymmetry(true);
		cp.setStructureCoverageThreshold(0.75); // Lower coverage for internal symm

		QuatSymmetryParameters symmParams = new QuatSymmetryParameters();
		QuatSymmetryResults symmetry = QuatSymmetryDetector.calcGlobalSymmetry(
				pdb, symmParams, cp);

		// D2 combined internal and quaternary symmetry
		assertEquals("D3", symmetry.getSymmetry());
		assertEquals("A6", symmetry.getStoichiometry().toString());

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
		assertEquals("A3", symmetry.getStoichiometry().toString());

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
		assertEquals("A5", results.get(0).getStoichiometry().toString());

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
		cp.setSequenceIdentityThreshold(0.95);
		cp.setSequenceCoverageThreshold(0.95);
		cp.setClustererMethod(SubunitClustererMethod.SEQUENCE);
		QuatSymmetryParameters symmParams = new QuatSymmetryParameters();

		QuatSymmetryResults symmetry = QuatSymmetryDetector.calcGlobalSymmetry(
				pdb, symmParams, cp);

		assertEquals("C2", symmetry.getSymmetry());
		assertEquals("A2", symmetry.getStoichiometry().toString());
		assertFalse(symmetry.isPseudoStoichiometric());
		assertEquals(SubunitClustererMethod.SEQUENCE, symmetry.getSubunitClusters().get(0).getClustererMethod());
		
	}
}
