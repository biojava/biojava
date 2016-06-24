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
package org.biojava.nbio.structure.symmetry.core;

import static org.junit.Assert.*;

import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.cluster.SubunitClustererMethod;
import org.biojava.nbio.structure.cluster.SubunitClustererParameters;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetryDetector;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetryParameters;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetryResults;
import org.junit.Test;

/**
 * Test the algorithm for symmetry detection on a variety of structures with
 * different symmetries.
 * 
 * @author Peter Rose
 * @author Aleix Lafita
 *
 */
public class TestQuatSymmetryDetector {

	/**
	 * An NMR multi-model entry: 1B4C
	 * 
	 * @throws StructureException
	 * @throws IOException
	 */
	@Test
	public void testNMR() throws IOException, StructureException {

		AtomCache cache = new AtomCache();
		cache.setUseMmCif(true);
		Structure pdb = cache.getStructure("BIO:1b4c:1");

		SubunitClustererParameters clusterParams = new SubunitClustererParameters();
		QuatSymmetryParameters symmParams = new QuatSymmetryParameters();
		QuatSymmetryResults symmetry = QuatSymmetryDetector.calcGlobalSymmetry(
				pdb, symmParams, clusterParams);

		// C2 symmetry non pseudosymmetric
		assertEquals("C2", symmetry.getSymmetry());
		assertEquals("A2", symmetry.getSubunits().getStoichiometry());
		assertFalse(symmetry.getSubunits().isPseudoStoichiometric());

	}

	/**
	 * Hemoglobin has both symmetry and pseudosymmetry: 4HHB
	 * 
	 * @throws StructureException
	 * @throws IOException
	 */
	@Test
	public void testPseudosymmetry() throws IOException, StructureException {

		AtomCache cache = new AtomCache();
		cache.setUseMmCif(true);
		Structure pdb = cache.getStructure("BIO:4hhb:1");

		// Non-pseudosymmetry
		SubunitClustererParameters clusterParams = new SubunitClustererParameters();
		QuatSymmetryParameters symmParams = new QuatSymmetryParameters();
		QuatSymmetryResults symmetry = QuatSymmetryDetector.calcGlobalSymmetry(
				pdb, symmParams, clusterParams);

		// C2 symmetry
		assertEquals("C2", symmetry.getSymmetry());
		assertEquals("A2B2", symmetry.getSubunits().getStoichiometry());

		// Use pseudosymmetry (structural clustering)
		clusterParams.setClustererMethod(SubunitClustererMethod.STRUCTURE);
		symmetry = QuatSymmetryDetector.calcGlobalSymmetry(pdb, symmParams,
				clusterParams);

		// D2 pseudo-symmetry
		assertEquals("D2", symmetry.getSymmetry());
		assertEquals("A4", symmetry.getSubunits().getStoichiometry());
	}

	/**
	 * A structure with no global symmetry, but local symmetry: 4P2C
	 * 
	 * @throws IOException
	 * @throws StructureException
	 */
	@Test
	public void testLocal() throws IOException, StructureException {

		AtomCache cache = new AtomCache();
		cache.setUseMmCif(true);
		Structure pdb = cache.getStructure("BIO:4p2c:1");

		// Global Symmetry
		SubunitClustererParameters clusterParams = new SubunitClustererParameters();
		QuatSymmetryParameters symmParams = new QuatSymmetryParameters();
		QuatSymmetryResults symmetry = QuatSymmetryDetector.calcGlobalSymmetry(
				pdb, symmParams, clusterParams);

		// C1 global symmetry
		assertEquals("C1", symmetry.getSymmetry());
		assertEquals("AB5C5", symmetry.getSubunits().getStoichiometry());

		// Local symmetry
		List<QuatSymmetryResults> localSymm = QuatSymmetryDetector
				.calcLocalSymmetries(pdb, symmParams, clusterParams);

		// C5 local symmetry excluding chain A
		assertEquals(localSymm.size(), 3);
		assertEquals("C5", localSymm.get(0).getSymmetry());
		assertEquals("C5", localSymm.get(1).getSymmetry());
		assertEquals("C5", localSymm.get(2).getSymmetry());

		// Two A5 and one A5B5 stoichiometries as local symmetry
		List<String> stoich = localSymm.stream().map(
				t -> t.getSubunits().getStoichiometry()).collect(Collectors.toList());

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

		AtomCache cache = new AtomCache();
		cache.setUseMmCif(true);
		Structure pdb = cache.getStructure("BIO:4e3e:1");

		// Internal symmetry analysis
		SubunitClustererParameters cp = new SubunitClustererParameters();
		cp.setClustererMethod(SubunitClustererMethod.INTERNAL_SYMMETRY);
		cp.setCoverageThreshold(0.8); // Lower coverage for internal symm

		QuatSymmetryParameters symmParams = new QuatSymmetryParameters();
		QuatSymmetryResults symmetry = QuatSymmetryDetector.calcGlobalSymmetry(
				pdb, symmParams, cp);

		// D2 combined internal and quaternary symmetry
		assertEquals("D3", symmetry.getSymmetry());
		assertEquals("A6", symmetry.getSubunits().getStoichiometry());

	}
}
