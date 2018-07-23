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
package org.biojava.nbio.structure.test.cluster;

import static org.junit.Assert.*;

import java.io.IOException;
import java.util.List;
import java.util.function.Function;
import java.util.stream.Collectors;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.cluster.SubunitCluster;
import org.biojava.nbio.structure.cluster.SubunitClusterer;
import org.biojava.nbio.structure.cluster.SubunitClustererMethod;
import org.biojava.nbio.structure.cluster.SubunitClustererParameters;
import org.biojava.nbio.structure.symmetry.core.Stoichiometry;
import org.junit.Test;

/**
 * Test the {@link SubunitClusterer} clustering correctness on different real
 * structures with different types of chain clustering difficulties.
 * 
 * @author Aleix Lafita
 *
 */
public class TestSubunitClustererExamples {

	/**
	 * Test modified residues: 1HIV
	 */
	@Test
	public void testPTMs() throws StructureException, IOException {

		Structure s = StructureIO.getStructure("1hiv");
		SubunitClustererParameters params = new SubunitClustererParameters();
		params.setClustererMethod(SubunitClustererMethod.SEQUENCE);

		List<SubunitCluster> clusters = SubunitClusterer.cluster(s, params).getClusters();

		// We expect a single cluster with length 99
		assertEquals(clusters.size(), 1);
		assertEquals(clusters.get(0).length(), 99);
	}


	/**
	 * A couple of structures with large number of subunits
	 *
	 * @throws IOException
	 * @throws StructureException
	 */

	@Test
	public void testStoichiometry() throws IOException, StructureException {

		SubunitClustererParameters clusterParams = new SubunitClustererParameters(true);
		clusterParams.setClustererMethod(SubunitClustererMethod.SEQUENCE);
		clusterParams.setSequenceIdentityThreshold(0.75);

		Structure structure = StructureIO.getStructure("BIO:5WVK:1");

		Stoichiometry composition = SubunitClusterer.cluster(structure,clusterParams);
		// default alphabet and strategy
		assertEquals("A2B2C2D2E2F2G2H2I2J2K2L2M2N2OPQRSTUVWXYZabcdef",composition.toString());
		// nothing should change
		composition.setStrategy(Stoichiometry.StringOverflowStrategy.DOUBLE);
		assertEquals("A2B2C2D2E2F2G2H2I2J2K2L2M2N2OPQRSTUVWXYZabcdef",composition.toString());
		// nothing should change
		composition.setStrategy(Stoichiometry.StringOverflowStrategy.QUESTIONMARK);
		assertEquals("A2B2C2D2E2F2G2H2I2J2K2L2M2N2OPQRSTUVWXYZabcdef",composition.toString());

		// now make the alphabet shorter to force changes in representation
		composition.setAlphabet("ABCDEFGHIJKLMNOPQRSTUVWXYZ");

		composition.setStrategy(Stoichiometry.StringOverflowStrategy.DOUBLE);
		assertEquals("AA2AB2AC2AD2AE2AF2AG2AH2AI2AJ2AK2AL2AM2AN2AO1AP1AQ1AR1AS1AT1AU1AV1AW1AX1AY1AZ1" +
						"BA1BB1BC1BD1BE1BF1",
				composition.toString());

		composition.setStrategy(Stoichiometry.StringOverflowStrategy.QUESTIONMARK);
		assertEquals("A2B2C2D2E2F2G2H2I2J2K2L2M2N2OPQRSTUVWXYZ??????",
				composition.toString());

		// back to default
		composition.setStrategy(Stoichiometry.StringOverflowStrategy.CYCLE);
		assertEquals("A2B2C2D2E2F2G2H2I2J2K2L2M2N2OPQRSTUVWXYZABCDEF",
				composition.toString());

		// manipulations with composition subsets
		Stoichiometry a2 = composition.getComponent(0);
		assertEquals("A2",a2.toString());
		Stoichiometry b2 = composition.getComponent(1);
		assertEquals("B2",b2.toString());
		Stoichiometry c2 = composition.getComponent(2);
		assertEquals("C2",c2.toString());
		Stoichiometry z = composition.getComponent(25);
		assertEquals("Z",z.toString());

		Stoichiometry b2c2 = b2.combineWith(c2);
		assertEquals("B2C2",b2c2.toString());

		Stoichiometry a2z = a2.combineWith(z);
		assertEquals("A2Z",a2z.toString());

		// result must be ordered by subunit count
		Stoichiometry a2b2z = a2z.combineWith(b2);
		assertEquals("A2B2Z",a2b2z.toString());

		//overlapping clusters counted once
		Stoichiometry a2b2c2z = a2b2z.combineWith(b2c2);
		assertEquals("A2B2C2Z",a2b2c2z.toString());

		// custom string representation, comma-separated subunit counts
		Function<List<SubunitCluster>,String> numbersOnly =
				l-> l.stream().
					map(SubunitCluster::size).
					map(String::valueOf).
					collect(Collectors.joining(", "));

		a2b2c2z.setCustomStringGenerator(numbersOnly);
		assertEquals("2, 2, 2, 1", a2b2c2z.toString());

		composition.setCustomStringGenerator(numbersOnly);
		assertEquals("2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1",
				composition.toString());
	}


	/**
	 * Test pseudostoichiometry: 4HHB
	 */
	@Test
	public void testPseudostoichiometry() throws StructureException,
			IOException {

		Structure s = StructureIO.getStructure("4HHB");
		SubunitClustererParameters params = new SubunitClustererParameters();
		params.setClustererMethod(SubunitClustererMethod.SEQUENCE);
		params.setSequenceIdentityThreshold(0.95);

		List<SubunitCluster> clusters = SubunitClusterer.cluster(s, params).getClusters();

		// We expect two SEQUENCE clusters with length 141 and 146
		assertEquals(clusters.size(), 2);
		assertEquals(clusters.get(0).length(), 141);
		assertEquals(clusters.get(1).length(), 146);
		assertEquals(clusters.get(0).getClustererMethod(),
				SubunitClustererMethod.SEQUENCE);
		assertEquals(clusters.get(1).getClustererMethod(),
				SubunitClustererMethod.SEQUENCE);

		params.setClustererMethod(SubunitClustererMethod.SEQUENCE_STRUCTURE);
		params.setRMSDThreshold(3.0);

		clusters = SubunitClusterer.cluster(s, params).getClusters();

		// We expect a single STRUCTURE cluster with length 140
		assertEquals(clusters.size(), 1);
		assertEquals(clusters.get(0).length(), 140, 2);
		assertEquals(clusters.get(0).getClustererMethod(),
				SubunitClustererMethod.STRUCTURE);
	}

	/**
	 * Test internally symmetric: 4E3E bioassembly 1
	 */
	@Test
	public void testInternalSymmetry() throws StructureException, IOException {

		Structure s = StructureIO.getStructure("BIO:4E3E:1");

		SubunitClustererParameters params = new SubunitClustererParameters();
		params.setClustererMethod(SubunitClustererMethod.SEQUENCE);
		params.setSequenceCoverageThreshold(0.8);

		List<SubunitCluster> clusters = SubunitClusterer.cluster(s, params).getClusters();

		// We expect one SEQUENCE cluster with 3 Subunits of length 351
		assertEquals(clusters.size(), 1);
		assertEquals(clusters.get(0).size(), 3);
		assertEquals(clusters.get(0).length(), 351);
		assertEquals(clusters.get(0).getClustererMethod(),
				SubunitClustererMethod.SEQUENCE);

		params.setClustererMethod(SubunitClustererMethod.SEQUENCE_STRUCTURE);
		params.setStructureCoverageThreshold(0.8);
		params.setInternalSymmetry(true);
		params.setRMSDThreshold(3.0);

		clusters = SubunitClusterer.cluster(s, params).getClusters();

		// We expect a single INTERNAL_SYMMETRY cluster with 6 Subunits
		assertEquals(clusters.size(), 1);
		assertEquals(clusters.get(0).size(), 6);
		assertTrue(clusters.get(0).length() < 177);
		assertEquals(clusters.get(0).getClustererMethod(),
				SubunitClustererMethod.STRUCTURE);
	}
}
