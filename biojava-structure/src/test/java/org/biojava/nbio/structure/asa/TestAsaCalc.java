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
package org.biojava.nbio.structure.asa;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.io.mmcif.ChemCompGroupFactory;
import org.biojava.nbio.structure.io.mmcif.DownloadChemCompProvider;
import static org.junit.Assert.*;

import org.junit.Test;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * Testing of Accessible Surface Area calculations
 *
 *
 * @author Jose Duarte
 *
 */
public class TestAsaCalc {


	@Test
	public void testAsa3PIU() throws StructureException, IOException {

		// important: without this the tests can fail when running in maven (but not in IDE)
		// that's because it depends on the order on how tests were run - JD 2018-03-10
		ChemCompGroupFactory.setChemCompProvider(new DownloadChemCompProvider());

		Structure structure = StructureIO.getStructure("3PIU");


		AsaCalculator asaCalc = new AsaCalculator(structure,
				AsaCalculator.DEFAULT_PROBE_SIZE,
				1000, 1, false);

		double totResidues = 0;
		double totAtoms = 0;

		GroupAsa[] groupAsas = asaCalc.getGroupAsas();

		double[] asas = asaCalc.calculateAsas();

		for (double asa:asas) {
			totAtoms += asa;
		}

		for (GroupAsa groupAsa: groupAsas) {
			//System.out.println(groupAsa.getGroup().getPDBName() + " " + groupAsa.getGroup().getResidueNumber() + " " + groupAsa.getAsaU());
			totResidues+=groupAsa.getAsaU();

			assertTrue(groupAsa.getRelativeAsaU() <= 1.0);
		}

		assertEquals(totAtoms, totResidues, 0.000001);

		assertEquals(17462.0, totAtoms, 1.0);

	}

	@Test
	public void testNeighborIndicesFinding() throws StructureException, IOException {
		// important: without this the tests can fail when running in maven (but not in IDE)
		// that's because it depends on the order on how tests were run - JD 2018-03-10
		ChemCompGroupFactory.setChemCompProvider(new DownloadChemCompProvider());

		Structure structure = StructureIO.getStructure("3PIU");

		AsaCalculator asaCalc = new AsaCalculator(structure,
				AsaCalculator.DEFAULT_PROBE_SIZE,
				1000, 1, false);

		int[][] allNbsSh = asaCalc.findNeighborIndicesSpatialHashing();

		int[][] allNbs = asaCalc.findNeighborIndices();

		for (int indexToTest =0; indexToTest < asaCalc.getAtomCoords().length; indexToTest++) {
			//int indexToTest = 198;
			int[] nbsSh = allNbsSh[indexToTest];
			int[] nbs = allNbs[indexToTest];

			List<Integer> listOfMatchingIndices = new ArrayList<>();
			for (int i = 0; i < nbsSh.length; i++) {
				for (int j = 0; j < nbs.length; j++) {
					if (nbs[j] == nbsSh[i]) {
						listOfMatchingIndices.add(j);
						break;
					}
				}
			}
			
//		for (int i = 0; i<nbs.length; i++) {
//			double dist = asaCalc.getAtomCoords()[i].distance(asaCalc.getAtomCoords()[indexToTest]);
//			if (listOfMatchingIndices.contains(i)) {
//				System.out.printf("Matching     - indices %d-%d: %5.2f\n", indexToTest, i, dist);
//			} else {
//				System.out.printf("Not matching - indices %d-%d: %5.2f\n", indexToTest, i, dist);
//			}
//		}

			assertEquals(nbs.length, nbsSh.length);

			assertEquals(nbs.length, listOfMatchingIndices.size());
		}

	}

	@Test
	public void testPerformance() throws StructureException, IOException {
		// important: without this the tests can fail when running in maven (but not in IDE)
		// that's because it depends on the order on how tests were run - JD 2018-03-10
		ChemCompGroupFactory.setChemCompProvider(new DownloadChemCompProvider());

		Structure structure = StructureIO.getStructure("3HBX");
		Atom[] atoms = StructureTools.getAllAtomArray(structure);
		System.out.printf("Total of %d atoms. n(n-1)/2= %d \n", atoms.length, atoms.length*(atoms.length-1)/2);

		int nThreads = 1;
		int nSpherePoints = 100;

		// 1. WITH SPATIAL HASHING
		long start = System.currentTimeMillis();
		AsaCalculator asaCalc = new AsaCalculator(atoms,
				AsaCalculator.DEFAULT_PROBE_SIZE,
				nSpherePoints, nThreads);
		asaCalc.setUseSpatialHashingForNeighbors(true);

		double[] asas = asaCalc.calculateAsas();
		long end = System.currentTimeMillis();
		System.out.printf("ASA calculation took %6.2f s with spatial hashing\n", (end-start)/1000.0);

		double totAtoms = 0;
		for (double asa:asas) {
			totAtoms += asa;
		}
		double withSH = totAtoms;
		System.out.printf("Total ASA is %6.2f \n", totAtoms);


		// 2. WITHOUT SPATIAL HASHING
		start = System.currentTimeMillis();
		asaCalc = new AsaCalculator(atoms,
				AsaCalculator.DEFAULT_PROBE_SIZE,
				nSpherePoints, nThreads);
		asaCalc.setUseSpatialHashingForNeighbors(false);

		asas = asaCalc.calculateAsas();
		end = System.currentTimeMillis();
		System.out.printf("ASA calculation took %6.2f s without spatial hashing\n", (end-start)/1000.0);

		totAtoms = 0;
		for (double asa:asas) {
			totAtoms += asa;
		}
		double withoutSH = totAtoms;
		System.out.printf("Total ASA is %6.2f \n", totAtoms);


		assertEquals(withoutSH, withSH, 0.000001);

	}

	@Test
	public void testNoNeighborsIssue() throws StructureException, IOException {
		// important: without this the tests can fail when running in maven (but not in IDE)
		// that's because it depends on the order on how tests were run - JD 2018-03-10
		ChemCompGroupFactory.setChemCompProvider(new DownloadChemCompProvider());

		Structure structure = StructureIO.getStructure("1EMI");

		// chain B, atom with index 35 does not have any neighbors. The calculation should not fail
		Atom[] atoms = StructureTools.getAllNonHAtomArray(structure.getPolyChainByPDB("B"), false);

		AsaCalculator asaCalc = new AsaCalculator(atoms,
				AsaCalculator.DEFAULT_PROBE_SIZE,
				1000, 1);

		int[][] allNbsSh = asaCalc.findNeighborIndicesSpatialHashing();

		int[][] allNbs = asaCalc.findNeighborIndices();

		for (int indexToTest =0; indexToTest < asaCalc.getAtomCoords().length; indexToTest++) {
			//int indexToTest = 198;
			int[] nbsSh = allNbsSh[indexToTest];
			int[] nbs = allNbs[indexToTest];

			List<Integer> listOfMatchingIndices = new ArrayList<>();
			for (int i = 0; i < nbsSh.length; i++) {
				for (int j = 0; j < nbs.length; j++) {
					if (nbs[j] == nbsSh[i]) {
						listOfMatchingIndices.add(j);
						break;
					}
				}
			}

			assertEquals(nbs.length, nbsSh.length);

			assertEquals(nbs.length, listOfMatchingIndices.size());
		}

	}

	@Test
	public void testNoAtomsAsaCalc() {

		// in case of no atoms at all, the calculation should not fail and return an empty array
		Atom[] atoms = new Atom[0];

		AsaCalculator asaCalc = new AsaCalculator(atoms,
				AsaCalculator.DEFAULT_PROBE_SIZE,
				1000, 1);
		double[] asas = asaCalc.calculateAsas();
		assertNotNull(asas);
		assertEquals(0, asas.length);

	}
}
