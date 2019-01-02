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

import org.biojava.nbio.structure.*;
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

		for (int indexToTest =0; indexToTest < asaCalc.getAtomCoords().length; indexToTest++) {
			//int indexToTest = 198;

			Integer[] nbsSh = asaCalc.findNeighborIndicesSpatialHashing(indexToTest);

			Integer[] nbs = asaCalc.findNeighborIndices(indexToTest);

			int countNotInNbs = 0;
			List<Integer> listOfMatchingIndices = new ArrayList<>();
			for (int i = 0; i < nbsSh.length; i++) {
				boolean contained = false;
				for (int j = 0; j < nbs.length; j++) {
					if (nbs[j].equals(nbsSh[i])) {
						listOfMatchingIndices.add(j);
						contained = true;
						break;
					}
				}
				if (!contained) {
					countNotInNbs++;
				}
			}

			//System.out.println("In nbsSh but not in nbs: " + countNotInNbs);
			//System.out.println("Number of matching indices: " + listOfMatchingIndices.size());

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

		Structure structure = StructureIO.getStructure("4F5X");
		Chain c = structure.getPolyChainByPDB("W");
		Atom[] atoms = StructureTools.getAllAtomArray(c);
		System.out.printf("Total of %d atoms\n", atoms.length);

		int nThreads = 1;
		// 1. WITH SPATIAL HASHING

		long start = System.currentTimeMillis();
		AsaCalculator asaCalc = new AsaCalculator(atoms,
				AsaCalculator.DEFAULT_PROBE_SIZE,
				100, nThreads);
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
				100, nThreads);
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
}
