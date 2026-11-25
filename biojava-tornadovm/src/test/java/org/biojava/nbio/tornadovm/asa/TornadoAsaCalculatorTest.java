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
package org.biojava.nbio.tornadovm.asa;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.asa.AsaCalculator;
import org.biojava.nbio.structure.asa.GroupAsa;
import org.biojava.nbio.tornadovm.TestStructures;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;

class TornadoAsaCalculatorTest {

	@ParameterizedTest
	@ValueSource(strings = { "1SMT", "4HHB", "1CDG" })
	void perAtomAsasAreBitIdentical(String id) throws Exception {
		Atom[] atoms = StructureTools.getAllNonHAtomArray(TestStructures.load(id), false, 0);
		double[] expected = new AsaCalculator(atoms, AsaCalculator.DEFAULT_PROBE_SIZE, 1000, 1).calculateAsas();
		double[] actual = TestStructures.forced(
				() -> new TornadoAsaCalculator(atoms, AsaCalculator.DEFAULT_PROBE_SIZE, 1000, 1).calculateAsas());
		assertArrayEquals(expected, actual, 0.0);
	}

	/** The single-precision kernel (devices without FP64) with its exact re-checks gives the same result. */
	@ParameterizedTest
	@ValueSource(strings = { "1SMT", "4HHB", "1CDG" })
	void singlePrecisionPathIsBitIdentical(String id) throws Exception {
		Atom[] atoms = StructureTools.getAllNonHAtomArray(TestStructures.load(id), false, 0);
		double[] expected = new AsaCalculator(atoms, AsaCalculator.DEFAULT_PROBE_SIZE, 1000, 1).calculateAsas();
		System.setProperty("biojava.tornado.asa.precision", "fp32");
		try {
			double[] actual = TestStructures.forced(
					() -> new TornadoAsaCalculator(atoms, AsaCalculator.DEFAULT_PROBE_SIZE, 1000, 1).calculateAsas());
			assertArrayEquals(expected, actual, 0.0);
		} finally {
			System.clearProperty("biojava.tornado.asa.precision");
		}
	}

	@ParameterizedTest
	@ValueSource(ints = { 30, 100, 960 })
	void otherSpherePointCounts(int nPoints) throws Exception {
		Atom[] atoms = StructureTools.getAllNonHAtomArray(TestStructures.load("1SMT"), false, 0);
		double[] expected = new AsaCalculator(atoms, 1.4, nPoints, 1).calculateAsas();
		double[] actual = TestStructures.forced(() -> new TornadoAsaCalculator(atoms, 1.4, nPoints, 1).calculateAsas());
		assertArrayEquals(expected, actual, 0.0);
	}

	@ParameterizedTest
	@ValueSource(strings = { "4HHB" })
	void groupAsasMatch(String id) throws Exception {
		Structure s = TestStructures.load(id);
		GroupAsa[] expected = new AsaCalculator(s, 1.4, 1000, 1, true).getGroupAsas();
		GroupAsa[] actual = TestStructures.forced(() -> new TornadoAsaCalculator(s, 1.4, 1000, 1, true).getGroupAsas());
		assertEquals(expected.length, actual.length);
		for (int i = 0; i < expected.length; i++) {
			assertEquals(expected[i].getGroup().getResidueNumber(), actual[i].getGroup().getResidueNumber());
			assertEquals(expected[i].getAsaU(), actual[i].getAsaU(), 0.0);
		}
	}
}
