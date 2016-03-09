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
package org.biojava.nbio.phylo;

import org.forester.evoinference.matrix.distance.BasicSymmetricalDistanceMatrix;
import org.forester.phylogeny.Phylogeny;
import org.junit.Test;

import static org.junit.Assert.*;

/**
 * Test for the evaluation of distance trees.
 *
 * @author Aleix Lafita
 *
 */
public class TestDistanceTreeEvaluator {

	@Test
	public void testErrorFree() throws Exception {

		// Create a perfect additive distance matrix
		BasicSymmetricalDistanceMatrix DM = new BasicSymmetricalDistanceMatrix(
				3);

		// dAB = 0.8, dBC = 0.4, dAC = 0.8
		for (int i = 0; i < 3; i++) {
			char id = (char) ('A' + i);
			DM.setIdentifier(i, id + "");
			for (int j = i; j < 3; j++) {
				if (i == j) {
					DM.setValue(i, j, 0.0);
				} else if (i == 0) {
					DM.setValue(i, j, 0.8);
				} else {
					DM.setValue(i, j, 0.4);
				}
			}
		}

		// Calculate error free tree and get the cv
		Phylogeny tree = TreeConstructor.distanceTree(
				ForesterWrapper.cloneDM(DM), TreeConstructorType.NJ);
		double cv = DistanceTreeEvaluator.evaluate(tree, DM);

		// Assert error free
		assertEquals(0.0, cv, 0.001);
	}

	@Test
	public void testErrorEstimation() throws Exception {

		// Matrix taken from forester test
		BasicSymmetricalDistanceMatrix m = new BasicSymmetricalDistanceMatrix(4);
		m.setIdentifier(0, "A");
		m.setIdentifier(1, "B");
		m.setIdentifier(2, "C");
		m.setIdentifier(3, "D");
		m.setRow("0 1 0 1", 0);
		m.setRow("1 0 0 1", 1);
		m.setRow("0 0 0 1", 2);
		m.setRow("1 1 1 0", 3);

		// Calculate error free tree and get the cv
		Phylogeny tree = TreeConstructor.distanceTree(
				ForesterWrapper.cloneDM(m), TreeConstructorType.NJ);
		double cv = DistanceTreeEvaluator.evaluate(tree, m);

		// Assert error is about 30%
		assertEquals(0.3, cv, 0.05);

	}
}
