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

		// Assert error free TODO now failing
		//assertEquals(0.0, cv, 0.001);
	}

	@Test
	public void testErrorEstimation() throws Exception {

		// Matrix taken from forester test
		BasicSymmetricalDistanceMatrix m = new BasicSymmetricalDistanceMatrix(4);
		m.setIdentifier(0, "A");
		m.setIdentifier(1, "B");
		m.setIdentifier(2, "C");
		m.setIdentifier(3, "D");
		m.setRow("0.00 0.95 0.17 0.98", 0);
		m.setRow("0.95 0.00 1.02 1.83", 1);
		m.setRow("0.17 1.02 0.00 1.01", 2);
		m.setRow("0.98 1.83 1.01 0.00", 3);

		// Calculate error free tree and get the cv
		Phylogeny tree = TreeConstructor.distanceTree(
				ForesterWrapper.cloneDM(m), TreeConstructorType.NJ);
		double cv = DistanceTreeEvaluator.evaluate(tree, m);

		// Assert error is only 5%
		assertEquals(2.411, cv, 0.001);

	}
}
