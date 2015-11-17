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
		BasicSymmetricalDistanceMatrix DM = new BasicSymmetricalDistanceMatrix(3);
		
		// dAB = 0.8, dBC = 0.4, dAC = 0.8
		for (int i = 0; i < 3; i++){
			char id = (char) ('A' + i);
			DM.setIdentifier(i, id + "");
			for (int j = i; j < 3; j++){
				if (i == j){
					DM.setValue(i, j, 0.0);
				}
				else if (i == 0) {
					DM.setValue(i, j, 0.8);
					DM.setValue(j, i, 0.8);
				} else {
					DM.setValue(i, j, 0.4);
					DM.setValue(j, i, 0.4);
				}
			}
		}
		
		BasicSymmetricalDistanceMatrix cloneDM = ForesterWrapper.cloneDM(DM);
		
		//Calculate error free tree and get the cv
		Phylogeny tree = TreeConstructor.distanceTree(cloneDM, TreeConstructorType.NJ);
		double cv = DistanceTreeEvaluator.evaluate(tree, DM);

		// Assert error free
		assertEquals(0.0, cv, 0.001);

	}
}
