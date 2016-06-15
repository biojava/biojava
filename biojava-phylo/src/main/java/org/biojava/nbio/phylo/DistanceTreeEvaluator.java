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

import org.forester.evoinference.matrix.distance.DistanceMatrix;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Check the accuracy of a Distance Tree by least squares error (LSE) of the
 * Tree branch lengths and the original Distance Matrix.
 *
 * @author Scooter Willis
 * @author Aleix Lafita
 *
 */
public class DistanceTreeEvaluator {

	private static final Logger logger = LoggerFactory
			.getLogger(DistanceTreeEvaluator.class);

	/** Prevent instantiation */
	private DistanceTreeEvaluator() {
	}

	/**
	 * Evaluate the goodness of fit of a given tree to the original distance
	 * matrix. The returned value is the coefficient of variation, i.e. the
	 * square root of the LS error normalized by the mean.
	 * <p>
	 * This measure can also give an estimate of the quality of the distance
	 * matrix, because a bad fit may mean that the distance is non-additive.
	 *
	 * @param tree
	 *            Phylogenetic Distance Tree to evaluate
	 * @param matrix
	 *            Distance Matrix with the original distances
	 * @return the square root of the average tree LS error normalized by the
	 *         average tree distance (coefficient of variation, CV).
	 */
	public static double evaluate(Phylogeny tree, DistanceMatrix matrix) {
		int numSequences = matrix.getSize();
		List<PhylogenyNode> externalNodes = tree.getExternalNodes();
		HashMap<String, PhylogenyNode> externalNodesHashMap = new HashMap<String, PhylogenyNode>();
		Set<PhylogenyNode> path = new HashSet<PhylogenyNode>();

		for (PhylogenyNode node : externalNodes) {
			externalNodesHashMap.put(node.getName(), node);
		}
		int count = 0;
		double averageMatrixDistance = 0.0;
		double averageTreeDistance = 0.0;
		double averageTreeErrorDistance = 0.0;
		for (int row = 0; row < numSequences - 1; row++) {
			String nodeName1 = matrix.getIdentifier(row);
			PhylogenyNode node1 = externalNodesHashMap.get(nodeName1);
			markPathToRoot(node1, path);
			for (int col = row + 1; col < numSequences; col++) {
				count++;
				String nodeName2 = matrix.getIdentifier(col);
				PhylogenyNode node2 = externalNodesHashMap.get(nodeName2);
				double distance = matrix.getValue(col, row);
				averageMatrixDistance = averageMatrixDistance + distance;
				PhylogenyNode commonParent = findCommonParent(node2, path);
				if (commonParent != null) {
					double treeDistance = getNodeDistance(commonParent, node1)
							+ getNodeDistance(commonParent, node2);

					averageTreeDistance += treeDistance;
					averageTreeErrorDistance += (distance - treeDistance)
							* (distance - treeDistance);
					logger.info("{} {} Distance: {}Tree: {} difference: {}",
							nodeName1, nodeName2, distance, treeDistance,
							Math.abs(distance - treeDistance));
				} else {
					logger.warn("Unable to find common parent with {} {}",
							node1, node2);
				}
			}
			path.clear();
		}
		averageMatrixDistance /= count;
		averageTreeDistance /= count;
		averageTreeErrorDistance /= count;

		logger.info("Average matrix distance: {}", averageMatrixDistance);
		logger.info("Average tree distance: {}", averageTreeDistance);
		logger.info("Average LS error: {}", averageTreeErrorDistance);

		return Math.sqrt(averageTreeErrorDistance) / averageMatrixDistance;
	}

	private static double getNodeDistance(PhylogenyNode parentNode,
			PhylogenyNode childNode) {
		double distance = 0.0;
		while (childNode != parentNode) {
			distance = distance + childNode.getDistanceToParent();
			childNode = childNode.getParent();
		}

		return distance;
	}

	private static PhylogenyNode findCommonParent(PhylogenyNode node,
			Set<PhylogenyNode> path) {
		while (!path.contains(node)) {
			node = node.getParent();
		}
		return node;
	}

	private static void markPathToRoot(PhylogenyNode node,
			Set<PhylogenyNode> path) {
		path.add(node);
		while (!node.isRoot()) {
			node = node.getParent();
			path.add(node);
		}
	}
}
