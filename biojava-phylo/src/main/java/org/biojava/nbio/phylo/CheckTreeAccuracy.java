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
/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava.nbio.phylo;

import org.forester.evoinference.matrix.distance.BasicSymmetricalDistanceMatrix;
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
 *
 * @author Scooter
 */
public class CheckTreeAccuracy {

	private static final Logger logger = LoggerFactory.getLogger(CheckTreeAccuracy.class);

    public static DistanceMatrix copyMatrix(DistanceMatrix matrix) {

        DistanceMatrix distanceMatrix = new BasicSymmetricalDistanceMatrix(matrix.getSize());
        for (int i = 0; i < matrix.getSize(); i++) {
            distanceMatrix.setIdentifier(i, matrix.getIdentifier(i));
        }

        for (int col = 0; col < matrix.getSize(); col++) {
            for (int row = 0; row < matrix.getSize(); row++) {
                distanceMatrix.setValue(col, row, matrix.getValue(col, row));
            }
        }

        return distanceMatrix;

    }

    public void process(Phylogeny tree, DistanceMatrix matrix) {
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
                    double treeDistance = getNodeDistance(commonParent, node1) + getNodeDistance(commonParent, node2);

                    averageTreeDistance = averageTreeDistance + treeDistance;
                    averageTreeErrorDistance = averageTreeErrorDistance + Math.abs(distance - treeDistance);
                    logger.info("{} {} Distance: {}Tree: {} difference: {}", nodeName1, nodeName2, distance, treeDistance, Math.abs(distance - treeDistance));
                } else {
                    logger.warn("Unable to find common parent with {} {}", node1, node2);
                }
            }
            path.clear();
        }

        logger.info("Average matrix distance: {}", averageMatrixDistance / count);
        logger.info("Average tree distance: {}", averageTreeDistance / count);
        logger.info("Average error: {}", averageTreeErrorDistance / count);
    }

    public double getNodeDistance(PhylogenyNode parentNode, PhylogenyNode childNode) {
        double distance = 0.0;
        while (childNode != parentNode) {
            distance = distance + childNode.getDistanceToParent();
            childNode = childNode.getParent();
        }

        return distance;
    }

    public PhylogenyNode findCommonParent(PhylogenyNode node, Set<PhylogenyNode> path) {
        while (!path.contains(node)) {
            node = node.getParent();
        }
        return node;
    }

    public void markPathToRoot(PhylogenyNode node, Set<PhylogenyNode> path) {
        path.add(node);
        while (!node.isRoot()) {
            node = node.getParent();
            path.add(node);
        }
    }
}
