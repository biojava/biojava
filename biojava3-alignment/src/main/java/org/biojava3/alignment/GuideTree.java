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
 * Created on July 1, 2010
 * Author: Mark Chapman
 */

package org.biojava3.alignment;

import java.util.Enumeration;
import java.util.List;
import java.util.Vector;
import javax.swing.tree.TreeNode;

import org.biojava3.alignment.Alignments.PairwiseScorer;
import org.biojava3.alignment.template.GapPenalty;
import org.biojava3.alignment.template.PairwiseSequenceScorer;
import org.biojava3.alignment.template.Profile;
import org.biojava3.alignment.template.SubstitutionMatrix;
import org.biojava3.core.sequence.AccessionID;
import org.biojava3.core.sequence.template.Compound;
import org.biojava3.core.sequence.template.Sequence;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogenyinference.BasicSymmetricalDistanceMatrix;
import org.forester.phylogenyinference.NeighborJoining;

public class GuideTree<S extends Sequence<C>, C extends Compound> {

    private List<S> sequences;
    private List<PairwiseSequenceScorer<S, C>> scorers;
    private GapPenalty gapPenalty;
    private SubstitutionMatrix<C> subMatrix;
    private BasicSymmetricalDistanceMatrix distances;
    private Phylogeny phylogeny;
    private Node root;

    public GuideTree(List<S> sequences, PairwiseScorer type, GapPenalty gapPenalty, SubstitutionMatrix<C> subMatrix) {
        this.sequences = sequences;
        Alignments.runScorers(scorers = Alignments.getScorerList(sequences, type, gapPenalty, subMatrix));
        this.gapPenalty = gapPenalty;
        this.subMatrix = subMatrix;
        distances = new BasicSymmetricalDistanceMatrix(sequences.size());
        for (int i = 0, n = 0; i < sequences.size(); i++) {
            AccessionID id = sequences.get(i).getAccession();
            distances.setIdentifier(i, (id == null) ? Integer.toString(i + 1) : id.getID());
            for (int j = i+1; j < sequences.size(); j++) {
                PairwiseSequenceScorer<S, C> scorer = scorers.get(n++);
                distances.setValue(i, j, (double)(scorer.getMaxScore() - scorer.getScore()) / (scorer.getMaxScore()
                        - scorer.getMinScore()));
            }
        }
        phylogeny = NeighborJoining.createInstance().execute(distances);
        root = new Node(phylogeny.getRoot());
    }

    public int[] getAllPairsScores() {
        int[] scores = new int[scorers.size()];
        int n = 0;
        for (PairwiseSequenceScorer<S, C> scorer : scorers) {
            scores[n++] = scorer.getScore();
        }
        return scores;
    }

    public double[][] getDistanceMatrix() {
        double[][] matrix = new double[distances.getSize()][distances.getSize()];
        for (int i = 0; i < matrix.length; i++) {
            for (int j = i+1; j < matrix.length; j++) {
                matrix[i][j] = matrix[j][i] = distances.getValue(i, j);
            }
        }
        return matrix;
    }

    public Node getRoot() {
        return root;
    }

    public int[][] getScoreMatrix() {
        int[][] matrix = new int[sequences.size()][sequences.size()];
        for (int i = 0, n = 0; i < matrix.length; i++) {
            matrix[i][i] = scorers.get(i).getMaxScore();
            for (int j = i+1; j < matrix.length; j++) {
                matrix[i][j] = matrix[j][i] = scorers.get(n++).getScore();
            }
        }
        return matrix;
    }

    @Override
    public String toString() {
        return phylogeny.toString();
    }

    public class Node implements TreeNode {

        private PhylogenyNode node;

        private Node(PhylogenyNode node) {
            this.node = node;
        }

        public Node getChild1() {
            return new Node(node.getChildNode1());
        }

        public Node getChild2() {
            return new Node(node.getChildNode2());
        }

        public double getLengthToParent() {
            return node.getDistanceToParent();
        }

        public Node getParent() {
            return new Node(node.getParent());
        }

        public Profile<S, C> getProfile() {
            // TODO GuideTree.Node.getProfile parallel version
            return isLeaf() ? new SimpleProfile<S, C>(sequences.get(distances.getIndex(node.getNodeName()))) :
                new SimpleProfileProfileAligner<S, C>(getChild1().getProfile(), getChild2().getProfile(), gapPenalty,
                        subMatrix).getPair();
        }

        public boolean isLeaf() {
            return node.isExternal();
        }

        // methods for TreeNode

        @Override
        public Enumeration<Node> children() {
            Vector<Node> children = new Vector<Node>();
            for (PhylogenyNode n : node.getDescendants()) {
                children.add(new Node(n));
            }
            return children.elements();
        }

        @Override
        public boolean getAllowsChildren() {
            return !isLeaf();
        }

        @Override
        public TreeNode getChildAt(int childIndex) {
            return new Node(node.getChildNode(childIndex));
        }

        @Override
        public int getChildCount() {
            return node.getNumberOfDescendants();
        }

        @Override
        public int getIndex(TreeNode child) {
            if (child instanceof GuideTree<?, ?>.Node) {
                Node c = (Node) child;
                int i = 0;
                for (PhylogenyNode n : node.getDescendants()) {
                    if (n.equals(c.node)) {
                        return i;
                    } else {
                        i++;
                    }
                }
            }
            return -1;
        }

    }

}
