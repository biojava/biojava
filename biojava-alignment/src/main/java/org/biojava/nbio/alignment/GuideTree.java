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

package org.biojava.nbio.alignment;

import org.biojava.nbio.core.alignment.SimpleProfile;
import org.biojava.nbio.alignment.template.GuideTreeNode;
import org.biojava.nbio.alignment.template.PairwiseSequenceScorer;
import org.biojava.nbio.core.alignment.template.Profile;
import org.biojava.nbio.core.alignment.template.ProfilePair;
import org.biojava.nbio.core.sequence.AccessionID;
import org.biojava.nbio.core.sequence.template.Compound;
import org.biojava.nbio.core.sequence.template.Sequence;
import org.forester.evoinference.distance.NeighborJoining;
import org.forester.evoinference.matrix.distance.BasicSymmetricalDistanceMatrix;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;

import javax.swing.tree.TreeNode;
import java.util.*;
import java.util.concurrent.Future;

/**
 * Implements a data structure for a guide tree used during progressive multiple sequence alignment.  Leaf
 * {@link Node}s correspond to single {@link Sequence}s.  Internal {@link Node}s correspond to multiple sequence
 * alignments.  The root {@link Node} corresponds to the full multiple sequence alignment.
 *
 * @author Mark Chapman
 * @param <S> each {@link Sequence} in the tree is of type S
 * @param <C> each element of a {@link Sequence} is a {@link Compound} of type C
 */
public class GuideTree<S extends Sequence<C>, C extends Compound> implements Iterable<GuideTreeNode<S, C>> {

    private List<S> sequences;
    private List<PairwiseSequenceScorer<S, C>> scorers;
    private BasicSymmetricalDistanceMatrix distances;
    private String newick;
    private Node root;

    /**
     * Creates a guide tree for use during progressive multiple sequence alignment.
     *
     * @param sequences the {@link List} of {@link Sequence}s to align
     * @param scorers list of sequence pair scorers, one for each pair of sequences given
     */
    public GuideTree(List<S> sequences, List<PairwiseSequenceScorer<S, C>> scorers) {
        this.sequences = Collections.unmodifiableList(sequences);
        this.scorers = Collections.unmodifiableList(scorers);
        distances = new BasicSymmetricalDistanceMatrix(sequences.size());
        BasicSymmetricalDistanceMatrix distclone = new BasicSymmetricalDistanceMatrix(sequences.size());
        for (int i = 0, n = 0; i < sequences.size(); i++) {
            AccessionID id = sequences.get(i).getAccession();
            String str = (id == null) ? Integer.toString(i + 1) : id.getID();
            distances.setIdentifier(i, str);
            distclone.setIdentifier(i, str);
            for (int j = i+1; j < sequences.size(); j++) {
                double dist = scorers.get(n++).getDistance();
                distances.setValue(i, j, dist);
                distclone.setValue(i, j, dist);
            }
        }
        // TODO UPGMA and other hierarchical clustering routines
        Phylogeny phylogeny = NeighborJoining.createInstance().execute(distclone);
        newick = phylogeny.toString();
        root = new Node(phylogeny.getRoot(), null);
    }

    /**
     * Returns a sequence pair score for all {@link Sequence} pairs in the given {@link List}.
     *
     * @return list of sequence pair scores
     */
    public double[] getAllPairsScores() {
    	double[] scores = new double[scorers.size()];
        int n = 0;
        for (PairwiseSequenceScorer<S, C> scorer : scorers) {
            scores[n++] = scorer.getScore();
        }
        return scores;
    }

    /**
     * Returns the distance matrix used to construct this guide tree.  The scores have been normalized.
     *
     * @return the distance matrix used to construct this guide tree
     */
    public double[][] getDistanceMatrix() {
        double[][] matrix = new double[distances.getSize()][distances.getSize()];
        for (int i = 0; i < matrix.length; i++) {
            for (int j = i+1; j < matrix.length; j++) {
                matrix[i][j] = matrix[j][i] = distances.getValue(i, j);
            }
        }
        return matrix;
    }

    /**
     * Returns the root {@link Node} which corresponds to the full multiple sequence alignment.
     *
     * @return the root node
     */
    public Node getRoot() {
        return root;
    }

    /**
     * Returns the similarity matrix used to construct this guide tree.  The scores have not been normalized.
     *
     * @return the similarity matrix used to construct this guide tree
     */
    public double[][] getScoreMatrix() {
    	double[][] matrix = new double[sequences.size()][sequences.size()];
        for (int i = 0, n = 0; i < matrix.length; i++) {
            matrix[i][i] = scorers.get(i).getMaxScore();
            for (int j = i+1; j < matrix.length; j++) {
                matrix[i][j] = matrix[j][i] = scorers.get(n++).getScore();
            }
        }
        return matrix;
    }

    /**
     * Returns the {@link Sequence}s which make up the leaves of this tree.
     *
     * @return the sequences which make up the leaves of this tree
     */
    public List<S> getSequences() {
        return sequences;
    }

    // method for Iterable

    /**
     * Returns a post-order {@link Iterator} that traverses the tree from leaves to root.
     */
    @Override
    public Iterator<GuideTreeNode<S, C>> iterator() {
        return new PostOrderIterator();
    }

    // method from Object

    @Override
    public String toString() {
        return newick;
    }

    /**
     * Implements a data structure for the node in a guide tree used during progressive multiple sequence alignment.
     */
    public class Node implements GuideTreeNode<S, C> {

        private GuideTreeNode<S, C> parent, child1, child2;
        private double distance;
        private String name;
        private boolean isLeaf, isVisited;
        private Profile<S, C> profile;
        private Future<ProfilePair<S, C>> profileFuture;

        private Node(PhylogenyNode node, Node parent) {
            this.parent = parent;
            distance = node.getDistanceToParent();
            name = node.getName();
            if(isLeaf = node.isExternal()) {
                profile = new SimpleProfile<S, C>(sequences.get(distances.getIndex(name)));
            } else {
                child1 = new Node(node.getChildNode1(), this);
                child2 = new Node(node.getChildNode2(), this);
            }
        }

        // methods for GuideTreeNode

        @Override
        public GuideTreeNode<S, C> getChild1() {
            return child1;
        }

        @Override
        public GuideTreeNode<S, C> getChild2() {
            return child2;
        }

        @Override
        public double getDistanceToParent() {
            return distance;
        }

        @Override
        public String getName() {
            return name;
        }

        @Override
        public Profile<S, C> getProfile() {
            return profile;
        }

        @Override
        public Future<ProfilePair<S, C>> getProfileFuture() {
            return profileFuture;
        }

        @Override
        public void setProfile(Profile<S, C> profile) {
            this.profile = profile;
            profileFuture = null;
        }

        @Override
        public void setProfileFuture(Future<ProfilePair<S, C>> profileFuture) {
            this.profileFuture = profileFuture;
            profile = null;
        }

        // methods for TreeNode

        @Override
        public Enumeration<GuideTreeNode<S, C>> children() {
            Vector<GuideTreeNode<S, C>> children = new Vector<GuideTreeNode<S, C>>();
            children.add(getChild1());
            children.add(getChild2());
            return children.elements();
        }

        @Override
        public boolean getAllowsChildren() {
            return !isLeaf();
        }

        @Override
        public GuideTreeNode<S, C> getChildAt(int childIndex) {
            if (childIndex == 1) {
                return getChild1();
            } else if (childIndex == 2) {
                return getChild2();
            }
            throw new IndexOutOfBoundsException();
        }

        @Override
        public int getChildCount() {
            return 2;
        }

        @Override
        public int getIndex(TreeNode child) {
            return getChildAt(1) == child ? 1 : (getChildAt(2) == child ? 2 : -1);
        }

        @Override
        public GuideTreeNode<S, C> getParent() {
            return parent;
        }

        @Override
        public boolean isLeaf() {
            return isLeaf;
        }

        // helper methods for iterator

        private void clearVisited() {
            isVisited = false;
            if (!isLeaf()) {
                ((Node) getChild1()).clearVisited();
                ((Node) getChild2()).clearVisited();
            }
        }

        private boolean isVisited() {
            return isVisited;
        }

        private void visit() {
            isVisited = true;
        }

    }

    // helper class that defines the default post-order (leaves to root) traversal
    private class PostOrderIterator implements Iterator<GuideTreeNode<S, C>> {

        private Stack<Node> nodes;

        private PostOrderIterator() {
            getRoot().clearVisited();
            nodes = new Stack<Node>();
            nodes.push(getRoot());
        }

        // methods for Iterator

        @Override
        public boolean hasNext() {
            return !nodes.isEmpty();
        }

        @Override
        public GuideTreeNode<S, C> next() {
            while (hasNext()) {
                Node next = nodes.peek(), child1 = (Node) next.getChild1(), child2 = (Node) next.getChild2();
                if (child1 != null && !child1.isVisited()) {
                    nodes.push(child1);
                } else if (child2 != null && !child2.isVisited()) {
                    nodes.push(child2);
                } else {
                    next.visit();
                    return nodes.pop();
                }
            }
            return null;
        }

        @Override
        public void remove() {
            throw new UnsupportedOperationException();
        }

    }

}
