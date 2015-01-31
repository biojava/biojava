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
 * Created on June 7, 2010
 * Author: Mark Chapman
 */

package org.biojava.nbio.alignment.template;

import org.biojava.nbio.core.sequence.template.Compound;
import org.biojava.nbio.core.sequence.template.Sequence;

import javax.swing.tree.TreeNode;
import java.util.concurrent.Future;

/**
 * Defines a data structure for the node in a guide tree used during progressive multiple sequence alignment.
 * 
 * @author Mark Chapman
 * @param <S> each {@link Sequence} in the tree is of type S
 * @param <C> each element of a {@link Sequence} is a {@link Compound} of type C
 */
public interface GuideTreeNode<S extends Sequence<C>, C extends Compound> extends TreeNode {

    /**
     * Returns the first child node of this node.  For leaf nodes (sequences), this will be null.
     *
     * @return the first child node of this node
     */
    GuideTreeNode<S, C> getChild1();

    /**
     * Returns the second child node of this node.  For leaf nodes (sequences), this will be null.
     *
     * @return the second child node of this node
     */
    GuideTreeNode<S, C> getChild2();

    /**
     * Returns the difference in height of this node and it's parent node.  A likely meaning of this distance is half
     * the percent difference between this node and it's sibling node.
     *
     * @return the difference in height of this node to it's parent node
     */
    double getDistanceToParent();

    /**
     * Returns the name of this node.  For leaf nodes (sequences), this will likely be the accession ID.
     *
     * @return the name of this node
     */
    String getName();

    /**
     * Returns the profile stored at this node.  If the node is a leaf, the profile is that of a single sequence.  If
     * not, this returns null until {@link #setProfile(Profile)} has been called.
     *
     * @return the profile stored at this node
     */
    Profile<S, C> getProfile();

    /**
     * Returns the profile future stored at this node, but does not force the calculation, yet.  This allows alignment
     * tasks for the entire tree to be queued in a post-order traversal before concurrent execution.
     *
     * @return the profile future stored at this node
     */
    Future<ProfilePair<S, C>> getProfileFuture();

    /**
     * Stores the given profile.
     *
     * @param profile new profile stored at this node
     */
    void setProfile(Profile<S, C> profile);

    /**
     * Stores the given profile future.  This allows concurrent execution of alignment tasks.
     *
     * @param profileFuture new profile to be calculated and then stored at this node
     */
    void setProfileFuture(Future<ProfilePair<S, C>> profileFuture);

}
