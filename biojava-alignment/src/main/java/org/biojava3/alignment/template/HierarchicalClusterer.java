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

package org.biojava3.alignment.template;

import org.biojava3.core.sequence.template.Compound;
import org.biojava3.core.sequence.template.Sequence;

/**
 * Defines a clustering algorithm that converts a distance matrix into a tree.
 *
 * @author Mark Chapman
 * @param <S> each {@link Sequence} in the tree is of type S
 * @param <C> each element of a {@link Sequence} is a {@link Compound} of type C
 */
public interface HierarchicalClusterer<S extends Sequence<C>, C extends Compound> {

    /**
     * Returns the distance matrix used in clustering.  May be calculated from another original source.
     *
     * @return the distance matrix input to clustering
     */
    float[][] getDistanceMatrix();

    /**
     * Returns the root node of the tree resulting from this clustering algorithm.
     *
     * @return the resulting tree output from clustering
     */
    GuideTreeNode<S, C> getRoot();

}
