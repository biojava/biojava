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

package org.biojava.nbio.core.alignment.template;

import org.biojava.nbio.core.alignment.template.SequencePair;
import org.biojava.nbio.core.alignment.template.AlignedSequence;
import org.biojava.nbio.core.sequence.template.Compound;
import org.biojava.nbio.core.sequence.template.Sequence;

/**
 * Defines a mutable (editable) data structure for the results of pairwise sequence alignment.
 *
 * @author Mark Chapman
 * @author Paolo Pavan
 * @param <S> each element of the alignment {@link Profile} is of type S
 * @param <C> each element of an {@link AlignedSequence} is a {@link Compound} of type C
 */
public interface MutableSequencePair<S extends Sequence<C>, C extends Compound> extends MutableProfile<S, C>,
        SequencePair<S, C> {

    /**
     * Sets both {@link AlignedSequence}s of the pair.
     *
     * @param query becomes the first {@link AlignedSequence} of the pair
     * @param target becomes the second {@link AlignedSequence} of the pair
     * @throws IllegalArgumentException if query and target are different lengths
     */
    void setPair(AlignedSequence<S, C> query, AlignedSequence<S, C> target);

    /**
     * Sets the first {@link AlignedSequence} of the pair.
     *
     * @param query becomes the first {@link AlignedSequence} of the pair
     * @throws IllegalArgumentException if (new) query and (old) target are different lengths
     */
    void setQuery(AlignedSequence<S, C> query);

    /**
     * Sets the second {@link AlignedSequence} of the pair.
     *
     * @param target becomes the second {@link AlignedSequence} of the pair
     * @throws IllegalArgumentException if (old) query and (new) target are different lengths
     */
    void setTarget(AlignedSequence<S, C> target);

}
