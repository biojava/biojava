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

import org.biojava.nbio.core.sequence.template.Compound;
import org.biojava.nbio.core.sequence.template.Sequence;

/**
 * Defines a data structure for the results of pairwise sequence alignment.
 *
 * @author Mark Chapman
 * @author Paolo Pavan
 * @param <S> each element of the alignment {@link Profile} is of type S
 * @param <C> each element of an {@link AlignedSequence} is a {@link Compound} of type C
 */
public interface SequencePair<S extends Sequence<C>, C extends Compound> extends Profile<S, C> {

    /**
     * Returns the {@link Compound} in query sequence at given column index in alignment.
     *
     * @param alignmentIndex column index in alignment
     * @return the query sequence element
     * @throws IndexOutOfBoundsException if alignmentIndex < 1 or alignmentIndex > {@link #getLength()}
     */
    C getCompoundInQueryAt(int alignmentIndex);

    /**
     * Returns the {@link Compound} in target sequence at given column index in alignment.
     *
     * @param alignmentIndex column index in alignment
     * @return the target sequence element
     * @throws IndexOutOfBoundsException if alignmentIndex < 1 or alignmentIndex > {@link #getLength()}
     */
    C getCompoundInTargetAt(int alignmentIndex);

    /**
     * Returns the query index corresponding to a given alignment column.
     *
     * @param alignmentIndex column index in alignment
     * @return index in query {@link Sequence}
     * @throws IndexOutOfBoundsException if alignmentIndex < 1 or alignmentIndex > {@link #getLength()}
     */
    int getIndexInQueryAt(int alignmentIndex);

    /**
     * Returns the query index corresponding to a given target index.
     *
     * @param targetIndex index in target {@link Sequence}
     * @return index in query {@link Sequence}
     * @throws IndexOutOfBoundsException if targetIndex < 1 or targetIndex > {@link #getTarget()}.getLength()
     */
    int getIndexInQueryForTargetAt(int targetIndex);

    /**
     * Returns the target index corresponding to a given alignment column.
     *
     * @param alignmentIndex column index in alignment
     * @return index in target {@link Sequence}
     * @throws IndexOutOfBoundsException if alignmentIndex < 1 or alignmentIndex > {@link #getLength()}
     */
    int getIndexInTargetAt(int alignmentIndex);

    /**
     * Returns the target index corresponding to a given query index.
     *
     * @param queryIndex index in query {@link Sequence}
     * @return index in target {@link Sequence}
     * @throws IndexOutOfBoundsException if queryIndex < 1 or queryIndex > {@link #getQuery()}.getLength()
     */
    int getIndexInTargetForQueryAt(int queryIndex);

    /**
     * Returns the number of indices for which both the query and target sequences have an identical {@link Compound}.
     *
     * @return the number of identical indices
     */
    int getNumIdenticals();

    /**
     * Returns the number of indices for which both the query and target sequences have a similar {@link Compound}.
     *
     * @return the number of similar indices
     */
    int getNumSimilars();

    /**
     * Returns the first {@link AlignedSequence} of the pair.
     *
     * @return the first {@link AlignedSequence} of the pair
     */
    AlignedSequence<S, C> getQuery();

    /**
     * Returns the second {@link AlignedSequence} of the pair.
     *
     * @return the second {@link AlignedSequence} of the pair
     */
    AlignedSequence<S, C> getTarget();

}
