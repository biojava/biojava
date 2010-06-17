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
 * Defines an {@link Aligner} which builds a score matrix during computation.
 *
 * @author Mark Chapman
 * @param <S> each element of the alignment {@link Profile} is of type S
 * @param <C> each element of an {@link AlignedSequence} is a {@link Compound} of type C
 */
public interface MatrixAligner<S extends Sequence<C>, C extends Compound> extends Aligner<S, C> {

    /**
     * Returns the entire score matrix built during alignment.  The first dimension has the length of the first (query)
     * {@link Sequence}; the second has the length of the second (target) {@link Sequence}.
     *
     * @return the score matrix
     */
    short[][] getScoreMatrix();

    /**
     * Returns a single value from within the score matrix.
     *
     * @param queryIndex index in the first (query) {@link Sequence}
     * @param targetIndex index in the second (target) {@link Sequence}
     * @return score at given point in score matrix
     * @throws IndexOutOfBoundsException if queryIndex < 0, queryIndex > query length, targetIndex < 0, or
     *     targetIndex > target length
     */
    short getScoreMatrixAt(int queryIndex, int targetIndex);

    /**
     * Returns a depiction of the score matrix as a {@link String}.  This may include additional description such as
     * labeling each axes with the {@link Compound}s of each {@link Sequence}.
     *
     * @return the score matrix as a character sequence
     */
    String getScoreMatrixAsString();

}
