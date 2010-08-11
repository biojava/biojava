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
     * sequence + 1; the second has the length of the second (target) sequence + 1; the third has length equal to the
     * number of scores stored per pairing of an element from each {@link Sequence}.
     *
     * @return the score matrix
     */
    short[][][] getScoreMatrix();

    /**
     * Returns a depiction of the score matrix as a {@link String}.  This may include additional description such as
     * labels for each dimension: element from query sequence, element from target sequence, and meaning of each score.
     *
     * @return the score matrix as a character sequence
     */
    String getScoreMatrixAsString();

}
