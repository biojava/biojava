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

import org.biojava.nbio.core.alignment.template.ProfilePair;
import org.biojava.nbio.core.sequence.template.Compound;
import org.biojava.nbio.core.sequence.template.Sequence;

/**
 * Defines an {@link Aligner} for a pair of {@link Profile}s.
 *
 * @author Mark Chapman
 * @param <S> each {@link Sequence} in the pair of alignment {@link Profile}s is of type S
 * @param <C> each element of an {@link AlignedSequence} is a {@link Compound} of type C
 */
public interface ProfileProfileAligner<S extends Sequence<C>, C extends Compound> extends Aligner<S, C>,
        ProfileProfileScorer<S, C> {

    /**
     * Returns alignment profile pair.
     *
     * @return alignment profile pair
     */
    ProfilePair<S, C> getPair();

}
