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

import org.biojava.nbio.core.alignment.template.ProfilePair;
import org.biojava.nbio.core.alignment.template.Profile;
import org.biojava.nbio.core.sequence.template.Compound;
import org.biojava.nbio.core.sequence.template.Sequence;

/**
 * Defines a mutable (editable) data structure for a {@link ProfilePair}.
 *
 * @author Mark Chapman
 * @author Paolo Pavan
 * @param <S> each element of the alignment {@link Profile} is of type S
 * @param <C> each element of an {@link AlignedSequence} is a {@link Compound} of type C
 */
public interface MutableProfilePair<S extends Sequence<C>, C extends Compound> extends MutableProfile<S, C>,
        ProfilePair<S, C> {

    /**
     * Sets both {@link Profile}s of the pair.
     *
     * @param query becomes the first {@link Profile} of the pair
     * @param target becomes the second {@link Profile} of the pair
     * @throws IllegalArgumentException if query and target are different lengths
     */
    void setPair(Profile<S, C> query, Profile<S, C> target);

    /**
     * Sets the first {@link Profile} of the pair.
     *
     * @param query becomes the first {@link Profile} of the pair
     * @throws IllegalArgumentException if (new) query and (old) target are different lengths
     */
    void setQuery(Profile<S, C> query);

    /**
     * Sets the second {@link Profile} of the pair.
     *
     * @param target becomes the second {@link Profile} of the pair
     * @throws IllegalArgumentException if (old) query and (new) target are different lengths
     */
    void setTarget(Profile<S, C> target);

}
