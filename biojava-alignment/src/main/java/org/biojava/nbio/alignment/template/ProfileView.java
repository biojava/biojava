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

/**
 * Defines a data structure for a view of sequence alignment.
 *
 * @author Mark Chapman
 * @param <S> each element of the alignment {@link Profile} is of type S
 * @param <C> each element of an {@link AlignedSequence} is a {@link Compound} of type C
 */
public interface ProfileView<S extends Sequence<C>, C extends Compound> extends Profile<S, C> {

    /**
     * Returns the column index of the viewed {@link Profile} corresponding to the final element in this view
     *
     * @return column index of this view's final element
     */
    int getEnd();

    /**
     * Returns the column index of the viewed {@link Profile} corresponding to the first element in this view
     *
     * @return column index of this view's first element
     */
    int getStart();

    /**
     * Returns the entire {@link Profile} being viewed
     *
     * @return the entire alignment profile
     */
    Profile<S, C> getViewedProfile();

}
