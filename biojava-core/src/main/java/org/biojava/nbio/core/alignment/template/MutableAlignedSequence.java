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

import org.biojava.nbio.core.alignment.template.AlignedSequence;
import org.biojava.nbio.core.sequence.location.template.Location;
import org.biojava.nbio.core.sequence.template.Compound;
import org.biojava.nbio.core.sequence.template.Sequence;

/**
 * Defines a mutable (editable) data structure for an {@link AlignedSequence}.
 *
 * @author Mark Chapman
 * @author Paolo Pavan
 * @param <C> each element of the {@link AlignedSequence} is a {@link Compound} of type C
 */
public interface MutableAlignedSequence<S extends Sequence<C>, C extends Compound> extends AlignedSequence<S, C> {

    /**
     * Sets the position of the {@link AlignedSequence} to the given {@link Location} (start, gaps, end).
     *
     * @param location new location for this sequence
     * @throws IllegalArgumentException if location is invalid
     */
    void setLocationInAlignment(Location location);

    /**
     * Slides a part of the {@link AlignedSequence}.
     *
     * @param location portion of sequence moved in alignment coordinates
     * @param shift amount the alignment index changes for each contained element
     * @throws IllegalArgumentException if location is invalid or the shift causes a collision with stationary elements
     */
    void shiftAtAlignmentLocation(Location location, int shift);

    /**
     * Slides a part of the {@link AlignedSequence}.
     *
     * @param location portion of sequence moved in sequence coordinates
     * @param shift amount the alignment index changes for each contained element
     * @throws IllegalArgumentException if location is invalid or the shift causes a collision with stationary elements
     */
    void shiftAtSequenceLocation(Location location, int shift);

}
