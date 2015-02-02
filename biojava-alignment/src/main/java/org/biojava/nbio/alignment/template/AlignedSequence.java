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

import org.biojava.nbio.core.sequence.location.template.Location;
import org.biojava.nbio.core.sequence.location.template.Point;
import org.biojava.nbio.core.sequence.template.Compound;
import org.biojava.nbio.core.sequence.template.Sequence;

/**
 * Defines a data structure for a {@link Sequence} within an alignment.
 *
 * @author Mark Chapman
 * @param <C> each element of the {@link Sequence} is a {@link Compound} of type C
 */
public interface AlignedSequence<S extends Sequence<C>, C extends Compound> extends Sequence<C> {

    /**
     * Defines an alignment step in order to pass alignment information from an {@link Aligner} to a constructor.
     */
    enum Step { COMPOUND, GAP }

    /**
     * Nullifies cached arrays/objects.
     */
    void clearCache();

    /**
     * Returns the column index within an alignment corresponding to the given index in the original {@link Sequence}.
     * Both indices are 1-indexed and inclusive.
     *
     * @param sequenceIndex index in the original {@link Sequence}
     * @return column index within an alignment
     * @throws IndexOutOfBoundsException if sequenceIndex < 1 or sequenceIndex >
     *         {@link #getOriginalSequence()}.{@link #getLength()}
     */
    int getAlignmentIndexAt(int sequenceIndex);

    /**
     * Returns the {@link Point} within an alignment of the last element of the original {@link Sequence}.
     *
     * @return position within an alignment of final original {@link Sequence} element
     */
    Point getEnd();

    /**
     * Returns the {@link Location} of the original {@link Sequence} within an alignment.  This provides access to
     * additional substructure beyond start and end points.
     *
     * @return location within an alignment
     */
    Location getLocationInAlignment();

    /**
     * Returns number of gaps in the sequence.  This could be determined from the {@link Location} information or from
     * gap {@link Compound}s, which may not necessarily result in the same number.
     *
     * @return number of gaps in the sequence
     */
    int getNumGaps();

    /**
     * Returns the original {@link Sequence} before alignment.
     *
     * @return the original sequence
     */
    S getOriginalSequence();

    /**
     * Returns the maximum number of elements contributed to a column of an alignment by this {@link Sequence}.  If
     * this {@link Sequence} is circular, this number is >= 1.  If not, this overlap count is definitely 1.
     *
     * @return the most elements contributed to any alignment column
     */
    int getOverlapCount();

    /**
     * Returns the index in the original {@link Sequence} corresponding to the given index within an alignment.  Both
     * indices are 1-indexed and inclusive.
     *
     * @param alignmentIndex column index within an alignment
     * @return index in the original {@link Sequence}
     * @throws IndexOutOfBoundsException if alignmentIndex < 1 or alignmentIndex > {@link #getLength()}
     */
    int getSequenceIndexAt(int alignmentIndex);

    /**
     * Returns the {@link Point} within an alignment of the first element of the original {@link Sequence}.
     *
     * @return position within an alignment of first original {@link Sequence} element
     */
    Point getStart();

    /**
     * Returns true if this {@link Sequence} wraps around from the last alignment column back to the first.  This makes
     * overlap possible, but does not require an overlap count > 1.
     *
     * @return true for circular alignment elements
     */
    boolean isCircular();

    /**
     * Returns true if this {@link Sequence} has a gap at a particular alignment column.
     *
     * @param alignmentIndex column index within an alignment
     * @return true if this column has a gap
     * @throws IndexOutOfBoundsException if alignmentIndex < 1 or alignmentIndex > {@link #getLength()}
     */
    boolean isGap(int alignmentIndex);

}
