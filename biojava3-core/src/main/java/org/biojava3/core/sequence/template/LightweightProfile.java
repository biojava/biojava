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
 * Created on November 19, 2010
 * Author: Mark Chapman
 */

package org.biojava3.core.sequence.template;

import java.util.List;

/**
 * Defines a minimal data structure for reading and writing a sequence alignment.  The full {@code Profile} data
 * structure in the alignment module provides additional functionality.
 *
 * @author Mark Chapman
 * @param <S> each element of the alignment profile is of type S
 * @param <C> each element of an {@link Sequence} is a {@link Compound} of type C
 */
public interface LightweightProfile<S extends Sequence<C>, C extends Compound> {

    /**
     * List of output formats.
     */
    enum StringFormat {
        ALN,
        CLUSTALW,
        FASTA,
        GCG,
        MSF,
        PDBWEB
    }

    /**
     * Returns {@link Sequence} at given index.
     *
     * @param listIndex index of sequence in profile
     * @return desired sequence
     * @throws IndexOutOfBoundsException if listIndex < 1 or listIndex > number of sequences
     */
    S getAlignedSequence(int listIndex);

    /**
     * Returns a {@link List} containing the individual {@link Sequence}s of this alignment.
     *
     * @return list of aligned sequences
     */
    List<S> getAlignedSequences();

    /**
     * Returns the {@link Compound} elements of the original {@link Sequence}s at the given column.
     *
     * @param alignmentIndex column index within an alignment
     * @return the sequence elements
     * @throws IndexOutOfBoundsException if alignmentIndex < 1 or alignmentIndex > {@link #getLength()}
     */
    List<C> getCompoundsAt(int alignmentIndex);

    /**
     * Returns {@link CompoundSet} of all {@link Sequence}s
     *
     * @return set of {@link Compound}s in contained sequences
     */
    CompoundSet<C> getCompoundSet();

    /**
     * Returns the number of columns in the alignment profile.
     *
     * @return the number of columns
     */
    int getLength();

    /**
     * Returns the number of rows in this profile.  If any {@link Sequence}s are circular and overlap within the
     * alignment, the returned size will be greater than the number of sequences, otherwise the numbers will be equal.
     *
     * @return number of rows
     */
    int getSize();

    /**
     * Returns a simple view of the alignment profile.  This shows each sequence on a separate line (or multiple lines,
     * if circular) and nothing more.  This should result in {@link #getSize()} lines with {@link #getLength()}
     * {@link Compound}s per line.
     *
     * @return a simple view of the alignment profile
     */
    String toString();

    /**
     * Returns a formatted view of the alignment profile.  This shows the start and end indices of the profile for each
     * group of lines of the given width.  Each line may also be labeled.
     *
     * @param width limit on the line length
     * @return a formatted view of the alignment profile
     */
    String toString(int width);

    /**
     * Returns a formatted view of the alignment profile.  Details depend on the format given.
     *
     * @param format output format
     * @return a formatted view of the alignment profile
     */
    String toString(StringFormat format);

}
