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
import org.biojava.nbio.core.sequence.template.Compound;
import org.biojava.nbio.core.sequence.template.CompoundSet;
import org.biojava.nbio.core.sequence.template.Sequence;

import java.util.List;

/**
 * Defines a data structure for the results of sequence alignment.  Every {@link List} returned is unmodifiable unless
 * the class implements the subinterface {@link MutableProfile}.
 *
 * @author Mark Chapman
 * @param <S> each element of the alignment {@link Profile} is of type S
 * @param <C> each element of an {@link AlignedSequence} is a {@link Compound} of type C
 */
public interface Profile<S extends Sequence<C>, C extends Compound> extends Iterable<AlignedSequence<S, C>> {

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
     * Returns {@link AlignedSequence} at given index.
     *
     * @param listIndex index of sequence in profile
     * @return desired sequence
     * @throws IndexOutOfBoundsException if listIndex < 1 or listIndex > number of sequences
     */
    AlignedSequence<S, C> getAlignedSequence(int listIndex);

    /**
     * Searches for the given {@link Sequence} within this alignment profile.  Returns the corresponding
     * {@link AlignedSequence}.
     *
     * @param sequence an original {@link Sequence}
     * @return the corresponding {@link AlignedSequence}
     */
    AlignedSequence<S, C> getAlignedSequence(S sequence);

    /**
     * Returns a {@link List} containing the individual {@link AlignedSequence}s of this alignment.
     *
     * @return list of aligned sequences
     */
    List<AlignedSequence<S, C>> getAlignedSequences();

    /**
     * Returns a {@link List} containing some of the individual {@link AlignedSequence}s of this alignment.
     *
     * @param listIndices indices of sequences in profile
     * @return list of aligned sequences
     */
    List<AlignedSequence<S, C>> getAlignedSequences(int... listIndices);

    /**
     * Returns a {@link List} containing some of the individual {@link AlignedSequence}s of this alignment.
     *
     * @param sequences original {@link Sequence}s
     * @return list of aligned sequences
     */
    List<AlignedSequence<S, C>> getAlignedSequences(S... sequences);

    /**
     * Returns the {@link Compound} at row of given sequence and column of alignment index.  If the given sequence has
     * overlap, this will return the {@link Compound} from the top row of the sequence.
     *
     * @param listIndex index of sequence in profile
     * @param alignmentIndex column index within an alignment
     * @return the sequence element
     * @throws IndexOutOfBoundsException if listIndex < 1, listIndex > number of sequences, alignmentIndex < 1, or
     *     alignmentIndex > {@link #getLength()}
     */
    C getCompoundAt(int listIndex, int alignmentIndex);

    /**
     * Returns the {@link Compound} at row of given sequence and column of alignment index.  If the given sequence has
     * overlap, this will return the {@link Compound} from the top row of the sequence.
     *
     * @param sequence either an {@link AlignedSequence} or an original {@link Sequence}
     * @param alignmentIndex column index within an alignment
     * @return the sequence element
     * @throws IndexOutOfBoundsException if alignmentIndex < 1 or alignmentIndex > {@link #getLength()}
     */
    C getCompoundAt(S sequence, int alignmentIndex);

    /**
     * Returns the number of each {@link Compound} in the given column for all compounds in {@link CompoundSet}.
     *
     * @param alignmentIndex column index within an alignment
     * @return list of counts
     * @throws IndexOutOfBoundsException if alignmentIndex < 1 or alignmentIndex > {@link #getLength()}
     */
    int[] getCompoundCountsAt(int alignmentIndex);

    /**
     * Returns the number of each {@link Compound} in the given column only for compounds in the given list.
     *
     * @param alignmentIndex column index within an alignment
     * @param compounds list of compounds to count
     * @return corresponding list of counts
     * @throws IndexOutOfBoundsException if alignmentIndex < 1 or alignmentIndex > {@link #getLength()}
     */
    int[] getCompoundCountsAt(int alignmentIndex, List<C> compounds);

    /**
     * Returns the {@link Compound} elements of the original {@link Sequence}s at the given column.
     *
     * @param alignmentIndex column index within an alignment
     * @return the sequence elements
     * @throws IndexOutOfBoundsException if alignmentIndex < 1 or alignmentIndex > {@link #getLength()}
     */
    List<C> getCompoundsAt(int alignmentIndex);

    /**
     * Returns {@link CompoundSet} of all {@link AlignedSequence}s
     *
     * @return set of {@link Compound}s in contained sequences
     */
    CompoundSet<C> getCompoundSet();

    /**
     * Returns the fraction of each {@link Compound} in the given column for all compounds in {@link CompoundSet}.
     *
     * @param alignmentIndex column index within an alignment
     * @return list of fractional weights
     * @throws IndexOutOfBoundsException if alignmentIndex < 1 or alignmentIndex > {@link #getLength()}
     */
    float[] getCompoundWeightsAt(int alignmentIndex);

    /**
     * Returns the fraction of each {@link Compound} in the given column only for compounds in the given list.
     *
     * @param alignmentIndex column index within an alignment
     * @param compounds list of compounds to count
     * @return corresponding list of fractional weights
     * @throws IndexOutOfBoundsException if alignmentIndex < 1 or alignmentIndex > {@link #getLength()}
     */
    float[] getCompoundWeightsAt(int alignmentIndex, List<C> compounds);

    /**
     * Returns the indices in the original {@link Sequence}s corresponding to the given column.  All indices are
     * 1-indexed and inclusive.
     *
     * @param alignmentIndex column index within an alignment
     * @return the sequence indices
     * @throws IndexOutOfBoundsException if alignmentIndex < 1 or alignmentIndex > {@link #getLength()}
     */
    int[] getIndicesAt(int alignmentIndex);

    /**
     * Searches for the given {@link Compound} within this alignment profile.  Returns column index nearest to the
     * start of the alignment profile, or -1 if not found.
     *
     * @param compound search element
     * @return index of column containing search element nearest to the start of the alignment profile
     */
    int getIndexOf(C compound);

    /**
     * Searches for the given {@link Compound} within this alignment profile.  Returns column index nearest to the end
     * of the alignment profile, or -1 if not found.
     *
     * @param compound search element
     * @return index of column containing search element nearest to the end of the alignment profile
     */
    int getLastIndexOf(C compound);

    /**
     * Returns the number of columns in the alignment profile.
     *
     * @return the number of columns
     */
    int getLength();

    /**
     * Returns a {@link List} containing the original {@link Sequence}s used for alignment.
     *
     * @return list of original sequences
     */
    List<S> getOriginalSequences();

    /**
     * Returns the number of rows in this profile.  If any {@link AlignedSequence}s are circular and overlap within the
     * alignment, the returned size will be greater than the number of sequences, otherwise the numbers will be equal.
     *
     * @return number of rows
     */
    int getSize();

    /**
     * Returns a {@link ProfileView} windowed to contain only the given {@link Location}.  This only includes the
     * {@link AlignedSequence}s which overlap the location.
     *
     * @param location portion of profile to view
     * @return a windowed view of the profile
     * @throws IllegalArgumentException if location is invalid
     */
    ProfileView<S, C> getSubProfile(Location location);

    /**
     * Returns true if any {@link AlignedSequence} has a gap at the given index.
     *
     * @param alignmentIndex column index within an alignment
     * @return true if any {@link AlignedSequence} has a gap at the given index
     * @throws IndexOutOfBoundsException if alignmentIndex < 1 or alignmentIndex > {@link #getLength()}
     */
    boolean hasGap(int alignmentIndex);

    /**
     * Returns true if any {@link AlignedSequence} is circular.  If so, sequences may simply wrap around from the end
     * to the start of the alignment or they may contribute multiple overlapping lines to the profile.
     *
     * @return true if any {@link AlignedSequence} is circular
     */
    boolean isCircular();

    /**
     * Returns a simple view of the alignment profile.  This shows each sequence on a separate line (or multiple lines,
     * if circular) and nothing more.  This should result in {@link #getSize()} lines with {@link #getLength()}
     * {@link Compound}s per line.
     *
     * @return a simple view of the alignment profile
     */
    @Override
	String toString();

    /**
     * Returns a formatted view of the alignment profile.  This shows the start and end indices of the profile and each
     * sequence for each group of lines of the given width.  Each line may also be labeled.
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
