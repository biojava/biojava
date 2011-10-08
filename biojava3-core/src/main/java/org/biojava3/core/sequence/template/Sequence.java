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
 * Created on 01-21-2010
 *
 * @author Richard Holland
 *
 *
 */
package org.biojava3.core.sequence.template;

import java.util.List;

/**
 * Main interface for defining a collection of Compounds and accessing them
 * using biological indexes
 *
 * @author Richard Holland
 * @author Andy Yates
 * @author Scooter Willis
 * @param <C> Compound a Sequence holds
 */
public interface Sequence<C extends Compound> extends Iterable<C>, Accessioned {

    /**
     * Returns the length of the Sequence
     */
    public int getLength();

    /**
     * Returns the Compound at the given biological index
     *
     * @param position Biological index (1 to n)
     * @return Compound at the specified position
     */
    public C getCompoundAt(int position);

    /**
     * Scans through the Sequence looking for the first occurrence of the given
     * compound
     *
     * @param compound Compounds to look for
     * @return Index of the first position of the compound in the sequence (1 to n)
     */
    public int getIndexOf(C compound);

/**
     * Scans through the Sequence looking for the last occurrence of the given
     * compound
     *
     * @param compound Compounds to look for
     * @return Index of the last position of the compound in the sequence (1 to n)
     */
    public int getLastIndexOf(C compound);

    /**
     * Returns the String representation of the Sequence
     */
    public String getSequenceAsString();

    /**
     * Returns the Sequence as a List of compounds
     */
    public List<C> getAsList();

    /**
     * Returns a portion of the sequence from the different positions. This is
     * indexed from 1
     *
     * @param start Biological index start; must be greater than 0
     * @param end Biological end; must be less than length + 1
     * @return A SequenceView of the offset
     */
    public SequenceView<C> getSubSequence(Integer start, Integer end);

    /**
     * Gets the compound set used to back this Sequence
     */
    public CompoundSet<C> getCompoundSet();

    /**
     * Returns the number of times we found a compound in the Sequence
     *
     * @param compounds Vargs of the compounds to count
     * @return Number of times a compound was found
     */
    public int countCompounds(C... compounds);

    /**
     * Does the <em>right thing</em> to get the inverse of the current
     * Sequence. This means either reversing the Sequence and optionally
     * complementing the Sequence.
     */
    public SequenceView<C> getInverse();
}
