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
package org.biojava.nbio.core.sequence.template;

import java.util.List;
import java.util.Set;


public interface CompoundSet<C extends Compound> {

	/**
	 * Returns the maximum size of a compound String this set holds
	 */
    int getMaxSingleCompoundStringLength();

	/**
	 * Returns true if all String representations of Compounds are of the
	 * same length.
	 */
    boolean isCompoundStringLengthEqual();

	/**
	 * Return null if not recognised. Throw IllegalArgumentException if string
	 * is longer than maximum allowed by {@link #getStringForCompound(Compound)}.
	 */
    C getCompoundForString(String string);

	String getStringForCompound(C compound);

	boolean compoundsEquivalent(C compoundOne, C compoundTwo);

	boolean isValidSequence(Sequence<C> sequence);

	Set<C> getEquivalentCompounds(C compound);

	boolean hasCompound(C compound);

	List<C> getAllCompounds();

	boolean isComplementable();
}
