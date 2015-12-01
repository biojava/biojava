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
 */
package org.biojava.nbio.structure.secstruc;

import org.biojava.nbio.structure.ResidueNumber;
import org.biojava.nbio.structure.ResidueRangeAndLength;

/**
 * A secondary structure element (SSE) is an object representing a block of
 * sequential residues that share the same secondary structure type.
 * 
 * @author Aleix Lafita
 * @since 4.1.1
 *
 */
public class SecStrucElement {

	private SecStrucType type;
	private ResidueRangeAndLength range;
	private int index;

	/**
	 * Create a new SSE object. The start and end residue numbers cannot be the
	 * same.
	 * 
	 * @param type
	 *            object describing the type of SS
	 * @param start
	 *            first residue of the SSE
	 * @param end
	 *            final residue of the SSE
	 * @param length
	 *            number of residues included in the SSE
	 * @param index
	 * @param chainID
	 *            the chain ID
	 */
	public SecStrucElement(SecStrucType type, ResidueNumber start,
			ResidueNumber end, int length, int index, String chainID) {

		this.type = type;
		this.index = index;
		range = new ResidueRangeAndLength(chainID, start, end, length);
	}

	/**
	 * Returns the {@link SecStrucType} of this element.
	 * 
	 * @return
	 */
	public SecStrucType getType() {
		return type;
	}

	/**
	 * Returns the index of the SSE for its type. This is, the sequential
	 * position of this SSE relative to the other SSE of the same type.
	 * 
	 * @return
	 */
	public int getIndex() {
		return index;
	}
	
	/**
	 * Return the length (number of residues) in the SSE.
	 * 
	 * @return
	 */
	public int getLength() {
		return range.getLength();
	}

	/**
	 * Returns the ID of this element. The ID is the concatenation of the type
	 * letter and the numerical element identifier (e.g. H1, S1, ...).
	 * 
	 * @return
	 */
	public String getId() {
		return type.toString() + index + "";
	}

	/**
	 * Returns the residue range of this SSE.
	 * 
	 * @return
	 */
	public ResidueRangeAndLength getRange() {
		return range;
	}

	@Override
	public String toString() {
		return getId() + ": " + range.toString();
	}

}
