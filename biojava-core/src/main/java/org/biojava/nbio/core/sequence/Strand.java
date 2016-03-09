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
package org.biojava.nbio.core.sequence;


/**
 * Provides a way of representing the strand of a sequence, location
 * hit or feature.
 *
 * @author ayates
 */
public enum Strand {

	POSITIVE("+", 1), NEGATIVE("-", -1), UNDEFINED(".", 0);
	private final String stringRepresentation;
	private final int numericRepresentation;

	private Strand(String stringRepresentation, int numericRepresentation) {
		this.stringRepresentation = stringRepresentation;
		this.numericRepresentation = numericRepresentation;
	}

	public int getNumericRepresentation() {
		return numericRepresentation;
	}

	public String getStringRepresentation() {
		return stringRepresentation;
	}

	public Strand getReverse() {
		switch (this) {
			case POSITIVE:
				return NEGATIVE;
			case NEGATIVE:
				return POSITIVE;
			default:
				return UNDEFINED;
		}
	}
}
