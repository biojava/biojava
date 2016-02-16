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

/**
 * This enum contains all of the secondary structure types found in the DSSP
 * output. It also contains some methods to operate with the SS types.
 * <p>
 * When compared, the types are sorted in the declaration order of the enum,
 * which is the DSSP preference of type assignment.
 * 
 * @author Andreas Prlic
 * @author Aleix Lafita
 *
 */
public enum SecStrucType {

	helix4("alpha Helix", 'H'), 
	extended("Extended", 'E'), 
	bridge("Bridge", 'B'), 
	helix3("3-10 Helix", 'G'), 
	helix5("pi Helix", 'I'), 
	turn("Turn", 'T'), 
	bend("Bend", 'S'), 
	coil("Coil", ' ');

	public final Character type;
	public final String name;

	private SecStrucType(String name, Character stype) {
		this.name = name;
		this.type = stype;
	}

	/**
	 * Converts a Character representing a Secondary Structure type into the
	 * corresponding enum object.
	 * 
	 * @param stype
	 *            the character representing the SS type
	 * @return SecStrucType or null if the character is invalid
	 */
	public static SecStrucType fromCharacter(Character stype) {

		for (SecStrucType c : SecStrucType.values()) {
			if (c.type.equals(stype)) {
				return c;
			}
		}
		return null;
	}

	@Override
	public String toString() {
		return type.toString();
	}

	/**
	 * Helix type can be 3-10 helix, pi-helix or alpha-helix.
	 * 
	 * @return true if the type is any of the helix types, false otherwise
	 */
	public boolean isHelixType() {
		if (type.equals(helix4.type) || type.equals(helix3.type)
				|| type.equals(helix5.type))
			return true;
		else
			return false;
	}

	/**
	 * A Beta-Strand is an extended set of sequential Bridges that, together
	 * with other Beta-Strands, is part of a Beta-Sheet.
	 * 
	 * @return true if the type is a Beta-Strand
	 */
	public boolean isBetaStrand() {
		if (type.equals(extended.type))
			return true;
		else
			return false;
	}

}