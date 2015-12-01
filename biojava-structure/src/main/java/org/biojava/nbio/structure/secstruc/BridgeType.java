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
 * A bridge is formed by two non-overlapping stretches of three residues each
 * (i-1,i,i+1) and (j-1,j,j+1), where i<j.
 * <p>
 * Depending on two basic patterns, a Bridge can be either of type parallel (H
 * bonds in {(i-1,j) and (j,i+1)} OR {(j-1,i) and (i,j-1)}) or antiparallel (H
 * bonds in {(i1,j) and (j,i} OR {(i-1,j+1) and (j-1,i+1)})
 * 
 * @author Andreas Prlic
 * @author Aleix Lafita
 *
 */
public enum BridgeType {

	parallel("parallel", 'p'), 
	antiparallel("antiparallel", 'a');

	public final Character type;
	public final String name;

	private BridgeType(String name, Character stype) {
		this.name = name;
		this.type = stype;
	}

	public static BridgeType fromCharacter(Character stype) {

		for (BridgeType c : BridgeType.values()) {
			if (c.type.equals(stype))
				return c;
		}
		return null;
	}

	@Override
	public String toString() {
		return type.toString();
	}

}
