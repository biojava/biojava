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
 * Container that represents a beta Bridge between two residues. It contains the
 * two partner indices and the type of the bridge. For consistency, partner1 is
 * always the small index.
 * 
 * @author Aleix Lafita
 *
 */
public class BetaBridge {

	BridgeType type;
	int partner1;
	int partner2;

	public BetaBridge(int i, int j, BridgeType t) {
		partner1 = Math.min(i, j);
		partner2 = Math.max(i, j);
		type = t;
	}

	@Override
	public boolean equals(Object o) {

		if (!(o instanceof BetaBridge))
			return false;

		BetaBridge b = (BetaBridge) o;
		if (type != b.type)
			return false;
		if (partner1 != b.partner1)
			return false;
		if (partner2 != b.partner2)
			return false;
		return true;
	}

	public BridgeType getType() {
		return type;
	}

	public int getPartner1() {
		return partner1;
	}

	public int getPartner2() {
		return partner2;
	}

}
