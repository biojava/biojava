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
 * Created on 2012-11-20
 *
 */

package org.biojava.bio.structure.align.symm.protodomain;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Identifier;
import org.biojava.bio.structure.StructureIdentifier;
import org.biojava.bio.structure.align.model.AFPChain;

/**
 * Contains the information required to uniquely identify an alignment and information resulting from the alignment.
 * Although the implementation will not enforce it, if one aligned structure is a query being aligned against multiple other structures, it should be assigned as position 1.
 * @author dmyerstu
 *
 */
public class AFPChainAndAtoms {

	private final AFPChain afpChain;
	private final Atom[] ca1;
	private final Atom[] ca2;

	public String getName1() {
		return afpChain.getName1();
	}
	
	public StructureIdentifier getIdentifier1() {
		return Identifier.loadIdentifier(afpChain.getName1(), ResourceList.get().getCache());
	}

	public StructureIdentifier getIdentifier2() {
		return Identifier.loadIdentifier(afpChain.getName2(), ResourceList.get().getCache());
	}

	public String getName2() {
		return afpChain.getName2();
	}

	public AFPChainAndAtoms(AFPChain afpChain, Atom[] ca1, Atom[] ca2) {
		this.afpChain = afpChain;
		this.ca1 = ca1;
		this.ca2 = ca2;
	}

	public AFPChain getAfpChain() {
		return afpChain;
	}

	public Atom[] getCa1() {
		return ca1;
	}

	public Atom[] getCa2() {
		return ca2;
	}

}
