/*
 *                  BioJava development code
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
 * Created on Jan. 22, 2016
 *
 */
package org.biojava.nbio.structure.io;

import java.util.List;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.io.mmcif.ChemCompGroupFactory;
import org.biojava.nbio.structure.io.mmcif.model.ChemComp;
import org.biojava.nbio.structure.io.mmcif.model.ChemCompAtom;

public class ChargeAdder {

	private Structure structure = null;

	public ChargeAdder(Structure structure) {
		this.structure = structure;
	}

	public void addCharges() {
		// Only adds charges to the first MODEL
		for(Chain c: structure.getChains()){
			for(Group g: c.getAtomGroups()){
				ChemComp thisChemComp = ChemCompGroupFactory.getChemComp(g.getPDBName());
				List<ChemCompAtom> chemAtoms = thisChemComp.getAtoms();
				List<Atom> protAtoms = g.getAtoms();
				for(int i=0; i<protAtoms.size();i++){
					ChemCompAtom cca = chemAtoms.get(i);
					Atom a = protAtoms.get(i);
					// Get the charge 
					short charge = Short.parseShort(cca.getCharge());
					a.setCharge(charge);
				}
			}
		}
	}
}
