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
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * A class to add appropriate charge information to a structure.
 * @author Anthony Bradley
 *
 */
public class ChargeAdder {

	private static final Logger logger = LoggerFactory.getLogger(ChargeAdder.class);

	/**
	 * Function to add the charges to a given structure.
	 */
	public static void addCharges(Structure structure) {
		// Loop through the models
		for(int i=0; i<structure.nrModels(); i++){
			for(Chain c: structure.getChains(i)){
				for(Group g: c.getAtomGroups()){
					ChemComp thisChemComp = ChemCompGroupFactory.getChemComp(g.getPDBName());
					List<ChemCompAtom> chemAtoms = thisChemComp.getAtoms();
					for(ChemCompAtom chemCompAtom : chemAtoms) {
						Atom atom = g.getAtom(chemCompAtom.getAtom_id());	
						String stringCharge = chemCompAtom.getCharge();
						short shortCharge = 0;
						if (stringCharge!=null){
							if(!stringCharge.equals("?")){
								try{
									shortCharge = Short.parseShort(stringCharge);
								}
								catch(NumberFormatException e){
									logger.warn("Number format exception. Parsing '"+stringCharge+"' to short");
								}
							}
							else{
								logger.warn("? charge on atom "+chemCompAtom.getAtom_id()+" in group "+thisChemComp.getId());
							}
						}
						else{
							logger.warn("Null charge on atom "+chemCompAtom.getAtom_id()+" in group "+thisChemComp.getId());
						}
						if(atom!=null){
							atom.setCharge(shortCharge);
						}
						// Now do the same for alt locs
						for (Group altLoc : g.getAltLocs()) {
							Atom altAtom = altLoc.getAtom(chemCompAtom.getAtom_id());
							if(altAtom!=null){
								altAtom.setCharge(shortCharge);
							}
						}
					}
				}

			}
		}
	}
}
