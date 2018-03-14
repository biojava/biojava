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
 * Created on Sep 12, 2007
 *
 */
package org.biojava.nbio.structure.io;

import org.biojava.nbio.structure.*;

import java.util.ArrayList;
import java.util.List;

/** Converts full atom representations to Calpha only ones.
 *
 * @author Andreas Prlic
 * @version %I% %G%
 */
public class CAConverter {


	/** Convert a List of chain objects to another List of chains, containing Representative atoms only.
	 *
	 * @param chains list of chains
	 * @return a list of chains
	 * @since Biojava 4.1.0
	 */
	public static List<Chain> getRepresentativeAtomsOnly(List<Chain> chains){
		List<Chain> newChains = new ArrayList<Chain>();

		for (Chain chain : chains){
			Chain newChain = getRepresentativeAtomsOnly(chain);
			newChains.add(newChain);
		}

		return newChains;
	}

	/** 
	 * Convert a Chain to a new Chain containing C-alpha atoms only.
	 *
	 * @param chain to convert
	 * @return a new chain containing Amino acids with C-alpha only.
	 * @since Biojava 4.1.0
	 */
	public static Chain getRepresentativeAtomsOnly(Chain chain){

		Chain newChain = new ChainImpl();
		newChain.setId(chain.getId());
		newChain.setName(chain.getName());
		newChain.setEntityInfo(chain.getEntityInfo());
		newChain.setSwissprotId(chain.getSwissprotId());

		List<Group> groups = chain.getAtomGroups();

		grouploop:
			for (Group g: groups){
				List<Atom> atoms = g.getAtoms();

				if ( ! (g instanceof AminoAcid))
					continue;

				for (Atom a : atoms){

					if ( a.getName().equals(StructureTools.CA_ATOM_NAME) && a.getElement()==Element.C){
						// we got a CA atom in this group!
						AminoAcid n = new AminoAcidImpl();
						n.setPDBName(g.getPDBName());
						n.setResidueNumber(g.getResidueNumber());
						n.addAtom(a);
						newChain.addGroup(n);
						continue grouploop;

					}
				}

			}
		return newChain;
	}
}
