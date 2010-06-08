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
 * Created on Jun 6, 2010
 * Author: Jianjiong Gao 
 *
 */

package org.biojava3.protmod.parser;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.Structure;

import org.biojava3.protmod.Component;
import org.biojava3.protmod.ModificationCategory;
import org.biojava3.protmod.ModificationCondition;
import org.biojava3.protmod.ModifiedCompound;
import org.biojava3.protmod.ModifiedCompoundFactory;
import org.biojava3.protmod.ProteinModification;

/**
 * Identify modified residues in a 3-D structure.
 * 
 * @author Jianjiong Gao
 * @since 3.0
 */
public class ModifiedResidueParser implements ProteinModificationParser {
	/**
	 * Parse chemically modified residues in a structure.
	 * @param structure query {@link Structure}.
	 * @param potentialModifications query {@link ProteinModification}s.
	 * @param modelnr model number.
	 * @return an list of {@link ModifiedCompound}s, or null if the
	 *  nodelnr is larger than the number of models in the structure.
	 * @throws IllegalArgumentException if null structure, or null or 
	 *  empty potentialModifications, or potentialModifications contain 
	 *  modifications other than ModifiedResidues.
	 */
	@Override
	public List<ModifiedCompound> parse(final Structure structure, 
			final Set<ProteinModification> potentialModifications,
			final int modelnr) {
		if (structure==null) {
			throw new IllegalArgumentException("Null structure.");
		}
		
		if (potentialModifications==null || potentialModifications.isEmpty()) {
			throw new IllegalArgumentException("Null or empty potentialModifications.");
		}
		
		for (ProteinModification mod:potentialModifications) {
			if (mod.getCategory()!=ModificationCategory.CHEMICAL_MODIFICATION) {
				throw new IllegalArgumentException("Only CHEMICAL_MODIFICATION is allowed.");
			}
		}
		
		if (modelnr >= structure.nrModels())
			return null;
		
		List<ModifiedCompound> ret = new ArrayList<ModifiedCompound>();
		
		// TODO: how to deal with multi-model structure?
		List<Chain> chains = structure.getChains(modelnr);
		for (Chain chain : chains) {
			List<Group> residues = chain.getSeqResGroups();
			int sizeRes = residues.size();
			
			// for all amino acid
			for (int iRes=0; iRes<sizeRes; iRes++) {
				Group residue = residues.get(iRes);
				String pdbccId = residue.getPDBName();
				ProteinModification mod = ProteinModification.getByPdbccId(pdbccId);
				
				if (mod==null || !potentialModifications.contains(mod)) {
					continue;
				}
				
				ModificationCondition condition = mod.getCondition();
				Component comp = condition.getComponents().get(0);
				
				// TODO: is this the correct way to determine N/C-terminal?
				if ((comp.isNTerminal() && iRes==0) ||           // N-terminal
						(comp.isCTerminal() && iRes==sizeRes-1)) {    // C-terminal
					continue;
				}
				
				ModifiedCompound modRes = ModifiedCompoundFactory
					.createModifiedResidue(mod, residue);
				ret.add(modRes);
			}
		}
		
		return ret;
	}
}
