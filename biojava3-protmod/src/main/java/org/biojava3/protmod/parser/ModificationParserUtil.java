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
 * Created on Jun 11, 2010
 * Author: Jianjiong Gao 
 *
 */

package org.biojava3.protmod.parser;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.StructureException;

import org.biojava3.protmod.Component;
import org.biojava3.protmod.ModificationCondition;
import org.biojava3.protmod.ProteinModification;

final class ModificationParserUtil {
	
	private ModificationParserUtil() {};
	
	/**
	 * 
	 * @param chain {@link Chain}.
	 * @param modifications a set of {@link ProteinModification}s.
	 * @return map from component to list of corresponding residues
	 *  in the chain.
	 */
	static Map<Component, List<Group>> getModifiableResidues(
			final Chain chain, 
			final Set<ProteinModification> modifications) {
		if (chain==null || modifications==null) {
			throw new IllegalArgumentException("Null argument(s).");
		}
		
		Set<Component> comps = new HashSet<Component>();
		for (ProteinModification mod : modifications) {
			ModificationCondition condition = mod.getCondition();
			for (Component comp : condition.getComponents()) {
				comps.add(comp);
			}
		}
		
		Map<Component, List<Group>> mapCompRes = 
			new HashMap<Component, List<Group>>();
		
		List<Group> residues = chain.getSeqResGroups();
		
		if (residues==null || residues.isEmpty()) {
			return mapCompRes;
		}
		
		// for all residue
		for (Group res : residues) {
			String pdbccId = res.getPDBName();
			Component comp = Component.of(pdbccId);
			if (!comps.contains(comp)) {
				continue;
			}
			List<Group> groups = mapCompRes.get(comp);
			if (groups==null) {
				groups = new ArrayList<Group>();
				mapCompRes.put(comp, groups);
			}
			groups.add(res);
		}
		
		// for N-terminal
		Group res = residues.get(0);
		Component comp = Component.of(res.getPDBName(), true, false);
		if (comps.contains(comp)) {
			List<Group> groups = new ArrayList<Group>(1);
			groups.add(res);
			mapCompRes.put(comp, groups);
		}
		
		// for C-terminal
		res = residues.get(residues.size()-1);
		comp = Component.of(res.getPDBName(), false, true);
		if (comps.contains(comp)) {
			List<Group> groups = new ArrayList<Group>(1);
			groups.add(res);
			mapCompRes.put(comp, groups);
		}

		return mapCompRes;
	}
	
	/**
	 * Find the nearest Atoms between a pair of {@link Group}s.
	 * @param group1
	 * @param group2
	 * @return a pair of Atoms.
	 * @throws StructureException ...
	 */
	static Atom[] findNearestAtoms(Group group1, Group group2)
			throws StructureException {		
		double nearestDistance = Double.MAX_VALUE;
		Atom[] ret = new Atom[2];
		
		Iterator<Atom> it1 = group1.iterator();
		while (it1.hasNext()) {
			Atom atom1 = it1.next();
			Iterator<Atom> it2 = group2.iterator();
			while (it2.hasNext()) {
				Atom atom2 = it2.next();
				double dis = Calc.getDistance(atom1, atom2);
				if (dis < nearestDistance) {
					nearestDistance = dis;
					ret[0] = atom1;
					ret[1] = atom2;
				}
			}
		}
		
		if (ret[0]==null) {
			return null;
		}
		
		return ret;
	}
}
