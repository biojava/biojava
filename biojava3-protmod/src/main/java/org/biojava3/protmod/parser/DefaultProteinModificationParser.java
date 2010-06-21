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
 * Created on Jun 12, 2010
 * Author: Jianjiong Gao 
 *
 */

package org.biojava3.protmod.parser;

import java.util.ArrayList;
import java.util.Collections;
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
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;

import org.biojava3.protmod.Component;
import org.biojava3.protmod.ComponentType;
import org.biojava3.protmod.ModificationCondition;
import org.biojava3.protmod.ModifiedCompound;
import org.biojava3.protmod.ModifiedCompoundImpl;
import org.biojava3.protmod.ProteinModification;

/**
 * Identify attachment modification in a 3-D structure.
 * 
 * @author Jianjiong Gao
 * @since 3.0
 */
public class DefaultProteinModificationParser
implements ProteinModificationParser {
	
	private double bondLengthTolerance = 0.4;
	
	/**
	 * 
	 * @param bondLengthTolerance tolerance of error (in Angstroms) of the
	 *  covalent bond length, when calculating the atom distance threshold.
	 */
	public void setbondLengthTolerance(final double bondLengthTolerance) {
		if (bondLengthTolerance<0) {
			throw new IllegalArgumentException("bondLengthTolerance " +
					"must be positive.");
		}
		this.bondLengthTolerance = bondLengthTolerance;
	}
	
	/**
	 * {@inheritDoc}
	 * @throws IllegalArgumentException if null structure or
	 *  potentialModification, or , or if the nodelnr is less then 0 or
	 *  larger than or equal to the number of models in the structure.
	 */
	@Override
	public List<ModifiedCompound> parse(final Structure structure, 
			final Set<ProteinModification> potentialModifications,
			final int modelnr) {
		if (structure==null) {
			throw new IllegalArgumentException("Null structure.");
		}
		
		if (potentialModifications==null) {
			throw new IllegalArgumentException("Null potentialModifications.");
		}
		
		if (modelnr >= structure.nrModels() || modelnr < 0) {
			throw new IllegalArgumentException("modelnr should be between 0 to "
					+ structure.nrModels() + " for the structure "
					+ structure.getName());
		}
		
		List<ModifiedCompound> ret = new ArrayList<ModifiedCompound>();
		
		if (potentialModifications.isEmpty()) {
			return ret;
		}
		
		List<Chain> chains = structure.getChains(modelnr);
		for (Chain chain : chains) {
			Map<Component, List<Group>> mapCompGroups = 
				getModificationGroups(chain, potentialModifications);
			
			for (ProteinModification mod : potentialModifications) {
				ModificationCondition condition = mod.getCondition();
				List<Component> components = condition.getComponents();
				if (!mapCompGroups.keySet().containsAll(components)) {
					// not all components exist for this mod.
					continue;
				}
				
				int sizeComps = components.size();
				if (sizeComps==1) {
					// modified residue
					// TODO: is this the correct logic for CROSS_LINK_1?
					List<Group> residues = mapCompGroups.get(components.get(0));
					if (residues != null) {
						for (Group residue : residues) {
							ModifiedCompound modRes = new ModifiedCompoundImpl(mod, residue);
							ret.add(modRes);
						}
					}
				} else {
					// for multiple components
					
					// find linkages first
					List<List<Atom[]>> matchedAtomsOfLinkages =
							getMatchedAtomsOfLinkages(condition, mapCompGroups);
					
					if (matchedAtomsOfLinkages.size() != condition.linkageCount()) {
						continue;
					}
					
					assembleLinkages(matchedAtomsOfLinkages, mod, ret);
				}
			}
			
			// TODO: identify additional attached groups that are not 
			// directly attached to protein residues.
		}
		
		return ret;
	}
	
	/**
	 * 
	 * @param chain {@link Chain}.
	 * @param modifications a set of {@link ProteinModification}s.
	 * @return map from component to list of corresponding residues
	 *  in the chain.
	 */
	private Map<Component, List<Group>> getModificationGroups(
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
		
		{
			// ligands
			List<Group> groups = chain.getAtomLigands();
			
			for (Group group : groups) {
				String pdbccId = group.getPDBName().trim();
				Component comp = Component.of(pdbccId, ComponentType.LIGAND);
				if (!comps.contains(comp)) {
					continue;
				}
				List<Group> gs = mapCompRes.get(comp);
				if (gs==null) {
					gs = new ArrayList<Group>();
					mapCompRes.put(comp, gs);
				}
				gs.add(group);
			}
		}
		
		{
			// residues
			List<Group> residues = chain.getSeqResGroups();
			if (residues.isEmpty()) {
				return mapCompRes;
			}
			
			// for all residues
			for (Group group : residues) {
				String pdbccId = group.getPDBName().trim();
				Component comp = Component.of(pdbccId, ComponentType.AMINOACID);
				if (!comps.contains(comp)) {
					continue;
				}
				List<Group> gs = mapCompRes.get(comp);
				if (gs==null) {
					gs = new ArrayList<Group>();
					mapCompRes.put(comp, gs);
				}
				gs.add(group);
			}

			// for N-terminal
			Group res = residues.get(0);
			Component comp = Component.of(res.getPDBName(), ComponentType.AMINOACID, true, false);
			if (comps.contains(comp)) {
				List<Group> gs = Collections.singletonList(res);
				mapCompRes.put(comp, gs);
			}
			
			// for C-terminal
			res = residues.get(residues.size()-1);
			comp = Component.of(res.getPDBName(), ComponentType.AMINOACID, false, true);
			if (comps.contains(comp)) {
				List<Group> gs = Collections.singletonList(res);
				mapCompRes.put(comp, gs);
			}
		}

		return mapCompRes;
	}
	
	/**
	 * Get matched atoms for all linkages.	
	 */
	private List<List<Atom[]>> getMatchedAtomsOfLinkages(
			ModificationCondition condition, Map<Component, List<Group>> mapCompGroups) {
		List<List<Atom[]>> matchedAtomsOfLinkages = 
				new ArrayList<List<Atom[]>>(condition.linkageCount());
		
		List<Component> components = condition.getComponents();
		List<int[]> linkages = condition.getIndicesOfLinkedComponents();
		int nLink = linkages.size();
		
		
		for (int iLink=0; iLink<nLink; iLink++) {
			int[] linkage = linkages.get(iLink);
			Component comp1 = components.get(linkage[0]);
			Component comp2 = components.get(linkage[1]);
			List<Group> groups1 = mapCompGroups.get(comp1);
			List<Group> groups2 = mapCompGroups.get(comp2);
			
			List<String[]> atomNames = condition.getLinkedAtoms(linkage[0], linkage[1]);						
			
			for (String[] pair : atomNames) {
				List<Atom[]> list = new ArrayList<Atom[]>();
				
				for (Group g1 : groups1) {
					for (Group g2 : groups2) {
						if (g1 == g2) {
							continue;
						}
						Atom[] atoms = findLinkage(g1, g2, pair[0].split(","), pair[1].split(","));
//						Atom[] atoms = findNearestAtoms(g1, g2);								
						if (atoms!=null) {
							list.add(atoms);
						}
					}
				}
				
				if (list.isEmpty()) {
					// broken linkage
					break;
				}
				
				matchedAtomsOfLinkages.add(list);
			}
		}
		
		return matchedAtomsOfLinkages;
	}
	
	/**
	 * Find a linkage between two groups within tolerance of bond length,
	 * from potential atoms.
	 * @param group1
	 * @param group2
	 * @param nameOfAtomOnGroup1
	 * @param nameOfAtomOnGroup2
	 * @return an array of two Atoms that form bond between each other
	 *  if found; null, otherwise.
	 */
	private Atom[] findLinkage(final Group group1, final Group group2,
			String[] potentialNamesOfAtomOnGroup1, String[] potentialNamesOfAtomOnGroup2) {
		Atom[] ret = null;
		double minDistance = Double.POSITIVE_INFINITY;
		
		for (String namesOfAtomOnGroup1 : potentialNamesOfAtomOnGroup1) {
			for (String namesOfAtomOnGroup2 : potentialNamesOfAtomOnGroup2) {
				Atom[] atoms = findLinkage(group1, group2, namesOfAtomOnGroup1, namesOfAtomOnGroup2);
				if (atoms != null) {
					double distance;
					try {
						distance = Calc.getDistance(atoms[0], atoms[1]);
					} catch (StructureException e) {
						continue;
					}
					
					if (distance < minDistance) {
						minDistance = distance;
						ret = atoms;
					}
				}
			}
		}
		
		return ret;
	}
	
	/**
	 * Find a linkage between two groups within tolerance of bond length.
	 * @param group1
	 * @param group2
	 * @param nameOfAtomOnGroup1
	 * @param nameOfAtomOnGroup2
	 * @return an array of two Atoms that form bond between each other
	 *  if found; null, otherwise.
	 */
	private Atom[] findLinkage(final Group group1, final Group group2,
			String nameOfAtomOnGroup1, String nameOfAtomOnGroup2) {
		Atom[] ret = new Atom[2];
		double distance;
		
		try {
			ret[0] = group1.getAtom(nameOfAtomOnGroup1);
			ret[1] = group2.getAtom(nameOfAtomOnGroup2);
			distance = Calc.getDistance(ret[0], ret[1]);
		} catch (StructureException e) {
			return null;
		}
		
		if (ret[0]==null || ret[1]==null) {
			return null;
		}
		
		float radiusOfAtom1 = ret[0].getElement().getCovalentRadius();
		float radiusOfAtom2 = ret[1].getElement().getCovalentRadius();
		
		if (Math.abs(distance-radiusOfAtom1 -radiusOfAtom2)
				> bondLengthTolerance) {
			return null;
		}
		
		return ret;
	}
	
	/**
	 * Find the N nearest pairs of Atoms between a pair of {@link Group}s.
	 * This function is used for DEBUG only.
	 * 
	 * @param group1
	 * @param group2
	 * @return a list of a pair of Atoms if found, null otherwise.
	 */
	private Atom[] findNearestAtoms(Group group1, Group group2) {		
		double nearestDistance = Double.MAX_VALUE;
		Atom[] ret = new Atom[2];
		
		Iterator<Atom> it1 = group1.iterator();
		while (it1.hasNext()) {
			Atom atom1 = it1.next();
			Iterator<Atom> it2 = group2.iterator();
			while (it2.hasNext()) {
				Atom atom2 = it2.next();
				double dis;
				try {
					dis = Calc.getDistance(atom1, atom2);
				} catch (StructureException e) {
					continue;
				}
				if (dis < nearestDistance) {
					nearestDistance = dis;
					ret[0] = atom1;
					ret[1] = atom2;
				}
			}
		}
		
		if (ret[0]==null || ret[1]==null) {
			return null;
		}
		
		float radiusOfAtom1 = ret[0].getElement().getCovalentRadius();
		float radiusOfAtom2 = ret[1].getElement().getCovalentRadius();
		
		if (Math.abs(nearestDistance-radiusOfAtom1 -radiusOfAtom2)
				> bondLengthTolerance) {
			return null;
		}
		
		return ret;
	}
	
	/**
	 * Assembly the matched linkages.
	 * @param matchedAtomsOfLinkages
	 * @param mod
	 * @param condition
	 * @param ret ModifiedCompound will be stored here.
	 */
	private void assembleLinkages(List<List<Atom[]>> matchedAtomsOfLinkages,
			ProteinModification mod, List<ModifiedCompound> ret) {
		ModificationCondition condition = mod.getCondition();
		List<int[]> indicesOfComponentsForLinkages = 
				getIndicesOfComponentsForLinkages(condition);
		
		int nLink = matchedAtomsOfLinkages.size();
		int[] indices = new int[nLink];
		Set<ModifiedCompound> identifiedCompounds = new HashSet<ModifiedCompound>();
		while (indices[0]<matchedAtomsOfLinkages.get(0).size()) {
			List<Atom[]> atomLinkages = new ArrayList<Atom[]>(nLink);
			List<Group> groups = new ArrayList<Group>();
			for (int iLink=0; iLink<nLink; iLink++) {
				Atom[] atoms = matchedAtomsOfLinkages.get(iLink).get(indices[iLink]);
				atomLinkages.add(atoms);
				groups.add(atoms[0].getParent());
				groups.add(atoms[1].getParent());
			}
			if (matchLinkages(indicesOfComponentsForLinkages, atomLinkages)) {
				// matched
				ModifiedCompound mc = new ModifiedCompoundImpl(mod, 
						new ArrayList<Group>(groups), atomLinkages);
				if (!identifiedCompounds.contains(mc)) {
					ret.add(mc);
					identifiedCompounds.add(mc);
				}
			}
			
			// indices++ (e.g. [0,0,1]=>[0,0,2]=>[1,2,0])
			int i = nLink-1;
			while (i>=0) {
				if (i==0 || indices[i]<matchedAtomsOfLinkages.get(i).size()-1) {
					indices[i]++;
					break;
				} else {
					indices[i] = 0;
					i--;
				}
			}
		}
	}
	
	private List<int[]> getIndicesOfComponentsForLinkages(ModificationCondition condition) {
		List<int[]> ret = new ArrayList<int[]>();
		List<int[]> indicesOfLinkedComps = condition.getIndicesOfLinkedComponents();
		for (int[] indices : indicesOfLinkedComps) {
			int n = condition.getLinkedAtoms(indices[0], indices[1]).size();
			for (int i=0; i<n; i++) {
				ret.add(indices);
			}
		}
		return ret;
	}
	
	/**
	 * 
	 * @param condition
	 * @param atomLinkages
	 * @return true if atomLinkages satisfy the condition; false, otherwise.
	 */
	private boolean matchLinkages(List<int[]> indicesOfComponentsForLinkages, 
			List<Atom[]> atomLinkages) {
		int nLink = indicesOfComponentsForLinkages.size();
		if (nLink != atomLinkages.size()) {
			return false;
		}
		for (int i=0; i<nLink-1; i++) {
			int[] ix1 = indicesOfComponentsForLinkages.get(i);
			Atom[] atoms1 = atomLinkages.get(i);
			for (int j=i+1; j<nLink; j++) {
				int[] ix2 = indicesOfComponentsForLinkages.get(j);
				Atom[] atoms2 = atomLinkages.get(j);
				if (((ix1[0]==ix2[0]) != (atoms1[0].getParent()==atoms2[0].getParent()))
						 ||((ix1[0]==ix2[1]) != (atoms1[0].getParent()==atoms2[1].getParent()))
						 ||((ix1[1]==ix2[0]) != (atoms1[1].getParent()==atoms2[0].getParent()))
						 ||((ix1[1]==ix2[1]) != (atoms1[1].getParent()==atoms2[1].getParent()))) {
					return false;
				}
			}
		}
		
		return true;
	}
}
