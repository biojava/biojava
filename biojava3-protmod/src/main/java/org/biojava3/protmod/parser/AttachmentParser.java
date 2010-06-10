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
import java.util.HashMap;
import java.util.Iterator;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.GroupType;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;

import org.biojava3.protmod.Component;
import org.biojava3.protmod.ComponentType;
import org.biojava3.protmod.ModificationCategory;
import org.biojava3.protmod.ModificationCondition;
import org.biojava3.protmod.ModificationLinkage;
import org.biojava3.protmod.ModifiedCompound;
import org.biojava3.protmod.ModifiedCompoundFactory;
import org.biojava3.protmod.ProteinModification;

/**
 * Identify attachment modification in a 3-D structure.
 * 
 * @author Jianjiong Gao
 * @since 3.0
 */
public class AttachmentParser implements ProteinModificationParser {
	final double bondLengthTolerance;
	
	/**
	 * 
	 * @param bondLengthTolerance tolerance of error (in Angstroms) of the
	 *  covalent bond length, when calculating the atom distance threshold.
	 */
	public AttachmentParser(final double bondLengthTolerance) {
		if (bondLengthTolerance<0) {
			throw new IllegalArgumentException("bondLengthTolerance " +
					"must be positive.");
		}
		this.bondLengthTolerance = bondLengthTolerance;
	}
	
	/**
	 * Parse attachement modification in a structure.
	 * @param structure query {@link Structure}.
	 * @param potentialModifications query {@link ProteinModification}s.
	 * @param modelnr model number.
	 * @return an list of {@link ModifiedCompound}s, or null if the
	 *  nodelnr is larger than the number of models in the structure.
	 * @throws IllegalArgumentException if null structure, or null or 
	 *  empty potentialModifications, or potentialModifications contain 
	 *  modifications other than ATTACHMENT.
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
			if (mod.getCategory()!=ModificationCategory.ATTACHMENT) {
				throw new IllegalArgumentException("Only ATTACHMENT is allowed.");
			}
		}
		
		if (modelnr >= structure.nrModels())
			return null;
		
		List<ModifiedCompound> ret = new ArrayList<ModifiedCompound>();
		
		List<Chain> chains = structure.getChains(modelnr);
		for (Chain chain : chains) {
			Map<Component, List<Group>> mapCompRes = 
					modifiableResidues(chain, potentialModifications);
			
			List<Group> groups = chain.getAtomGroups(GroupType.HETATM);
			
			// for all heta
			for (Group group : groups) {
				String pdbccId = group.getPDBName();
				Component comp = Component.of(pdbccId);
				if (comp==null) {
					continue;
				}
				
				Set<ProteinModification> mods = ProteinModification.getByComponent(comp);
				if (mods==null) {
					continue;
				}
				
				mods = new HashSet<ProteinModification>(mods);
				mods.retainAll(potentialModifications);
				
				if (mods.isEmpty()) {
					continue;
				}
				
				for (ProteinModification mod : mods) {
					ModificationCondition condition = mod.getCondition();
					comp = condition.getComponents().get(0);
					if (comp.getType() != ComponentType.AMINOACID) {
						throw new IllegalStateException("No residue involved in modification");
					}
					
					List<Group> residues = mapCompRes.get(comp);
					if (residues == null || residues.isEmpty()) {
						continue;
					}
					
					Atom atomOnAttachedGroup = null;
					
					//* use atom specified by condition
					ModificationLinkage bond = condition.getBonds().get(0);
					String nameOfAtomOnResidue = bond.getAtom1();
					String nameOfAtomOnAttachedGroup = bond.getAtom2();
					
					try {
						atomOnAttachedGroup = group.getAtom(nameOfAtomOnAttachedGroup);
					} catch (StructureException e) {
						e.printStackTrace();
					}
					if (atomOnAttachedGroup==null) {
						System.err.println("Atom does not exist.");
						continue;
					}//*/
					
					double clostestDistance = Double.POSITIVE_INFINITY;
					Group clostestResidue = null;
					Atom closestAtomOnResidue = null;
					for (Group residue : residues) {
						//* use atom specified in condition
						Atom atomOnResidue = null;
						try {
							atomOnResidue = residue.getAtom(nameOfAtomOnResidue);
						} catch (StructureException e) {
							//e.printStackTrace();
							continue;
						}
						
						if (atomOnResidue==null) {
							//System.err.println("Atom does not exist.");
							continue;
						}
						
						double distance;
						try {
							distance = Calc.getDistance(atomOnAttachedGroup, atomOnResidue);
						} catch (StructureException e) {
							e.printStackTrace();
							continue;
						}
						
						if (distance < clostestDistance) {
							clostestDistance = distance;
							clostestResidue = residue;
							closestAtomOnResidue = atomOnResidue;
						}//*/
						
						/* find the atom
						Atom[] atoms;
						double distance;
						try {
							atoms = findNearestAtoms(residue, group);
							if (atoms==null) {
								continue;
							}
							distance = Calc.getDistance(atoms[0], atoms[1]);
						} catch (StructureException e) {
							e.printStackTrace();
							continue;
						}
						
						if (distance < clostestDistance) {
							clostestDistance = distance;
							clostestResidue = residue;
							closestAtomOnResidue = atoms[0];
							atomOnAttachedGroup = atoms[1];
						}//*/
					}
					
					if (Double.isInfinite(clostestDistance)) {
						continue;
					}
					
					float radiusOfAtomOnResidue = 
						closestAtomOnResidue.getElement().getCovalentRadius();					
					float radiusOfAtomOnAttachedGroup = 
						atomOnAttachedGroup.getElement().getCovalentRadius();
					if (Math.abs(clostestDistance-radiusOfAtomOnResidue
							-radiusOfAtomOnAttachedGroup) < bondLengthTolerance) {
						ModifiedCompound attachment = ModifiedCompoundFactory
								.createAttachmentModification(mod, clostestResidue, 
										closestAtomOnResidue, group, atomOnAttachedGroup);
						ret.add(attachment);
					}
				}
			}
			
			// TODO: identify additional attached groups that are not 
			// directly attached to protein residues.
		}
		
		return ret;
	}
	
	/**
	 * Find the nearest Atoms between a pair of {@link Group}s.
	 * @param group1
	 * @param group2
	 * @return a pair of Atoms.
	 * @throws StructureException ...
	 */
	private Atom[] findNearestAtoms(Group group1, Group group2)
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
	
	/**
	 * 
	 * @param chain {@link Chain}.
	 * @param modifications a set of {@link ProteinModification}s.
	 * @return map from component to list of corresponding residues
	 *  in the chain.
	 */
	private Map<Component, List<Group>> modifiableResidues(
			final Chain chain, 
			final Set<ProteinModification> modifications) {
		Map<Component, List<Group>> mapCompRes = new HashMap<Component, List<Group>>();
		
		List<Group> residues = chain.getSeqResGroups();
		
		if (residues==null || residues.isEmpty()) {
			return mapCompRes;
		}
		
		// for all residue
		for (Group res : residues) {
			String pdbccId = res.getPDBName();
			Component comp = Component.of(pdbccId);
			if (comp==null) {
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
		if (comp!=null) {
			List<Group> groups = new ArrayList<Group>(1);
			groups.add(res);
			mapCompRes.put(comp, groups);
		}
		
		// for C-terminal
		res = residues.get(residues.size()-1);
		comp = Component.of(res.getPDBName(), false, true);
		if (comp!=null) {
			List<Group> groups = new ArrayList<Group>(1);
			groups.add(res);
			mapCompRes.put(comp, groups);
		}

		return mapCompRes;
	}
}
