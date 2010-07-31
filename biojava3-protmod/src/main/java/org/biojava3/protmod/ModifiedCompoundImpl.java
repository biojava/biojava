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
 * Created on Jun 5, 2010
 * Author: Jianjiong Gao 
 *
 */

package org.biojava3.protmod;

import java.io.Serializable;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.biojava.bio.structure.PDBResidueNumber;

/**
 * 
 * @author Jianjiong Gao
 * @since 3.0
 */
public class ModifiedCompoundImpl
implements ModifiedCompound, Serializable {
	private static final long serialVersionUID = -3809110608450334460L;

    private static class SerializationProxy implements Serializable {
        private static final long serialVersionUID = -5345107617243742042L;
        private final String modificationId;
    	private final Set<PDBResidueNumber> residues;
    	private final Set<PDBResidueNumber> ligands;
    	private final List<PDBAtom[]> atomLinkages;
        SerializationProxy(ModifiedCompoundImpl actualObj) {
            this.modificationId = actualObj.modification.getId();
            this.residues = actualObj.residues;
            this.ligands = actualObj.ligands;
            if (actualObj.atomLinkages == null)
            	this.atomLinkages = null;
            else
            	this.atomLinkages = actualObj.getAtomLinkages();
        }

        private Object readResolve() {
        	ProteinModification modification = ProteinModification.getById(modificationId);
            return new ModifiedCompoundImpl(modification, residues, ligands, atomLinkages);
        }
    }

    private Object writeReplace() {
        return new SerializationProxy(this);
    }

    private void readObject(java.io.ObjectInputStream stream)
    		throws java.io.InvalidObjectException {
        throw new java.io.InvalidObjectException("Proxy required");
    }
	
	private final ProteinModification modification;
	private final Set<PDBResidueNumber> residues;
	private final Set<PDBResidueNumber> ligands;
	private Map<Set<PDBResidueNumber>, Set<Set<PDBAtom>>> atomLinkages;
	
	/**
	 * Create a ModifiedCompoundImpl that has only one involved component.
	 * Use this constructor for a modified residue.
	 * @param modification {@link ProteinModification}.
	 * @param modifiedResidue modified group.
	 * @return a {@link ModifiedCompound}.
	 * @throws IllegalArgumentException if either argument is null.
	 */
	public ModifiedCompoundImpl (
			final ProteinModification modification,
			final PDBResidueNumber modifiedResidue) {
		if (modification==null || modifiedResidue==null) {
			throw new IllegalArgumentException("Null argument(s)");
		}
		this.modification = modification;
		
		residues = new HashSet<PDBResidueNumber>(1);
		residues.add(modifiedResidue);
		
		ligands = null;

		// it is possible that components be added by addLinkage later
		atomLinkages = null;
	}
	
	/**
	 * Create a ModifiedCompoundImpl for an attachment modification.
	 * @param modification {@link ProteinModification}.
	 * @param residue modified residue.
	 * @param atom1 atom on the residue.
	 * @param ligand attached ligand.
	 * @param atom2 atom on the ligand.
	 * @return a {@link ModifiedCompound}.
	 * @throws IllegalArgumentException if any argument is null.
	 */
	public ModifiedCompoundImpl (
			final ProteinModification modification,
			final PDBResidueNumber residue, 
			final String atom1,
			final PDBResidueNumber ligand, 
			final String atom2) {
		this (modification, 
				new PDBResidueNumber[] {residue},
				new PDBResidueNumber[] {ligand},
				Collections.singletonList(
					new PDBAtom[] {
						new PDBAtom(residue, atom1), 
						new PDBAtom(ligand, atom2)
						}
				)
			);
	}
	
	/**
	 * 
	 * @param modification {@link ProteinModification}.
	 * @param residues an array of residues.
	 * @param ligands an array of ligands.
	 * @param linkages a set of pairs of linked atoms.
	 * @see PDBAtom
	 */
	public ModifiedCompoundImpl(final ProteinModification modification,
			final PDBResidueNumber[] residues,
			final PDBResidueNumber[] ligands,
			final List<PDBAtom[]> linkages) {
		this (modification, 
				new LinkedHashSet<PDBResidueNumber>(Arrays.asList(residues)),
				new LinkedHashSet<PDBResidueNumber>(Arrays.asList(ligands)),
				linkages);
	}
	
	/**
	 * 
	 * @param modification {@link ProteinModification}.
	 * @param residues a set of residues.
	 * @param ligands a set of ligands.
	 * @param linkages a set of pairs of linked atoms.
	 * @see PDBAtom
	 */
	public ModifiedCompoundImpl(final ProteinModification modification,
			final Set<PDBResidueNumber> residues,
			final Set<PDBResidueNumber> ligands,
			final List<PDBAtom[]> linkages) {
		if (modification==null) {
			throw new IllegalArgumentException("modification cannot be null");
		}

		if (residues==null||residues.isEmpty()) {
			throw new IllegalArgumentException("at least one linkage.");
		}
		
		this.modification = modification;
		this.residues = residues;
		this.ligands = ligands;
		
		if (linkages==null) {
			this.atomLinkages = null;
		} else {
			for (PDBAtom[] linkage : linkages) {
				if (linkage==null || linkage.length!=2) {
					throw new IllegalArgumentException("Error in "+
							modification.getId()+"a linkage must contain two atoms.");
				}
				
				addAtomLinkage(linkage[0], linkage[1]);
			}
		}
	}

	@Override
	public ProteinModification getModification() {
		return modification;
	}
	
	@Override
	public Set<PDBResidueNumber> getResidues() {
		return Collections.unmodifiableSet(residues);
	}
	
	@Override
	public Set<PDBResidueNumber> getLigands() {
		if (ligands == null) {
			return Collections.emptySet();
		}
		
		return Collections.unmodifiableSet(ligands);
	}
	
	@Override
	public boolean containsGroup(PDBResidueNumber group) {
		return residues.contains(group) ||
				(ligands!=null && ligands.contains(group));
	}

	@Override
	public List<PDBAtom[]> getAtomLinkages() {
		if (atomLinkages==null) {
			return Collections.emptyList();
		} else {
			List<PDBAtom[]> result = new ArrayList<PDBAtom[]>();
			for (Set<Set<PDBAtom>> linkages : atomLinkages.values()) {
				for (Set<PDBAtom> linkage : linkages) {
					PDBAtom[] atoms = linkage.toArray(new PDBAtom[0]);
					result.add(atoms);
				}
			}
			
			return result;
		}
	}
	
	@Override
	public boolean addGroup(PDBResidueNumber group, boolean isResidue) {
		if (containsGroup(group)) {
			return false;
		}
		
		if (isResidue) {
			residues.add(group);
		} else {
			ligands.add(group);
		}
		
		return true;
	}

	@Override
	public boolean addAtomLinkage(PDBAtom atom1, PDBAtom atom2) {
		if (atom1==null || atom2==null) {
			throw new IllegalArgumentException("Null argument(s)");
		}
		
		if (!containsGroup(atom1.getGroup()) || !containsGroup(atom2.getGroup())) {
			throw new IllegalStateException("Add new groups first before adding a linkage.");
		}
		
		Set<PDBResidueNumber> groups = new HashSet<PDBResidueNumber>(2);
		groups.add(atom1.getGroup());
		groups.add(atom2.getGroup());
		
		if (atomLinkages==null) {
			atomLinkages = new HashMap<Set<PDBResidueNumber>, Set<Set<PDBAtom>>>();
		}

		Set<Set<PDBAtom>> linkages = atomLinkages.get(groups);
		if (linkages == null) {
			linkages = new HashSet<Set<PDBAtom>>();
			atomLinkages.put(groups, linkages);
		}
		
		Set<PDBAtom> atoms = new HashSet<PDBAtom>(2);
		atoms.add(atom1);
		atoms.add(atom2);

		return linkages.add(atoms);
	}
	
	@Override
	public boolean areLinkedGroups(PDBResidueNumber group1, PDBResidueNumber group2) {
		Set<PDBResidueNumber> groups = new HashSet<PDBResidueNumber>(2);
		groups.add(group1);
		groups.add(group2);
		return atomLinkages.containsKey(groups);
	}
	
	/**
	 * 
	 * @param atom1 the first atom.
	 * @param atom2 the second atom.
	 * @return true if atom1 and atom2 are linked.
	 * @see PDBAtom
	 */
	public boolean areLinkedAtoms(PDBAtom atom1, PDBAtom atom2) {
		Set<PDBResidueNumber> groups = new HashSet<PDBResidueNumber>(2);
		groups.add(atom1.getGroup());
		groups.add(atom2.getGroup());
		Set<Set<PDBAtom>> linkages = atomLinkages.get(groups);
		if (linkages == null) {
			return false;
		}
		
		Set<PDBAtom> atoms = new HashSet<PDBAtom>(2);
		atoms.add(atom1);
		atoms.add(atom2);
		
		return linkages.contains(atoms);
	}
	
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append("Modification: ");
		sb.append(modification.toString());
		sb.append("\nResidues: ");
		sb.append(residues);
		sb.append("\nLigands: ");
		sb.append(ligands);
		if (!atomLinkages.isEmpty()) {
			sb.append("\nAtom linkages:");
			for (Set<Set<PDBAtom>> linkages : atomLinkages.values()) {
				for (Set<PDBAtom> linkage : linkages) {
					sb.append("\n\t");
					PDBAtom[] atoms = linkage.toArray(new PDBAtom[0]);
					sb.append(atoms[0]);
					sb.append("<=>");
					sb.append(atoms[1]);
				}
			}
		}

		return sb.toString();
	}
	
	/**
	 * @return true if same modification and same components; false, otherwise.
	 */
	@Override
	public boolean equals(Object obj) {
		if (!(obj instanceof ModifiedCompound)) {
			return false;
		}
		
		ModifiedCompound mci = (ModifiedCompound)obj;
		if (mci.getModification() != modification) {
			return false;
		}
		
		if (!residues.equals(mci.getResidues())) {
			return false;
		}
		
		if (ligands==null) {
			if (mci.getLigands()!=null) {
				return false;
			}
		} else {
			if (!ligands.equals(mci.getLigands()))
			return false;
		}
		
		return true;
		
		// Do not need to consider linkage, since they can be determined by
		// modification and groups.
	}
	
	@Override
	public int hashCode() {
		int result = 17;
		result = result * 32 + modification.hashCode();
		
		int sum = 0;
		for (PDBResidueNumber group : residues) {
			sum += group.hashCode();
		}
		
		result = result * 32 + sum;
		
		sum = 0;
		if (ligands != null) {
			for (PDBResidueNumber group : ligands) {
				sum += group.hashCode();
			}
		}
		
		result = result * 32 + sum;
		
		return result;
	}
}
