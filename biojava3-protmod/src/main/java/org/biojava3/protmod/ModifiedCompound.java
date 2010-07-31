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

import java.util.List;
import java.util.Set;

import org.biojava.bio.structure.PDBResidueNumber;

/**
 * Root interface for all modifications in structure.
 * 
 * @author Jianjiong Gao
 * @since 3.0
 */
public interface ModifiedCompound {
	/**
	 * 
	 * @return {@link ProteinModification} occurred on the residue.
	 */
	public ProteinModification getModification();
	
	/**
	 * 
	 * @return a set of the involved residues.
	 * @see PDBResidueNumber
	 */
	public Set<PDBResidueNumber> getResidues();
	
	/**
	 * 
	 * @return a set of the involved ligands.
	 * @see PDBResidueNumber
	 */
	public Set<PDBResidueNumber> getLigands();
	
	/**
	 * 
	 * @param group a group.
	 * @return true if the group is contained.
	 */
	public boolean containsGroup(PDBResidueNumber group);
	
	/**
	 * 
	 * @return a list of pairs of linked atoms.
	 * @see #getLinkedGroupPairs
	 * @see PDBAtom
	 */
	public List<PDBAtom[]> getAtomLinkages();
	
	/**
	 * 
	 * @param group1 the first group.
	 * @param group2 the second group.
	 * @return true if group1 and group2 are linked.
	 * @see PDBResidueNumber
	 */
	public boolean areLinkedGroups(PDBResidueNumber group1, PDBResidueNumber group2);
	
	/**
	 * 
	 * @param atom1 the first atom.
	 * @param atom2 the second atom.
	 * @return trye if atom1 and atom2 are linked.
	 * @see PDBAtom
	 */
	public boolean areLinkedAtoms(PDBAtom atom1, PDBAtom atom2);
	
	/**
	 * 
	 * @param group a involved group.
	 * @param isResidue true if the group is a residue; false if it is a ligand.
	 * @return true if this group was not already contained.
	 */
	public boolean addGroup(PDBResidueNumber group, boolean isResidue);
	
	/**
	 * Add a linkage. Add new the involved groups first using {@link addGroup}. 
	 * @param atom1 the first atom.
	 * @param atom2 the second atom.
	 * @return true if this linkage was not already contained.
	 * @see PDBAtom
	 */
	public boolean addAtomLinkage(PDBAtom atom1, PDBAtom atom2);
}
