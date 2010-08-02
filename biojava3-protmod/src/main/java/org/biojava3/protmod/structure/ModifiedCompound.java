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

package org.biojava3.protmod.structure;

import java.util.Set;

import org.biojava3.protmod.ComponentType;
import org.biojava3.protmod.ProteinModification;

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
	 * @return a set of involved group.
	 */
	public Set<StructureGroup> getGroups();
	
	/**
	 * 
	 * @param type the type of groups.
	 * @return a set of involved group of the type.
	 * @see ComponentType
	 */
	public Set<StructureGroup> getGroups(ComponentType type);
	
//	/**
//	 * 
//	 * @param group a group.
//	 * @return true if the group is contained.
//	 * @see StructureGroup
//	 */
//	public boolean containsGroup(StructureGroup group);
	
	/**
	 * 
	 * @return a set of atom linkages.
	 * @see #getLinkedGroupPairs
	 * @see StructureAtomLinkage
	 */
	public Set<StructureAtomLinkage> getAtomLinkages();
//	
//	/**
//	 * 
//	 * @param group1 the first group.
//	 * @param group2 the second group.
//	 * @return true if group1 and group2 are linked.
//	 * @see StructureGroup
//	 */
//	public boolean areLinkedGroups(StructureGroup group1, StructureGroup group2);
//	
//	/**
//	 * 
//	 * @param atom1 the first atom.
//	 * @param atom2 the second atom.
//	 * @return true if atom1 and atom2 are linked.
//	 * @see StructureAtom
//	 */
//	public boolean areLinkedAtoms(StructureAtom atom1, StructureAtom atom2);
//	
//	/**
//	 * 
//	 * @param group a involved group.
//	 * @param isResidue true if the group is a residue; false if it is a ligand.
//	 * @return true if this group was not already contained.
//	 * @see StructureGroup
//	 */
//	public boolean addGroup(StructureGroup group);
	
	/**
	 * Add a linkage. Add new the involved groups first using {@link addGroup}. 
	 * @param linkage an atom linkage.
	 * @return true if this linkage was not already contained.
	 * @see StructureAtomLinkage
	 */
	public boolean addAtomLinkage(StructureAtomLinkage linkage);
}
