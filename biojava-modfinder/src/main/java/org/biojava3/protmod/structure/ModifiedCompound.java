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

import java.util.Collection;
import java.util.Set;

import org.biojava3.protmod.ProteinModification;

/**
 * Root interface for all modifications in structure.
 * 
 * @author Jianjiong Gao
 * @since 3.0
 */
public interface ModifiedCompound {
	
	/** return a description of this compound
	 * 
	 * @return a description
	 */
	public String getDescription();
	
	public void setDescription(String desc);
	
	/**
	 * 
	 * @return {@link ProteinModificationBean} occurred on the residue.
	 */
	public ProteinModification getModification();
	
	public void setModification(ProteinModification modi);
	
	/**
	 * 
	 * @return a set of involved group.
	 */
	public Set<StructureGroup> getGroups();
	
	public void setGroups(Set<StructureGroup> groups);
	
	/**
	 * 
	 * @param isAminoAcid true if amino acids.
	 * @return a set of involved group of the type.
	 */
	public Set<StructureGroup> getGroups(boolean isAminoAcid);
	
	
	
	/**
	 * 
	 * @return a set of atom linkages.
	 * @see #getLinkedGroupPairs
	 * @see StructureAtomLinkage
	 */
	public Set<StructureAtomLinkage> getAtomLinkages();
	
	/** Set atom linkages
	 * 
	 * @return  
	 */
	public void setAtomLinkages(Set<StructureAtomLinkage> linkages);
	
	/**
	 * Add a linkage. Add new the involved groups first using {@link addGroup}. 
	 * @param linkage an atom linkage.
	 * @return true if this linkage was not already contained.
	 * @see StructureAtomLinkage
	 */
	public boolean addAtomLinkage(StructureAtomLinkage linkage);
	
	/**
	 * Add a collections of linkages.
	 * @param linkages an atom linkage.
	 */
	public void addAtomLinkages(Collection<StructureAtomLinkage> linkages);
	
	
	/**
	 * 
	 * @return true if groups from multiple chains were involved
	 */
	public boolean crossChains();
	
	
}
