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
 * Created on Jun 7, 2010
 * Author: Jianjiong Gao 
 *
 */

package org.biojava3.protmod;

import java.util.ArrayList;
import java.util.List;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Group;

/**
 * Factory class to create and access {@link ModifiedCompound}s.
 * 
 * @author Jianjiong Gao
 * @since 3.0
 */
public final class ModifiedCompoundFactory {
	private ModifiedCompoundFactory() {
		throw new AssertionError();
	}

	/**
	 * Create a ModifiedCompound representing a modified residue.
	 * @param modification {@link ProteinModification}.
	 * @param modifiedResidue modified {@link Group}.
	 * @return a {@link ModifiedCompound}.
	 * @throws IllegalArgumentException if the modification is not a 
	 *  CHEMICAL_MODIFICATION.
	 */
	public static ModifiedCompound createModifiedResidue(
			final ProteinModification modification,
			final Group modifiedResidue) {
		if (modification.getCategory()
				!=ModificationCategory.CHEMICAL_MODIFICATION) {
			throw new IllegalArgumentException("This is not a CHEMICAL_MODIFICATION.");
		}
		List<Group> residues = new ArrayList<Group>(1);
		residues.add(modifiedResidue);
		return new ModifiedCompoundImpl(modification, residues, null, null);
	}
	
	/**
	 * Create a ModifiedCompound representing a attachemnt modification.
	 * @param modification {@link ProteinModification}.
	 * @param modifiedResidue modified {@link Group}.
	 * @param atomOnResidue {@link Atom} on the modified residue.
	 * @param attachedGroup attached chemical {@link Group}.
	 * @param atomOnAttachedGroup {@link Atom} on the attached group.
	 * @return a {@link ModifiedCompound}.
	 * @throws IllegalArgumentException if the modification is not a 
	 *  ATTACHMENT, or any argument is null.
	 */
	public static ModifiedCompound createAttachmentModification(
			final ProteinModification modification,
			final Group modifiedResidue, 
			final Atom atomOnResidue,
			final Group attachedGroup, 
			final Atom atomOnAttachedGroup) {
		if (modification==null || modifiedResidue==null || atomOnResidue==null ||
				attachedGroup==null || atomOnAttachedGroup==null) {
			throw new IllegalArgumentException("Null argument(s).");
		}
		
		if (modification.getCategory()!=ModificationCategory.ATTACHMENT) {
			throw new IllegalArgumentException("This is not a ATTACHMENT.");
		}
		
		List<Group> residues = new ArrayList<Group>(1);
		residues.add(modifiedResidue);
		
		List<Group> groups = new ArrayList<Group>(1);
		groups.add(attachedGroup);
		
		List<Atom[]> atomBonds = new ArrayList<Atom[]>(1);
		atomBonds.add(new Atom[]{atomOnResidue, atomOnAttachedGroup});
		
		return new ModifiedCompoundImpl(modification, residues, groups, atomBonds);
	}
	
	/**
	 * Create a ModifiedCompound representing a cross link.
	 * @param modification {@link ProteinModification}.
	 * @param linkedResidues protein residues.
	 * @param attachedGroups attached chemical {@link Group}s, except residues. 
	 * @param atomBonds a list of atom pairs, which represent the 
	 *  atom bonds that links the residues and/or the attached groups.
	 *  Each element of the list is a array containing two atoms. 
	 * @return a {@link ModifiedCompound} representing a cross link.
	 * @throws IllegalArgumentException if null modification, or
	 *  if null or empty linkedresidues, or if null or empty atomBonds,
	 *  or if it is not a cross link.
	 */
	public static ModifiedCompound createCrossLink(
			final ProteinModification modification,
			final List<Group> linkedResidues, 
			final List<Group> attachedGroups, 
			final List<Atom[]> atomBonds) {
		if (modification==null) {
			throw new IllegalArgumentException("Null modification.");
		}
		
		if (linkedResidues==null || linkedResidues.isEmpty()) {
			throw new IllegalArgumentException("Null or empty linked residues.");
		}
		
		if (atomBonds==null || atomBonds.isEmpty()) {
			throw new IllegalArgumentException("Null or empty atom bonds.");
		}
		
		if (!modification.getCategory().isCrossLink()) {
			throw new IllegalArgumentException("This is not a cross link.");
		}
		
		return new ModifiedCompoundImpl(modification, linkedResidues, attachedGroups, atomBonds);
	}
}
