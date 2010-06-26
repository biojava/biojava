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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Group;

/**
 * 
 * @author Jianjiong Gao
 * @since 3.0
 */
public class ModifiedCompoundImpl implements ModifiedCompound {
	private final ProteinModification modification;
	private final Set<Group> groups;
	private final List<Atom[]> atomLinkages;
	
	/**
	 * Create a ModifiedCompoundImpl that has only one involved component.
	 * Use this constructor for a modified residue.
	 * @param modification {@link ProteinModification}.
	 * @param modifiedResidue modified {@link Group}.
	 * @return a {@link ModifiedCompound}.
	 * @throws IllegalArgumentException if the modification is not a 
	 *  CHEMICAL_MODIFICATION.
	 */
	public ModifiedCompoundImpl (
			final ProteinModification modification,
			final Group modifiedResidue) {
		this(modification, new Group[] {modifiedResidue}, null);
	}
	
	/**
	 * Create a ModifiedCompoundImpl representing a attachemnt modification.
	 * @param modification {@link ProteinModification}.
	 * @param modifiedResidue modified {@link Group}.
	 * @param atomOnResidue {@link Atom} on the modified residue.
	 * @param attachedGroup attached chemical {@link Group}.
	 * @param atomOnAttachedGroup {@link Atom} on the attached group.
	 * @return a {@link ModifiedCompound}.
	 * @throws IllegalArgumentException if the modification is not a 
	 *  ATTACHMENT, or any argument is null.
	 */
	public ModifiedCompoundImpl (
			final ProteinModification modification,
			final Group modifiedResidue, 
			final Atom atomOnResidue,
			final Group attachedGroup, 
			final Atom atomOnAttachedGroup) {
		this(modification, new Group[] {modifiedResidue, attachedGroup},
				new Atom[][]{new Atom[]{atomOnResidue, atomOnAttachedGroup}});
	}

	/**
	 * 
	 * @param modification {@link ProteinModification}.
	 * @param otherGroups involved chemical {@link Group}s. 
	 * @param atomLinkages an Nx2 array of atom pairs, which represent the 
	 *  atom bonds that links the residues and/or the attached groups. 
	 */
	public ModifiedCompoundImpl(final ProteinModification modification,
			final Group[] groups,
			final Atom[][] atomLinkages) {		
		this(modification, new LinkedHashSet<Group>(Arrays.asList(groups)),
				atomLinkages==null?null:Arrays.asList(atomLinkages));
	}	
	
	/**
	 * 
	 * @param modification {@link ProteinModification}.
	 * @param groups involved chemical {@link Group}s residues. 
	 * @param atomLinkages a list of atom pairs, which represent the 
	 *  atom bonds that links the residues and/or the attached groups.
	 *  Each element of the list is a array containing two atoms. 
	 */
	public ModifiedCompoundImpl(final ProteinModification modification,
			final Set<Group> groups, final List<Atom[]> atomLinkages) {
		if (modification==null) {
			throw new IllegalArgumentException("modification cannot be null");
		}
		
		this.modification = modification;
		this.groups = new LinkedHashSet<Group>(groups);
		if (atomLinkages==null) {
			this.atomLinkages = Collections.emptyList();
		} else {
			this.atomLinkages = new ArrayList<Atom[]>(atomLinkages);
		}
		
		checkGroupAndAtomLinkagessProper();
	}
	
	/**
	 * 
	 * @param groups {@link Group}s.
	 * @param atomLinkages pairs of {@link Atom}s.
	 */
	private void checkGroupAndAtomLinkagessProper() {
		if (groups==null||groups.isEmpty()) {
			throw new IllegalArgumentException("Error in "+
					modification.getId()+": at least one involved residue.");
		}
		
		if (groups.size()>1) {
			if (atomLinkages.isEmpty()) {
				throw new IllegalArgumentException("Atom bonds are supposed to be " +
						"specified for more than one component.");
			}
			
			Set<Group> bondedGroups = new HashSet<Group>();
			for (Atom[] atoms : atomLinkages) {
				for (int j=0; j<2; j++) {
					if (atoms[j]==null) {
						throw new IllegalArgumentException("Null bond.");
					}
					
					Group g = atoms[j].getParent();
					if (!groups.contains(g)) {
						throw new IllegalArgumentException("Atoms must be on the " +
								"involved amino acids or other chemical groups.");
					}
					
					bondedGroups.add(g);
				}
			}
			
			if (!bondedGroups.containsAll(groups)) {
				throw new IllegalArgumentException("Some of the components were" +
						"not connected to others.");
			}
		}
	}
	
	/**
	 * {@inheritDoc}
	 */
	@Override
	public ProteinModification getModification() {
		return modification;
	}
	
	/**
	 * {@inheritDoc}
	 */
	@Override
	public Set<Group> getGroups() {
		return Collections.unmodifiableSet(groups);
	}
	
	/**
	 * {@inheritDoc}
	 */
	@Override
	public List<Atom[]> getAtomLinkages() {
		return Collections.unmodifiableList(atomLinkages);
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	public boolean addAtomLinkage(final Atom atom1, final Atom atom2) {
		if (atom1==null || atom2==null) {
			throw new IllegalArgumentException("Null atom(s).");
		}
		
		atomLinkages.add(new Atom[]{atom1, atom2});
		
		groups.add(atom1.getParent());
		groups.add(atom2.getParent());
		
		return true;
	}
	
	/**
	 * 
	 */
	@Override
	public String toString() {
		String str = "Modification: " + modification.toString() +
			   "\nGroups: " + groups;
		if (!atomLinkages.isEmpty()) {
			str += "\nAtoms: " + atomLinkages;
		}

		return str;
	}
	
	/**
	 * @return true if same modification and same components; false, otherwise.
	 */
	@Override
	public boolean equals(Object obj) {
		if (!(obj instanceof ModifiedCompoundImpl)) {
			return false;
		}
		
		ModifiedCompoundImpl mci = (ModifiedCompoundImpl)obj;
		if (mci.getModification() != modification) {
			return false;
		}
		
		Set<Group> gs = mci.getGroups();
		return gs.containsAll(groups);
		
		// Do not need to consider linkage, since they can be determined by
		// modification and groups.
	}
	
	/**
	 * 
	 */
	public int hashCode() {
		int hashCode = modification.hashCode();
		for (Group group : groups) {
			hashCode += group.hashCode();
		}
		return hashCode;
	}
}
