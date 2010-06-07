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

import java.util.HashSet;
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
	private final List<Group> residues;
	private final List<Group> otherGroups;
	private final List<Atom[]> atomBonds;
	
	/**
	 * 
	 * @param modification {@link ProteinModification}.
	 * @param residues protein residues.
	 * @param otherGroups involved chemical {@link Group}s, except residues. 
	 * @param atomBonds a list of atom pairs, which represent the 
	 *  atom bonds that links the residues and/or the attached groups.
	 *  Each element of the list is a array containing two atoms. 
	 */
	public ModifiedCompoundImpl(final ProteinModification modification,
			final List<Group> residues, final List<Group> otherGroups,
			final List<Atom[]> atomBonds) {
		if (modification==null) {
			throw new IllegalArgumentException("modification cannot be null");
		}
		
		checkGroupAndAtomBondsProper(residues, otherGroups, atomBonds);
		
		this.modification = modification;
		this.residues = residues;
		this.otherGroups = otherGroups;
		this.atomBonds = atomBonds;
	}
	
	/**
	 * 
	 * @param groups {@link Group}s.
	 * @param atomBonds pairs of {@link Atom}s.
	 */
	private void checkGroupAndAtomBondsProper(final List<Group> residues,
			final List<Group> otherGroups, final List<Atom[]> atomBonds) {
		if (residues==null||residues.isEmpty()) {
			throw new IllegalArgumentException("At least one involved residue.");
		}
		
		Set<Group> gs = new HashSet<Group>(residues);
		if (otherGroups != null) {
			gs.addAll(otherGroups);
		}
		
		if (gs.size()==1) {
			if (atomBonds!=null && !atomBonds.isEmpty()) {
				throw new IllegalArgumentException("No atomBond is necessary for" +
						" one component.");
			}
		} else {
			if (atomBonds==null) {
				throw new IllegalArgumentException("Atom bonds are supposed to be " +
						"specified for more than one component.");
			}
			
			for (Atom[] atoms : atomBonds) {
				for (int j=0; j<2; j++) {
					if (atoms[j]==null) {
						throw new IllegalArgumentException("Null bond.");
					}
					if (!gs.contains(atoms[j].getParent())) {
						throw new IllegalArgumentException("Atoms must be on the " +
								"involved amino acids or other chemical groups.");
					}
				}
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
	public List<Group> getProteinResidues() {
		return residues;
	}
	
	/**
	 * {@inheritDoc}
	 */
	@Override
	public List<Group> getOtherGroups() {
		return otherGroups;
	}
	
	/**
	 * {@inheritDoc}
	 */
	@Override
	public List<Atom[]> getAtomBonds() {
		return atomBonds;
	}
}
