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

import java.util.Set;
import java.util.HashSet;

import org.biojava.bio.structure.AminoAcid;
import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Group;

public class ModifiedCompoundImpl {
	private final ProteinModification modification;
	private final Group[] groups;
	private final Atom[][] atomBonds;
	
	/**
	 * 
	 * @param modification {@link ProteinModification}.
	 * @param groups involved chemical {@link Group}s. 
	 * @param atomBonds an N by 2 2-dimensional array, which represent the 
	 *  atom bonds that links the {@link AminoAcid}s and/or 
	 *  the attached {@link Group}s. <i>N</i> is the number of bonds;
	 *  <i>2</i> represents the pair of atom that form a bond.
	 *  @throws IllegalArgumentException if modification is null.
	 */
	public ModifiedCompoundImpl(final ProteinModification modification,
			final Group[] groups, final Atom[][] atomBonds) {
		if (modification==null) {
			throw new IllegalArgumentException("modification cannot be null");
		}
		
		checkGroupAndAtomBondsProper(groups, atomBonds);
		
		this.modification = modification;
		this.groups = groups;
		this.atomBonds = atomBonds;
	}
	
	/**
	 * 
	 * @param groups {@link Group}s.
	 * @param atomBonds pairs of {@link Atom}s.
	 */
	private void checkGroupAndAtomBondsProper(final Group[] groups,
			final Atom[][] atomBonds) {
		if (groups==null||groups.length==0) {
			throw new IllegalArgumentException("At least one involved chemical group.");
		}
		
		if (atomBonds==null) {
			return;
		}
		
		int n = atomBonds.length;
		if (n>0&&atomBonds[0].length!=2) {
			throw new IllegalArgumentException("atomBonds must be N by 2.");
		}
		
		Set<Group> gs = new HashSet<Group>();
		for (Group g:groups) {
			if (g==null) {
				throw new IllegalArgumentException("Null group.");
			}
			gs.add(g);
		}
		
		for (int i=0; i<n; i++) {
			for (int j=0; j<2; j++) {
				if (atomBonds[i][j]==null) {
					throw new IllegalArgumentException("Null bond.");
				}
				if (!gs.contains(atomBonds[i][j].getParent())) {
					throw new IllegalArgumentException("Atoms must be on the " +
							"involved amino acids or other chemical groups.");
				}
			}
		}
	}
	
	/**
	 * {@inheritDoc}
	 */
	public ProteinModification getModification() {
		return modification;
	}
	
	/**
	 * {@inheritDoc}
	 */
	public Group[] getGroups() {
		return groups;
	}
	
	/**
	 * {@inheritDoc}
	 */
	public Atom[][] getAtomBonds() {
		return atomBonds;
	}
}
