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

import org.biojava.bio.structure.AminoAcid;
import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Group;

public class CrossLinkImpl extends ModifiedCompoundImpl
		implements CrossLink {
	private final int nAminoAcids;
	
	/**
	 * 
	 * @param modification {@link ProteinModification}.
	 * @param linkedAminoAcids involved {@link AminoAcid}s.
	 * @param attachedGroups involved chemical {@link Group}s. 
	 * @param atomBonds an N by 2 2-dimensional array, which represent the 
	 *  atom bonds that links the {@link AminoAcid}s and/or 
	 *  the attached {@link Group}s. <i>N</i> is the number of bonds;
	 *  <i>2</i> represents the pair of atom that form a bond.
	 *  @throws IllegalArgumentException if modification is null.
	 */
	public CrossLinkImpl(final ProteinModification modification,
			final AminoAcid[] linkedAminoAcids, final Group[] attachedGroups, 
			final Atom[][] atomBonds) {
		super(checkType(modification), concatGroups(linkedAminoAcids,
				attachedGroups), atomBonds);
		nAminoAcids = linkedAminoAcids.length;
	}
	
	/**
	 * 
	 * @param modification {@link ProteinModification}.
	 * @return the same {@link ProteinModification}.
	 * @throws IllegalArgumentException if the modification is not a 
	 *  ATTACHMENT.
	 */
	private static ProteinModification checkType(ProteinModification modification) {
		if (!modification.getCategory().isCrossLink()) {
			throw new IllegalArgumentException("This is not a CrossLink.");
		}
		return modification;
	}
	
	/**
	 * 
	 * @param aminoAcids {@link AminoAcid}s.
	 * @param groups chemical {@link Group}s.
	 * @return the concatenated Groups containing both aminoAcids and groups.
	 */
	private static Group[] concatGroups(final AminoAcid[] aminoAcids,
			final Group[] groups) {
		if (aminoAcids==null||aminoAcids.length==0) {
			throw new IllegalArgumentException("Null or empty aminoAcids.");
		}
		
		if (groups==null||groups.length==0) {
			return aminoAcids;
		}
		
		Group[] res = new Group[aminoAcids.length+groups.length];
		System.arraycopy(aminoAcids, 0, res, 0, aminoAcids.length);
		System.arraycopy(groups, 0, res, aminoAcids.length, groups.length);
		return res;
	}
	
	/**
	 * 
	 * @return the {@link AminoAcid}s that are involved.
	 */
	public AminoAcid[] getLinkedAminoAcids() {
		AminoAcid[] ret = new AminoAcid[nAminoAcids];
		System.arraycopy(getGroups(), 0, ret, 0, nAminoAcids);
		return ret;
	}
	
	/**
	 * 
	 * @return the attached {@link Group}s that are involved,
	 *  excluding {@link AminoAcid}s.
	 */
	public Group[] getAttachedGroups() {
		Group[] groups = getGroups();
		int nGroups = groups.length - nAminoAcids;
		if (nGroups==0) {
			// no attached group.
			return null;
		}
		
		Group[] ret = new Group[nGroups];
		System.arraycopy(groups, nAminoAcids, ret, 0, nGroups);
		return ret;
	}
}
