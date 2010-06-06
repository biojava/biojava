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

public class AttachmentImpl extends ModifiedCompoundImpl
		implements Attachment {	
	/**
	 * 
	 * @param modification {@link ProteinModification}.
	 * @param modifiedAminoAcid modified {@link AminoAcid}.
	 * @throws IllegalArgumentException if the modification is not a 
	 *  ATTACHMENT.
	 */
	public AttachmentImpl(final ProteinModification modification,
			final AminoAcid modifiedAminoAcid, final Atom atomOnAminoAcid,
			final Group attachedGroup, final Atom atomOnAttachedGroup) {
		super(checkType(modification), new Group[]{modifiedAminoAcid, attachedGroup},
				new Atom[][]{new Atom[]{atomOnAminoAcid, atomOnAttachedGroup}});
	}
	
	/**
	 * 
	 * @param modification {@link ProteinModification}.
	 * @return the same {@link ProteinModification}.
	 * @throws IllegalArgumentException if the modification is not a 
	 *  ATTACHMENT.
	 */
	private static ProteinModification checkType(ProteinModification modification) {
		if (modification.getCategory()!=ModificationCategory.ATTACHMENT) {
			throw new IllegalArgumentException("This is not a ATTACHMENT.");
		}
		return modification;
	}
	
	/**
	 * 
	 * @return the {@link AminoAcid} that is attached to.
	 */
	public AminoAcid getModifiedAminoAcid() {
		return (AminoAcid)getGroups()[0];
	}
	
	/**
	 * 
	 * @return the attached {@link Group}.
	 */
	public Group getAttachedGroup() {
		return getGroups()[1];
	}
	
	/**
	 * 
	 * @return the attached point on the {@link AminoAcid}.
	 */
	public Atom getAtomOnAminoAcid() {
		return getAtomBonds()[0][0];
	}
	
	/**
	 * 
	 * @return the attached point on the attached {@link Group}.
	 */
	public Atom getAtomOnAttachedGroup() {
		return getAtomBonds()[0][1];
	}
}
