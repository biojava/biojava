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

public class ModifiedResidueImpl extends ModifiedCompoundImpl 
		implements ModifiedResidue	{
	/**
	 * 
	 * @param modification {@link ProteinModification}.
	 * @param modifiedAminoAcid modified {@link AminoAcid}.
	 * @throws IllegalArgumentException if the modification is not a 
	 *  CHEMICAL_MODIFICATION.
	 */
	public ModifiedResidueImpl(final ProteinModification modification,
			final AminoAcid modifiedAminoAcid) {
		super(checkType(modification), new AminoAcid[]{modifiedAminoAcid}, null);
	}
	
	/**
	 * 
	 * @param modification {@link ProteinModification}.
	 * @return the same {@link ProteinModification}.
	 * @throws IllegalArgumentException if the modification is not a 
	 *  CHEMICAL_MODIFICATION.
	 */
	private static ProteinModification checkType(ProteinModification modification) {
		if (modification.getCategory()
				!=ModificationCategory.CHEMICAL_MODIFICATION) {
			throw new IllegalArgumentException("This is not a CHEMICAL_MODIFICATION.");
		}
		return modification;
	}
	
	/**
	 * {@inheritDoc}
	 */
	public AminoAcid getModifiedAminoAcid() {
		return (AminoAcid)getGroups()[0];
	}
}
