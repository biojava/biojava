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

/**
 * Parent interface for all modifications in structure.
 */
public interface ModifiedCompound {
	/**
	 * 
	 * @return {@link ProteinModification} occurred on the residue.
	 */
	public ProteinModification getModificationType();
	
	/**
	 * 
	 * @return the chemical {@link Group}s that are involved.
	 */
	public Group[] getGroups();
	
	/**
	 * 
	 * @return an N by 2 2-dimensional array, which represent the 
	 *  atom bonds that links the {@link AminoAcid}s and/or 
	 *  the attached {@link Group}s. <i>N</i> is the number of bonds;
	 *  <i>2</i> represents the pair of atom that form a bond.
	 */
	public Atom[][] getAtomBonds();
}
