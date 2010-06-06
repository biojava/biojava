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
 * Created on Jun 1, 2010
 * Author: Jianjiong Gao 
 *
 */

package org.biojava3.protmod;

import org.biojava.bio.structure.AminoAcid;
import org.biojava.bio.structure.Group;

/**
 * CrossLink between residues, e.g. disulfide bonds.
 */
public interface CrossLink extends ModifiedCompound {	
	/**
	 * 
	 * @return the {@link AminoAcid}s that are involved.
	 */
	public AminoAcid[] getLinkedAminoAcids();
	
	/**
	 * 
	 * @return the attached {@link Group}s that are involved.
	 */
	public Group[] getAttachedGroups();
}
