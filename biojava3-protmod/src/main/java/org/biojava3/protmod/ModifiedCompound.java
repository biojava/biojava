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

import java.util.List;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Group;

/**
 * Root interface for all modifications in structure.
 * 
 * @author Jianjiong Gao
 * @since 3.0
 */
public interface ModifiedCompound {
	/**
	 * 
	 * @return {@link ProteinModification} occurred on the residue.
	 */
	public ProteinModification getModification();
	
	/**
	 * 
	 * @return involved chemical {@link Group}s.
	 */
	public List<Group> getGroups();
	
	/**
	 * 
	 * @return a list of atom pairs, which represent the 
	 *  atom bonds that links the residues and/or the attached groups.
	 *  Each element of the list is a array containing two atoms.
	 */
	public List<Atom[]> getAtomLinkages();
}
