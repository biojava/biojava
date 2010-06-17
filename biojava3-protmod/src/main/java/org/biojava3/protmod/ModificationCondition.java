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
 * Created on Jun 4, 2010
 * Author: Jianjiong Gao 
 *
 */

package org.biojava3.protmod;

import java.util.List;

/**
 * Conditions of a protein modification, e.g. components and atoms.
 * 
 * @author Jianjiong Gao
 * @since 3.0
 */
public interface ModificationCondition {
	/**
	 * 
	 * @return the involved components.
	 */
	public List<Component> getComponents();
	
	/**
	 * 
	 * @return a list of indices of linked components, or null if not 
	 *  exist. Each element if an array of two integers, which are 
	 *  indices of a pair of linked components.
	 */
	public List<int[]> getIndicesOfLinkedComponents();
	
	/**
	 * 
	 * @param indexComponent1 index of component1, starting from 0.
	 * @param indexComponent2 index of component2, starting from 0.
	 * @return a list of linked atoms of the linked component. 
	 * 	Each element is an array of two Strings, which are the PDBCC
	 *  name of a pair of atoms that link the two component; or null
	 *  if the two components are not linked.
	 * 
	 * Note that two components could have multiple linkages.
	 */
	public List<String[]> getLinkedAtoms(int indexComponent1, int indexComponent2);
	
	/**
	 * 
	 * @return the total number of linkages.
	 */
	public int linkageCount();
}
