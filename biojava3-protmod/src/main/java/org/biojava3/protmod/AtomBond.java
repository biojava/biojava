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

/**
 * store the information of a bond link between two atoms
 * on two components.
 */
public interface AtomBond {
	/**
	 * 
	 * @return the first component.
	 */
	public Component getComponent1();
	
	/**
	 * 
	 * @return the second component.
	 */
	public Component getComponent2();
	
	/**
	 * 
	 * @return the atom name on the first component.
	 */
	public String getAtom1();
	
	/**
	 * 
	 * @return the atom name on the second component.
	 */
	public String getAtom2();
}
