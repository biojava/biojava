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

public interface Component {
	/**
	 * 
	 * @return Protein Data Bank ID.
	 */
	public String getPdbccId();
	
	/**
	 * 
	 * @return true if occurring on N terminal; false, otherwise.
	 */
	public boolean isNTerminal();

	/**
	 * 
	 * @return true if occurring on C terminal; false, other wise.
	 */
	public boolean isCTerminal();
	
	/**
	 * 
	 * @return the component type.
	 */
	public ComponentType getType();
}
