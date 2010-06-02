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

package org.biojava3.ptm;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Group;

public interface ChemicalBond {
	/**
	 * 
	 * @return the first of components that were linked.
	 */
	public Group component1();
	
	/**
	 * 
	 * @return the {@Atom} that is on the first component
	 *  and forms a bond with the second one.
	 */
	public Atom atomOnComp1();

	/**
	 * 
	 * @return the second of components that were linked.
	 */
	public Group component2();
	
	/**
	 * 
	 * @return the {@Atom} that is on the second component
	 *  and forms a bond with the first one.
	 */
	public Atom atomOnComp2();
}
