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
 * store the information of a linkage between two atoms
 * on two components.
 * 
 * @author Jianjiong Gao
 * @since 3.0
 */
public final class ModificationLinkage {
	private final Component comp1;
	private final Component comp2;
	private final String atom1;
	private final String atom2;
	
	/**
	 * 
	 * @param comp1 the first component.
	 * @param comp2 the second component.
	 * @param atom1 atom on the first component.
	 * @param atom2 atom on the second component.
	 */
	public ModificationLinkage(final Component comp1, final Component comp2,
			final String atom1, final String atom2) {
		if (comp1==null || comp2==null || atom1==null || atom2==null) {
			throw new IllegalArgumentException("Null argument(s).");
		}
		
		this.comp1 = comp1;
		this.comp2 = comp2;
		this.atom1 = atom1;
		this.atom2 = atom2;
	}
	
	/**
	 * 
	 * @return the first component.
	 */
	public Component getComponent1() {
		return comp1;
	}
	
	/**
	 * 
	 * @return the second component.
	 */
	public Component getComponent2() {
		return comp2;
	}
	
	/**
	 * 
	 * @return the atom name on the first component.
	 */
	public String getAtom1() {
		return atom1;
	}
	
	/**
	 * 
	 * @return the atom name on the second component.
	 */
	public String getAtom2() {
		return atom2;
	}
}
