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

import java.util.HashSet;
import java.util.Set;

/**
 * 
 */
public class ModificationConditionImpl implements ModificationCondition {
	private final Component[] components;
	private final AtomBond[] bonds;
	
	/**
	 * 
	 * @param components components involved.
	 * @param bonds atom bonds.
	 * @throws IllegalArgumentException if components is null or empty,
	 *  or bonds have component(s) that are not included. 
	 */
	public ModificationConditionImpl(final Component[] components,
			final AtomBond[] bonds) {
		if (components==null||components.length==0) {
			throw new IllegalArgumentException("Null or empty components.");
		}
		
		checkBondsProper(components, bonds);
		
		this.components = components;
		this.bonds = bonds;
	}
	
	/**
	 * 
	 * @param components involved components.
	 * @param bonds atom bonds.
	 */
	private void checkBondsProper(final Component[] components,
			final AtomBond[] bonds) {
		if (bonds==null||bonds.length==0) {
			return;
		}
		
		Set<Component> set = new HashSet<Component>(components.length);
		for (Component comp:components) {
			set.add(comp);
		}
		
		for (AtomBond bond:bonds) {
			if (!set.contains(bond.getComponent1())
					||!set.contains(bond.getComponent2())) {
				throw new IllegalArgumentException("Atoms must be on the " +
					"involved components.");
			}
		}
	}
	
	/**
	 * 
	 * @return the involved components.
	 */
	public Component[] getComponents() {
		return components;
	}
	
	/**
	 * 
	 * @return atom bonds between components.
	 */
	public AtomBond[] getBonds() {
		return bonds;
	}
}
