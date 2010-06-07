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
import java.util.List;
import java.util.Set;

/**
 * 
 * @author Jianjiong Gao
 * @since 3.0
 */
public class ModificationConditionImpl implements ModificationCondition {
	private final List<Component> components;
	private final List<AtomBond> bonds;
	
	/**
	 * 
	 * @param components components involved.
	 * @param bonds atom bonds.
	 * @throws IllegalArgumentException if components is null or empty,
	 *  or bonds have component(s) that are not included. 
	 */
	public ModificationConditionImpl(final List<Component> components,
			final List<AtomBond> bonds) {
		
		checkComponentsAndBondsProper(components, bonds);
		
		this.components = components;
		this.bonds = bonds;
	}
	
	/**
	 * 
	 * @param components involved components.
	 * @param bonds atom bonds.
	 */
	private void checkComponentsAndBondsProper(final List<Component> components,
			final List<AtomBond> bonds) {
		if (components==null||components.isEmpty()) {
			throw new IllegalArgumentException("Null or empty components.");
		}
		
		if (bonds==null||bonds.isEmpty()) {
			return;
		}
		
		Set<Component> set = new HashSet<Component>(components);
		
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
	@Override
	public List<Component> getComponents() {
		return components;
	}
	
	/**
	 * 
	 * @return atom bonds between components.
	 */
	@Override
	public List<AtomBond> getBonds() {
		return bonds;
	}
}
