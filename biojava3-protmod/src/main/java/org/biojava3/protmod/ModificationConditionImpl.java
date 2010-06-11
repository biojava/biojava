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

import java.util.Arrays;
import java.util.Collections;
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
	private final List<ModificationLinkage> linkages;
	
	public ModificationConditionImpl(final Component component1,
			final Component component2, final String atom1, final String atom2) {
		this(new Component[]{component1, component2},
			 new ModificationLinkage[] {
				new ModificationLinkage(component1, component2, atom1, atom2)});
	}
	
	public ModificationConditionImpl(final Component[] components,
			final ModificationLinkage[] linkages) {
		this(components==null?null:Arrays.asList(components),
				linkages==null?null:Arrays.asList(linkages));
	}
	
	/**
	 * 
	 * @param components components involved.
	 * @param bonds atom bonds.
	 * @throws IllegalArgumentException if components is null or empty,
	 *  or bonds have component(s) that are not included. 
	 */
	public ModificationConditionImpl(final List<Component> components,
			final List<ModificationLinkage> linkages) {
		
		checkComponentsAndBondsProper(components, linkages);
		
		this.components = components;
		this.linkages = linkages;
	}
	
	/**
	 * 
	 * @param components involved components.
	 * @param bonds atom bonds.
	 */
	private void checkComponentsAndBondsProper(final List<Component> components,
			final List<ModificationLinkage> bonds) {
		if (components==null||components.isEmpty()) {
			throw new IllegalArgumentException("Null or empty components.");
		}
		
		if (bonds==null||bonds.isEmpty()) {
			return;
		}
		
		Set<Component> set = new HashSet<Component>(components);
		
		for (ModificationLinkage bond:bonds) {
			if (!set.contains(bond.getComponent1())
					||!set.contains(bond.getComponent2())) {
				throw new IllegalArgumentException("Atoms must be on the " +
					"involved components.");
			}
		}
	}
	
	/**
	 * 
	 * {@inheritDoc}}
	 */
	@Override
	public List<Component> getComponents() {
		return Collections.unmodifiableList(components);
	}
	
	/**
	 * {@inheritDoc}}
	 */
	@Override
	public List<ModificationLinkage> getLinkages() {
		return linkages==null ? null : Collections.unmodifiableList(linkages);
	}
}
