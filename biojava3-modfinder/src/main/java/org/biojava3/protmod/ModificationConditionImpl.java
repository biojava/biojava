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
	
	public ModificationConditionImpl(final List<Component> components,
			final List<ModificationLinkage> linkages) {
		
		if ( components == null)
			throw new IllegalArgumentException("Can not create ModificationCondition, components == null!");
		
		if ( components.isEmpty())
			throw new IllegalArgumentException("Can not create ModificationCondition, components is empty!");
		
		
		if (components.size() > 1) {
			Set<Integer> indices = new HashSet<Integer>();
			for (ModificationLinkage linkage : linkages) {
				indices.add(linkage.getIndexOfComponent1());
				indices.add(linkage.getIndexOfComponent2());
			}
			
			// TODO: a more comprehensive check would be checking whether 
			// all components are connected
			if (indices.size()!=components.size()) {
				throw new IllegalStateException("All components " +
						"have to be linked. indices.size:" + indices.size() + " components size:" + components.size()); // TODO: is this true?
			}
		}
		
		this.components = Collections.unmodifiableList(components);
		if (linkages==null) {
			this.linkages = Collections.emptyList();
		} else {
			this.linkages = Collections.unmodifiableList(linkages);
		}
	}
	
	/**
	 * 
	 * {@inheritDoc}}
	 */
	@Override
	public List<Component> getComponents() {
		return components;
	}
	
	/**
	 * 
	 * {@inheritDoc}}
	 */
	@Override
	public List<ModificationLinkage> getLinkages() {
		return linkages;
	}
	
	/**
	 * 
	 */
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append("Components:");
		for (Component comp : components) {
            sb.append(comp).append(";");
		}
		sb.deleteCharAt(sb.length()-1);
		
		if (!linkages.isEmpty()) {
			sb.append("\nLinkages:");
			for (ModificationLinkage link : linkages) {
                sb.append(link).append(";");
			}
			sb.deleteCharAt(sb.length()-1);
		}
		
		return sb.toString();
	}
}
