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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * 
 * @author Jianjiong Gao
 * @since 3.0
 */
public class ModificationConditionImpl implements ModificationCondition {
	private final List<Component> components;
	
	// TODO: is it possible that two components have more than one linkage?
	private final Map<IntegerPair,List<String[]>> linkages;
	private final int linkageCount;
	
	private transient List<int[]> indicesOfLinkedComponents = null;
	
	/**
	 * Builder pattern for adding linkages.
	 */
	public static class Builder {
		private final List<Component> components;
		private Map<IntegerPair,List<String[]>> linkages;
		private int linkageCount;
		
		/**
		 * 
		 * @param components a array of {@link Component}s.
		 */
		public Builder(final Component[] components) {
			if (components==null || components.length==0) {
				throw new IllegalArgumentException("Nul or empty components");
			}
			this.components = Arrays.asList(components);
			linkages = null;
			linkageCount = 0;
		}
		
		/**
		 * 
		 * @param components a list of {@link Component}s.
		 */
		public Builder(final List<Component> components) {
			if (components==null || components.isEmpty()) {
				throw new IllegalArgumentException("Nul or empty components");
			}
			this.components = components;
			linkages = null;
		}
		
		public Builder addLinkage(final int indexComponent1, 
				final int indexComponent2, final String atom1,
				final String atom2) {
			if (components.size()==1) {
				throw new IllegalStateException("No linkage is needed for one component.");
			}
			
			if (indexComponent1<0 || indexComponent1>=components.size()
					|| indexComponent2<0 || indexComponent2>=components.size()) {
				throw new IllegalArgumentException("Indices of components must" +
						" be an integer between 0 and the length of components (excluded)");
			}
			
			if (indexComponent1 == indexComponent2) {
				throw new IllegalArgumentException("No linkage is allowed for an" +
						" identical component.");
			}
			
			if (atom1==null || atom2==null) {
				throw new IllegalArgumentException("Null atom(s).");
			}
			
			if (linkages==null) {
				linkages = new HashMap<IntegerPair,List<String[]>>();
			}

			IntegerPair ip = IntegerPair.of(indexComponent1, indexComponent2);
			List<String[]> atoms = linkages.get(ip);
			if (atoms == null) {
				atoms = new ArrayList<String[]>();
				linkages.put(ip, atoms);
			}
			
			atoms.add(new String[] {atom1, atom2});
			
			linkageCount++;
			
			return this;
		}
		
		public ModificationConditionImpl build() {
			return new ModificationConditionImpl(this);
		}
	}
	
	private ModificationConditionImpl(final Builder builder) {
		components = builder.components;
		linkages = builder.linkages;
		if (components.size()>1) {
			if (linkages==null) {
				throw new IllegalStateException("No linkage was added.");
			}
			
			Set<Integer> indices = new HashSet<Integer>();
			for (IntegerPair ip : linkages.keySet()) {
				indices.add(ip.i1);
				indices.add(ip.i2);
			}
			
			// TODO: a more comprehensive check would be checking whether 
			// all components are connected
			if (indices.size()!=components.size()) {
				throw new IllegalStateException("All components " +
						"have to be linked."); // TODO: is this true?
			}
		}
		linkageCount = builder.linkageCount; 
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
	public List<int[]> getIndicesOfLinkedComponents() {
		if (linkages == null) {
			return Collections.emptyList();
		}
		
		if (indicesOfLinkedComponents == null) {
			indicesOfLinkedComponents = new ArrayList<int[]>(linkages.size());
			for (IntegerPair ip : linkages.keySet()) {
				indicesOfLinkedComponents.add(new int[]{ip.i1, ip.i2});
			}
		}
		
		return Collections.unmodifiableList(indicesOfLinkedComponents);
	}
	
	/**
	 * {@inheritDoc}
	 */
	public List<String[]> getLinkedAtoms(int indexComponent1, int indexComponent2) {
		if (indexComponent1<0 || indexComponent1>=components.size()
				|| indexComponent2<0 || indexComponent2>=components.size()) {
			throw new IllegalArgumentException("Indices of components must" +
					" be an integer between 0 and the length of components (excluded)");
		}

		IntegerPair ip = IntegerPair.of(indexComponent1, indexComponent2);
		List<String[]> atoms = linkages.get(ip);
		
		if (atoms == null) {
			return Collections.emptyList();
		}
		
		if (indexComponent1 < indexComponent2) {
			return atoms;
		} else {
			List<String[]> rev = new ArrayList<String[]>(atoms.size());
			for (String[] atomPair : atoms) {
				rev.add(new String[]{atomPair[1], atomPair[0]});
			}
			return rev;
		}
	}
	
	/**
	 * {@inheritDoc}
	 */
	public int linkageCount() {
		return linkageCount;
	}
	
	/**
	 * 
	 */
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append("Components:");
		for (Component comp : components) {	
			sb.append(comp+";");
		}
		
		List<int[]> indicesOfLinkedComps = getIndicesOfLinkedComponents();
		if (!indicesOfLinkedComps.isEmpty()) {
			sb.append("\tLinkages:");
			for (int[] link : indicesOfLinkedComps) {
				Component comp1 = components.get(link[0]);
				Component comp2 = components.get(link[1]);
				List<String[]> atoms = getLinkedAtoms(link[0], link[1]);
				for (String[] atomPair : atoms) {
					sb.append(comp1.getPdbccId()+"["+atomPair[0]+"]<=>"+
							comp2.getPdbccId()+"["+atomPair[1]+"];");
				}
			}
		}
		
		return sb.toString();
	}
	
	/**
	 * 
	 */
	private static class IntegerPair {
		private final int i1;
		private final int i2;
		
		private static Map<Integer, Map<Integer, IntegerPair>> instances = null;
		
		private IntegerPair(int i1, int i2) {
			this.i1 = i1;
			this.i2 = i2;
		}
		
		static IntegerPair of(int i1, int i2) {
			if (instances==null) {
				instances = new HashMap<Integer, Map<Integer, IntegerPair>>();
			}
			
			int small, large;
			if (i1 < i2) {
				small = i1;
				large = i2;
			} else {
				small = i2;
				large = i1;
			}
			
			Map<Integer, IntegerPair> map = instances.get(small);
			if (map == null) {
				map = new HashMap<Integer, IntegerPair>();
				instances.put(small, map);
			}
			
			IntegerPair ip = map.get(large);
			if (ip == null) {
				ip = new IntegerPair(small, large);
				map.put(large, ip);
			}
			
			return ip;
		}
		
		@Override
		public boolean equals(Object obj) {
			if (!(obj instanceof IntegerPair)) {
				return false;
			}
			
			IntegerPair other = (IntegerPair)obj;
			
			return (i1==other.i1 && i2==other.i2)
					|| (i1==other.i2 && i2==other.i1);
		}
		
		@Override
		public int hashCode() {
			return i1*i1*i2*i2;
		}
	}
}
