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
 * Created on Jun 21, 2010
 * Author: Jianjiong Gao 
 *
 */

package org.biojava3.protmod;

import java.util.Collections;
import java.util.List;

public class ModificationLinkage {
	private final List<Component> components;
	private final int indexOfComponent1;
	private final int indexOfComponent2;
	private final List<String> pdbNameOfPotentialAtomsOnComponent1;
	private final List<String> pdbNameOfPotentialAtomsOnComponent2;
	private final String labelOfAtomOnComponent1;
	private final String labelOfAtomOnComponent2;
	
	/**
	 * 
	 * @param components {@link Component}s involved in a modification.
	 * @param indexOfComponent1 index of the first component.
	 * @param indexOfComponent2 index of the second component.
	 */
	public ModificationLinkage(
			final List<Component> components,
			final int indexOfComponent1,
			final int indexOfComponent2) {
		this(components, indexOfComponent1, null, 
				null, indexOfComponent2, null, null);
	}
	
	/**
	 * 
	 * @param components {@link Component}s involved in a modification.
	 * @param indexOfComponent1 index of the first component.
	 * @param labelOfAtomOnComponent1 label of the atom on the first
	 *  component.
	 * @param indexOfComponent2 index of the second component.
	 * @param labelOfAtomOnComponent2 label of the atom on the second
	 *  component.
	 */
	public ModificationLinkage(
			final List<Component> components,
			final int indexOfComponent1,
			final String pdbNameOfAtomsOnComponent1,
			final int indexOfComponent2,
			final String pdbNameOfAtomsOnComponent2) {
		this(components, indexOfComponent1, 
				Collections.singletonList(pdbNameOfAtomsOnComponent1),
				null, indexOfComponent2,
				Collections.singletonList(pdbNameOfAtomsOnComponent2),
				null);
	}
	
	/**
	 * 
	 * @param components {@link Component}s involved in a modification.
	 * @param indexOfComponent1 index of the first component.
	 * @param labelOfAtomOnComponent1 label of the atom on the first
	 *  component.
	 * @param indexOfComponent2 index of the second component.
	 * @param labelOfAtomOnComponent2 label of the atom on the second
	 *  component.
	 */
	public ModificationLinkage(
			final List<Component> components,
			final int indexOfComponent1,
			final List<String> pdbNameOfPotentialAtomsOnComponent1,
			final int indexOfComponent2,
			final List<String> pdbNameOfPotentialAtomsOnComponent2) {
		this(components, indexOfComponent1, pdbNameOfPotentialAtomsOnComponent1,
				null, indexOfComponent2, pdbNameOfPotentialAtomsOnComponent2, null);
	}
	
	/**
	 * 
	 * @param components {@link Component}s involved in a modification.
	 * @param indexOfComponent1 index of the first component.
	 * @param pdbNameOfPotentialAtomsOnComponent1 a list of PDB names of
	 *  potential atoms on the first component.
	 * @param labelOfAtomOnComponent1 label of the atom on the first
	 *  component.
	 * @param indexOfComponent2 index of the second component.
	 * @param pdbNameOfPotentialAtomsOnComponent2 a list of PDB names of
	 *  potential atoms on the second component.
	 * @param labelOfAtomOnComponent2 label of the atom on the second
	 *  component.
	 */
	public ModificationLinkage(
			final List<Component> components,
			final int indexOfComponent1,
			final List<String> pdbNameOfPotentialAtomsOnComponent1,
			final String labelOfAtomOnComponent1,
			final int indexOfComponent2,
			final List<String> pdbNameOfPotentialAtomsOnComponent2,
			final String labelOfAtomOnComponent2) {
		if (components == null) {
			throw new IllegalArgumentException("Null components");
		}
		
		if ( indexOfComponent1 < 0)
			throw new IllegalArgumentException("indexOfComponent1 has to be >= 0");
		if ( indexOfComponent1 >= components.size())
			throw new IllegalArgumentException("indexOfComponent1 has to be <= components.size()");
		
		if ( indexOfComponent2 < 0)
			throw new IllegalArgumentException("indexOfComponent2 has to be >= 0");
		if ( indexOfComponent2 >= components.size())
			throw new IllegalArgumentException("indexOfComponent2 [" + indexOfComponent2 + "] has to be <= components.size() [" + components.size()+"]");
		
		
		
		if (indexOfComponent1 == indexOfComponent2) {
			throw new IllegalArgumentException("No linkage is allowed for an" +
					" identical component.");
		}
		
		this.components = components;
		this.indexOfComponent1 = indexOfComponent1;
		this.indexOfComponent2 = indexOfComponent2;
		this.pdbNameOfPotentialAtomsOnComponent1 = pdbNameOfPotentialAtomsOnComponent1;
		this.pdbNameOfPotentialAtomsOnComponent2 = pdbNameOfPotentialAtomsOnComponent2;
		this.labelOfAtomOnComponent1 = labelOfAtomOnComponent1;
		this.labelOfAtomOnComponent2 = labelOfAtomOnComponent2;
	}

	/**
	 * 
	 * @return index of the first component.
	 */
	public int getIndexOfComponent1() {
		return indexOfComponent1;
	}
	
	/**
	 * 
	 * @return index of the second component.
	 */
	public int getIndexOfComponent2() {
		return indexOfComponent2;
	}
	
	/**
	 * 
	 * @return the first component.
	 */
	public Component getComponent1() {
		return components.get(indexOfComponent1);
	}
	
	/**
	 * 
	 * @return the second component.
	 */
	public Component getComponent2() {
		return components.get(indexOfComponent2);
	}
	
	/**
	 * 
	 * @return a list of PDB names of potential atoms on the first component.
	 */
	public List<String> getPDBNameOfPotentialAtomsOnComponent1() {
		return pdbNameOfPotentialAtomsOnComponent1;
	}
	
	/**
	 * 
	 * @return a list of PDB names of potential atoms on the second component.
	 */
	public List<String> getPDBNameOfPotentialAtomsOnComponent2() {
		return pdbNameOfPotentialAtomsOnComponent2;
	}
	
	/**
	 * 
	 * @return label of the atom on the first component.
	 */
	public String getLabelOfAtomOnComponent1() {
		return labelOfAtomOnComponent1;
	}
	
	/**
	 * 
	 * @return label of the atom on the second component.
	 */
	public String getLabelOfAtomOnComponent2() {
		return labelOfAtomOnComponent2;
	}
	
	/**
	 * 
	 */
	@Override
	public String toString() {
		Component comp1 = getComponent1();
		Component comp2 = getComponent2();
		List<String> atom1 = getPDBNameOfPotentialAtomsOnComponent1();
		List<String> atom2 = getPDBNameOfPotentialAtomsOnComponent2();
		if ( comp1 == null || comp2 == null) {
			return "ModificationLinkage: empty";
		}
		if ( comp1.getPdbccIds() != null && comp2.getPdbccIds() != null) {
			return "ModificationLinkage: " + comp1.getPdbccIds().toString()+":"+atom1+"<=>" + comp2.getPdbccIds()+atom2;
		} else {
			return "ModificationLinkage :"+atom1+"<=>" + atom2;
		}
	}
}
