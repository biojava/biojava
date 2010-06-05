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

import java.util.HashMap;
import java.util.Map;

/**
 * indicates the type of a component involved in a protein
 * modification, e.g., an amino acid or a ligand.
 */
public enum ComponentType {
	AMINOACID("AminoAcid"),
	LIGAND("Ligand")
	;
	
	ComponentType(String label) {
		this.label = label;
	}
	
	/**
	 * 
	 * @return the label of this ModificationCategory.
	 */
	public String label() {
		return label;
	}
	
	/**
	 * @return the label of this ModificationCategory.
	 */
	public String toString() {
		return label;
	}
	
	/**
	 * The variable is the same as the &ltType&gt; in the ptm_list XML file.
	 */
	private String label;
	
	/**
	 * 
	 * @param label the label of ModificationCategory.
	 * @return the ModificationCategory that has the label.
	 */
	public static ComponentType getByLabel(String label) {
		return mapLabelType.get(label);
	}
	
	private static Map<String, ComponentType> mapLabelType;	
	static {
		mapLabelType = new HashMap<String, ComponentType>();
		for (ComponentType type:ComponentType.values()) {
			mapLabelType.put(type.label, type);
		}
	}
}
