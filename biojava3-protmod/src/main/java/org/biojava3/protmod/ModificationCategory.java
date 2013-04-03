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
 * Created on May 30, 2010
 * Author: Jianjiong Gao 
 *
 */

package org.biojava3.protmod;

import java.util.HashMap;
import java.util.Map;

/**
 * define modification categories.
 * 
 * @author Jianjiong Gao
 * @since 3.0
 */
public enum ModificationCategory {
	ATTACHMENT("attachment", "Attachment of a chemical group"),
	CHEMICAL_MODIFICATION("modified residue", "Chemical modification of an residue"), 
	CROSS_LINK_1("crosslink1", "A complicated cross-link (usually found in chromophores)"),
	CROSS_LINK_2("crosslink2", "A cross-link of two residues"),
	CROSS_LINK_3("crosslink3", "A cross-link of three residues"),
	CROSS_LINK_4("crosslink4", "A cross-link of four residues"),
	CROSS_LINK_5("crosslink5", "A cross-link of five residues"),
	CROSS_LINK_6("crosslink6", "A cross-link of six residues"),
	CROSS_LINK_7("crosslink7", "A cross-link of seven residues"),
	CROSS_LINK_8_OR_LARGE("crosslink8 or large", "A cross-link of eight or more residues"), // 8 or high
	UNDEFINED("undefined", "Undefined category")
	;
	
	ModificationCategory(String label, String desc) {
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
	 * 
	 * @return the description
	 */
	public String description() {
		return desc;
	}
	
	/**
	 * @return the label of this ModificationCategory.
	 */
	public String toString() {
		return label;
	}
	
	/**
	 * 
	 * @return true if it is a CrossLink; false, otherwise.
	 */
	public boolean isCrossLink() {
		return this == CROSS_LINK_1
			|| this == CROSS_LINK_2
			|| this == CROSS_LINK_3
			|| this == CROSS_LINK_4
			|| this == CROSS_LINK_5
			|| this == CROSS_LINK_6
			|| this == CROSS_LINK_7
			|| this == CROSS_LINK_8_OR_LARGE;
	}
	
	/**
	 * The variable is the same as the &ltType&gt; in the ptm_list XML file.
	 */
	private String label;
	
	private String desc;
	
	/**
	 * 
	 * @param label the label of ModificationCategory.
	 * @return the ModificationCategory that has the label.
	 */
	public static ModificationCategory getByLabel(String label) {
		return mapLabelCat.get(label);
	}
	
	private static Map<String, ModificationCategory> mapLabelCat;	
	static {
		mapLabelCat = new HashMap<String, ModificationCategory>();
		for (ModificationCategory cat:ModificationCategory.values()) {
			mapLabelCat.put(cat.label, cat);
		}
	}	
}
