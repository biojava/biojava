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

package org.biojava3.ptm;

import java.util.HashMap;
import java.util.Map;

/**
 * define modification categories.
 */
public enum ModificationCategory {
	ATTACHMENT("Attachment"),
	CHEMICAL_MODIFICATION("ModifiedResidue"), 
	CROSS_OVER1("CrossOver1"),
	CROSS_OVER2("CrossOver2"),
	CROSS_OVER3("CrossOver3"),
	CROSS_OVER4("CrossOver4"),
	CROSS_OVER5("CrossOver5"),
	CROSS_OVER6("CrossOver6"),
	CROSS_OVER7("CrossOver7")
	;
	
	ModificationCategory(String label) {
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
