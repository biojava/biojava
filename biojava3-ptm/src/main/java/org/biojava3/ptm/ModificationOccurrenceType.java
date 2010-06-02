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

public enum ModificationOccurrenceType {
	NATURAL("Natural");
	
	ModificationOccurrenceType(String label) {
		this.label = label;
	}
	
	/**
	 * 
	 * @return the label of this ModificationOccurrenceType.
	 */
	public String label() {
		return label;
	}
	
	/**
	 * @return the label of this ModificationOccurrenceType.
	 */
	public String toString() {
		return label;
	}
	
	/**
	 * The variable is the same as the &ltOccurrence&gt; in the ptm_list XML file.
	 */
	private String label;
	
	/**
	 * 
	 * @param label the label of ModificationOccurrenceType.
	 * @return the ModificationOccurrenceType that has the label.
	 */
	public static ModificationOccurrenceType getByLabel(String label) {
		return mapLabelOcc.get(label);
	}
	
	private static Map<String, ModificationOccurrenceType> mapLabelOcc;	
	static {
		mapLabelOcc = new HashMap<String, ModificationOccurrenceType>();
		for (ModificationOccurrenceType cat:ModificationOccurrenceType.values()) {
			mapLabelOcc.put(cat.label, cat);
		}
	}
}
