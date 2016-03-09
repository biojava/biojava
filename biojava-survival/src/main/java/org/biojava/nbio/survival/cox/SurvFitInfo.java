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
 */
package org.biojava.nbio.survival.cox;

import java.util.LinkedHashMap;

/**
 * Contains info for graphing km figures
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class SurvFitInfo {

	private LinkedHashMap<String, StrataInfo> strataInfoHashMap = new LinkedHashMap<String, StrataInfo>();
	private LinkedHashMap<String, StrataInfo> unweightedStrataInfoHashMap = new LinkedHashMap<String, StrataInfo>();
	private boolean weighted = false;


	/**
	 *
	 * @return
	 */
	public LinkedHashMap<String, StrataInfo> getUnweightedStrataInfoHashMap() {
		return unweightedStrataInfoHashMap;
	}

	/**
	 *
	 * @param unweightedStrataInfoHashMap
	 */
	public void setUnweightedStrataInfoHashMap(LinkedHashMap<String, StrataInfo> unweightedStrataInfoHashMap) {
		this.unweightedStrataInfoHashMap = unweightedStrataInfoHashMap;
	}

	/**
	 * @return the strataInfoHashMap
	 */
	public LinkedHashMap<String, StrataInfo> getStrataInfoHashMap() {
		return strataInfoHashMap;
	}

	/**
	 * @param strataInfoHashMap the strataInfoHashMap to set
	 */
	public void setStrataInfoHashMap(LinkedHashMap<String, StrataInfo> strataInfoHashMap) {
		this.strataInfoHashMap = strataInfoHashMap;
	}

	/**
	 *
	 * @param siHashMap
	 * @param label
	 */
	public void addStrataInfoHashMap(LinkedHashMap<String, StrataInfo> siHashMap, String label) {
		for (String key : siHashMap.keySet()) {
			StrataInfo si = siHashMap.get(key);
			strataInfoHashMap.put(label + " " + key, si);
		}
	}

	@Override
	public String toString() {
		return strataInfoHashMap.toString(); //To change body of generated methods, choose Tools | Templates.
	}

	/**
	 * @return the weighted
	 */
	public boolean isWeighted() {
		return weighted;
	}

	/**
	 * @param weighted the weighted to set
	 */
	public void setWeighted(boolean weighted) {
		this.weighted = weighted;
	}
}
