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
package org.biojava.nbio.survival.kaplanmeier.metadata;

import org.biojava.nbio.survival.data.WorkSheet;

import java.util.ArrayList;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class MetaDataInfo {

	/**
	 *
	 */
	public String column = "";
	/**
	 *
	 */
	public boolean numeric = false;
	/**
	 *
	 */
	public DiscreteQuantizerInterface discreteQuantizer = null;
	ArrayList<String> discreteValues = new ArrayList<>();

	/**
	 *
	 * @param column
	 * @param numeric
	 * @param discreteQuantizer
	 */
	public MetaDataInfo(String column, boolean numeric, DiscreteQuantizerInterface discreteQuantizer) {
		this.column = column;
		this.numeric = numeric;
		this.discreteQuantizer = discreteQuantizer;
	}

	/**
	 *
	 * @param column
	 */
	public MetaDataInfo(String column) {
		this.column = column;
	}

	/**
	 *
	 * @param worksheet
	 * @throws Exception
	 */
	public void setDiscreteValues(WorkSheet worksheet) throws Exception {
		discreteValues = worksheet.getDiscreteColumnValues(column);
	}

	/**
	 *
	 * @return
	 */
	public ArrayList<String> getDiscreteValues() {
		return discreteValues;
	}

	/**
	 *
	 * @return
	 */
	public int getNumberDiscreteValues() {
		return discreteValues.size();
	}
}
