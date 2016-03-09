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

import org.apache.commons.math.stat.descriptive.DescriptiveStatistics;
import org.biojava.nbio.survival.data.WorkSheet;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class MeanQuantizer implements DiscreteQuantizerInterface {

	/**
	 *
	 * @param worksheet
	 * @param column
	 */
	@Override
	public void process(WorkSheet worksheet, String column) {
		DescriptiveStatistics ds = new DescriptiveStatistics();
		for (String row : worksheet.getRows()) {
			try {
				Double d = Double.parseDouble(worksheet.getCell(row, column));
				ds.addValue(d);
			} catch (Exception e) {
			}
		}
		Double mean = ds.getMean();
		for (String row : worksheet.getRows()) {
			try {
				Double d = Double.parseDouble(worksheet.getCell(row, column));
				if (d < mean) {
					worksheet.addCell(row, column, "L");
				} else {
					worksheet.addCell(row, column, "H");
				}
			} catch (Exception e) {
				try {
					worksheet.addCell(row, column, "");
				} catch (Exception e1) {
				}
			}
		}
	}
}
