/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava3.survival.kaplanmeier.metadata;

import org.biojava3.survival.data.WorkSheet;
import org.apache.commons.math.stat.descriptive.DescriptiveStatistics;

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
