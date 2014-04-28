/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava3.survival.kaplanmeier.metadata;

import org.biojava3.survival.data.WorkSheet;
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
    ArrayList<String> discreteValues = new ArrayList<String>();

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
