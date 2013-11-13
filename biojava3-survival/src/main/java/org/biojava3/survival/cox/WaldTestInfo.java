/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava3.survival.cox;

import org.biojava3.survival.cox.stats.ChiSq;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class WaldTestInfo {

    private int df = 0;
    double[][] solve;
    double[] bsum;

    /**
     *
     * @return
     */
    public double getTest() {
        return bsum[0];
    }

    @Override
    public String toString() {
        return "Wald test=" + getTest() + "df=" + df + " p-value=" + getPvalue(); //To change body of generated methods, choose Tools | Templates.
    }

    /**
     *
     * @return
     */
    public double getPvalue() {
        double pvalue = ChiSq.chiSq(getTest(), df);
        return pvalue;
    }

    /**
     * @return the df
     */
    public int getDf() {
        return df;
    }

    /**
     * @param df the df to set
     */
    public void setDf(int df) {
        this.df = df;
    }
}
