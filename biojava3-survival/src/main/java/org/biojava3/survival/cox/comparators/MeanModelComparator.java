/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava3.survival.cox.comparators;

import org.biojava3.survival.cox.CoxInfo;
import org.biojava3.survival.cox.CoxVariables;
import java.util.Comparator;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class MeanModelComparator implements Comparator<CoxVariables> {

    String variable = "";

    /**
     *
     * @param variable
     */
    public MeanModelComparator(String variable) {
      this.variable = variable;
    }

    public int compare(CoxVariables coxVariables1, CoxVariables coxVariables2) {
        CoxInfo ci1LTMean = coxVariables1.getCoxInfo("<MEAN");
        CoxInfo ci1GTMean = coxVariables1.getCoxInfo(">MEAN");
        
        if(ci1LTMean == null || ci1GTMean == null)
            return 0;
        
        double c1LTpvalue = ci1LTMean.getCoefficient(variable).getPvalue();
        double c1GTpvalue = ci1GTMean.getCoefficient(variable).getPvalue();
        
        double c1ratio = Math.min(c1LTpvalue, c1GTpvalue) / Math.max(c1LTpvalue, c1GTpvalue);
        
        CoxInfo ci2LTMean = coxVariables2.getCoxInfo("<MEAN");
        CoxInfo ci2GTMean = coxVariables2.getCoxInfo(">MEAN");
        
        double c2LTpvalue = ci2LTMean.getCoefficient(variable).getPvalue();
        double c2GTpvalue = ci2GTMean.getCoefficient(variable).getPvalue();
        
       double c2ratio = Math.min(c2LTpvalue, c2GTpvalue) / Math.max(c2LTpvalue, c2GTpvalue);

        if (c1ratio > c2ratio) {
            return 1;
        } else if (c1ratio < c2ratio) {
            return -1;
        } else {
            return 0;
        }
        //ascending order
        // return coxVariables1.compareTo(coxVariables2);
    }
}
