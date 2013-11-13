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
public class CoxVariablesOverallModelFitComparator implements Comparator<CoxVariables> {

    String variables = "";

    /**
     *
     * @param variables
     */
    public CoxVariablesOverallModelFitComparator(String variables) {
        this.variables = variables;
    }

    public int compare(CoxVariables coxVariables1, CoxVariables coxVariables2) {
        CoxInfo ci1 = coxVariables1.getCoxInfo(variables);
        CoxInfo ci2 = coxVariables2.getCoxInfo(variables);

        if (ci1.getScoreLogrankTestpvalue() < ci2.getScoreLogrankTestpvalue()) {
            return -1;
        } else if (ci1.getScoreLogrankTestpvalue() > ci2.getScoreLogrankTestpvalue()) {
            return 1;
        } else {
            return 0;
        }
        //ascending order
        // return coxVariables1.compareTo(coxVariables2);
    }
}
