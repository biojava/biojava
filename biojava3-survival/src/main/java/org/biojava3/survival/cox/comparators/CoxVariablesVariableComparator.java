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
public class CoxVariablesVariableComparator implements Comparator<CoxVariables> {

    String variables = "";
    String variable = "";

    /**
     *
     * @param variables
     * @param variable
     */
    public CoxVariablesVariableComparator(String variables, String variable) {
        this.variables = variables;
        this.variable = variable;
    }

    public int compare(CoxVariables coxVariables1, CoxVariables coxVariables2) {
        CoxInfo ci1 = coxVariables1.getCoxInfo(variables);
        CoxInfo ci2 = coxVariables2.getCoxInfo(variables);
        if(ci1 == null || ci2 == null)
            return 0;
        if (ci1.getCoefficientsList().get(variable).getPvalue() < ci2.getCoefficientsList().get(variable).getPvalue()) {
            return -1;
        } else if (ci1.getCoefficientsList().get(variable).getPvalue() > ci2.getCoefficientsList().get(variable).getPvalue()) {
            return 1;
        } else {
            return 0;
        }
        //ascending order
        // return coxVariables1.compareTo(coxVariables2);
    }
}
