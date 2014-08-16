/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava3.survival.cox.comparators;

import org.biojava3.survival.cox.CoxInfo;
import org.biojava3.survival.cox.CoxVariables;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class CoxVariablesVariableComparator implements CoxComparatorInterface {
 
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
        description = "Signatures ranked by model " + variables + " and variable " + variable + " p-value";
    }

    @Override
    public int compare(CoxVariables coxVariables1, CoxVariables coxVariables2) {
        if(coxVariables1.equals(coxVariables2))
            return 0;
        CoxInfo ci1 = coxVariables1.getCoxInfo(variables);
        CoxInfo ci2 = coxVariables2.getCoxInfo(variables);
        if(ci1 == null && ci2 == null)
            return 0;
        if(ci1 == null && ci2 != null)
            return 1;
        if(ci2 == null && ci1 != null)
            return -1;
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

    @Override
    public String getDescription() {
       return description;
    }
    String description = "Signatures ranked by p-value ";

    @Override
    public void setDescription(String description) {
       this.description = description;
    }

    @Override
    public String getModelVariables() {
        return variables;
    }

    @Override
    public String getSortVariable() {
        return variable;
    }


    
}
