/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava3.survival.cox.comparators;

import org.biojava3.survival.cox.SurvivalInfo;
import java.util.Comparator;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class SurvivalInfoValueComparator implements Comparator<SurvivalInfo> {

    String variable = "";
    
    /**
     *
     * @param variable
     */
    public SurvivalInfoValueComparator(String variable){
        this.variable = variable;
    }
    
    public int compare(SurvivalInfo t, SurvivalInfo t1) {
        double v = t.getContinuousVariable(variable);
        double v1 = t1.getContinuousVariable(variable);
        if(v < v1)
            return -1;
        else if(v > v1)
            return 1;
        else
            return 0;
    }
    
}
