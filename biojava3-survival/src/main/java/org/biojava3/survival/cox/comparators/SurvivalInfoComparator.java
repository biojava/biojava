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
public class SurvivalInfoComparator implements Comparator<SurvivalInfo> {

    public int compare(SurvivalInfo t, SurvivalInfo t1) {
        return t.getOrder() - t1.getOrder();
    }
    
}
