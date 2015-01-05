/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava3.survival.cox.comparators;


import java.util.Comparator;
import org.biojava3.survival.cox.CoxVariables;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public interface CoxComparatorInterface extends Comparator<CoxVariables> {
    public String getDescription();
    public void setDescription(String description);
    public String getModelVariables();
    public String getSortVariable(); 
}
