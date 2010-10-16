/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.biojava3.core.sequence.features;

import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class QuantityFeature extends AbstractFeature {

    private List<Number> quantities = new ArrayList<Number>();


    public QuantityFeature(String type,String source){
       super(type,source);
    }

    public void addQuantity(Number value){
        quantities.add(value);
    }
    
    /**
     * @return the quantities
     */
    public List<Number> getQuantities() {
        return quantities;
    }

    /**
     * @param quantities the quantities to set
     */
    public void setQuantities(List<Number> quantities) {
        this.quantities = quantities;
    }



}
