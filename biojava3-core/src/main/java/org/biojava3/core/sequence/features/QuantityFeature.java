/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 * Created on 01-21-2010
 */
package org.biojava3.core.sequence.features;

import java.util.ArrayList;
import java.util.List;

import org.biojava3.core.sequence.template.AbstractSequence;
import org.biojava3.core.sequence.template.Compound;

/**
 * It is common to have a numerical value or values associated with a feature. This can then
 * be used in heat maps or other visual indicators when viewing a sequence. Multiple quantities
 * could represent a time corse study and display a color gradient
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class QuantityFeature<S extends AbstractSequence<C>, C extends Compound> extends AbstractFeature<S, C> {

    private List<Number> quantities = new ArrayList<Number>();

    /**
     *
     * @param type
     * @param source
     */
    public QuantityFeature(String type, String source) {
        super(type, source);
    }

    /**
     *
     * @param value
     */
    public void addQuantity(Number value) {
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
