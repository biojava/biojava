/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava3.core.sequence.features;

import org.biojava3.core.sequence.template.AbstractSequence;
import org.biojava3.core.sequence.template.Compound;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class TextFeature<S extends AbstractSequence<C>, C extends Compound> extends AbstractFeature<S, C> {

    public TextFeature(String type, String source,String shortDescription,String description) {
        super(type, source);
        this.setDescription(description);
        this.setShortDescription(shortDescription);
    }
}
