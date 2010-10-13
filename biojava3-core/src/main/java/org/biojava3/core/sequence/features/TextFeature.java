/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava3.core.sequence.features;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class TextFeature extends AbstractFeature {

    public TextFeature(String type, String source,String shortDescription,String description) {
        super(type, source);
        this.setDescription(description);
        this.setShortDescription(shortDescription);
    }
}
