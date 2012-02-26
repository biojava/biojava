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
 * DNA Sequences produced by modern sequencers usually have quality informaion
 * attached to them. This feature allows to store the information directly in
 * the DNASequence
 *
 * @since 3.0.3
 * @author brandstaetter
 */
public class QualityFeature<S extends AbstractSequence<C>, C extends Compound> extends AbstractFeature<S, C> {

    private List<Number> qualities = new ArrayList<Number>();

    /**
     *
     * @param type
     * @param source
     */
    public QualityFeature(String type, String source) {
        super(type, source);
    }

    /**
     * @return the qualities
     */
    public List<Number> getQualities() {
        return qualities;
    }

    /**
     * @param qualities the qualities to set
     */
    public void setQualities(List<Number> qualities) {
        this.qualities = qualities;
    }

    /**
     *
     * @param bioindex
     * @return
     */
    public Number getQualityAt(int bioindex) {
        return qualities.get(bioindex - 1);
    }
}
