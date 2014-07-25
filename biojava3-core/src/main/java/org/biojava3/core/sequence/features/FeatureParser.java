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
 * @author Jacek Grzebyta <github:jgrzebyta>
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 * Created on 01-07-2014
 *
 */

package org.biojava3.core.sequence.features;

import org.biojava3.core.sequence.template.AbstractSequence;
import org.biojava3.core.sequence.template.Compound;

/**
 * @author jgrzebyta
 * @param <C>
 */
public interface FeatureParser<C extends Compound> {

    FeatureInterface<AbstractSequence<C>, C> getFeature(String keyword);

    void parseFeatures(AbstractSequence<C> sequence);
    
}
