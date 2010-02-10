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
 */

package org.biojava.bio.seq.impl;

import java.io.Serializable;

import org.biojava.bio.Annotation;
import org.biojava.bio.seq.FeatureRealizer;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceFactory;
import org.biojava.bio.symbol.SymbolList;

/**
 * A no-frills implementation of SequenceFactory that produces SimpleSequence
 * objects.
 *
 * @author Matthew Pocock
 * @author Thomas Down
 */
public class SimpleSequenceFactory implements SequenceFactory, Serializable {
    private FeatureRealizer realizer = org.biojava.bio.seq.impl.FeatureImpl.DEFAULT;
    
    /**
    *Returns the FeatureRealizer set by "setFeatureRealizer".
    */
    
    public FeatureRealizer getFeatureRealizer() {
	return realizer;
    }

    /**
     * Set the FeatureRealizer used by new sequences created by this
     * factory.
     */

    public void setFeatureRealizer(FeatureRealizer fr) {
	realizer = fr;
    }

    public Sequence createSequence(SymbolList symList,
				   String uri, String name, Annotation annotation) {
	return new org.biojava.bio.seq.impl.SimpleSequence(symList, uri, name, annotation, realizer);
    }
}
