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
 
import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.FeatureHolder;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.homol.Homology;
import org.biojava.bio.seq.homol.HomologyFeature;
import org.biojava.bio.symbol.IllegalAlphabetException;

/**
 * @author Matthew Pocock
 * @author Keith James
 */
public class SimpleHomologyFeature
    extends SimpleStrandedFeature
    implements HomologyFeature {
    private Homology homology;
  
    public Homology getHomology() {
        return this.homology;
    }
  
    public Feature.Template makeTemplate() {
        HomologyFeature.Template ft = new HomologyFeature.Template();
        fillTemplate(ft);
        return ft;
    }
  
    protected void fillTemplate(HomologyFeature.Template ft) {
        super.fillTemplate(ft);
        ft.homology = getHomology();
    }
  
    public SimpleHomologyFeature(Sequence sourceSeq,
                                 FeatureHolder parent,
				 HomologyFeature.Template template)
	throws IllegalArgumentException, IllegalAlphabetException 
    {
	super(sourceSeq, parent, template);
        this.homology = template.homology;
    }
    
    public String toString() {
        return super.toString() + " " + getHomology();
    }
}
