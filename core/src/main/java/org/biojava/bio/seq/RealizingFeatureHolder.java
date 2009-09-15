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

package org.biojava.bio.seq;

import org.biojava.bio.BioException;

/**
 * Interface for <code>FeatureHolder</code> objects which know how to
 * instantiate new child Features.  This interface should not be
 * needed in normal client programs, since they will use the
 * <code>createFeature</code> method of <code>FeatureHolder</code> to
 * add new features.  However, this method exposes the feature
 * realization infrastructure to child features.  </p>
 *
 * @see org.biojavax.bio.seq.RichFeatureRelationshipHolder
 * @author Thomas Down
 */

public interface RealizingFeatureHolder extends FeatureHolder {
    /**
     * Realize a feature template.  This will be a template which has
     * been passed to the <code>createFeature</code> method of either
     * this <code>FeatureHolder</code> or one of our child Features. 
     */

    public Feature realizeFeature(FeatureHolder parent,
				  Feature.Template template)
	throws BioException;
}
