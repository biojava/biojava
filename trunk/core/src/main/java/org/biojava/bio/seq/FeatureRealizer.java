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
 * Interface for translators which map from Feature.Template
 * instances to real Feature objects. <strong>NOTE</strong>
 * this interface is intended for use primarily by Sequence
 * implementors.  Normal users can generally ignore it.
 *
 * <p>
 * There is no requirement that Sequence implementations
 * use FeatureRealizers to construct their features, but common
 * implementations such as SimpleSequence and ViewSequence do.
 * </p>
 *
 * @author Thomas Down
 */

public interface FeatureRealizer {
    /**
     * Install a feature on the specified sequence.
     *
     * @param seq The sequence to which the feature will be added.
     * @param template A description of the desired feature.
     * @param parent The FeatureHolder which is to be the Feature's
     *               immediate parent.
     * @return A newly constructed feature, to be added to the sequence.
     * @throws BioException If the feature could not be constructed.
     */

    public Feature realizeFeature(Sequence seq,
				  FeatureHolder parent,
				  Feature.Template template)
            throws BioException;
}


