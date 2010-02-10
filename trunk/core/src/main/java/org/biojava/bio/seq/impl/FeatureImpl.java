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

import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.bio.molbio.RestrictionSite;
import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.FeatureRealizer;
import org.biojava.bio.seq.FramedFeature;
import org.biojava.bio.seq.RemoteFeature;
import org.biojava.bio.seq.SimpleFeatureRealizer;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.seq.homol.HomologyFeature;
import org.biojava.bio.seq.homol.SimilarityPairFeature;
import org.biojava.utils.StaticMemberPlaceHolder;

/**
 * Wrap up default sets of Feature implementations.
 *
 * @author Thomas Down
 * @author Greg Cox
 * @author Keith James
 * @since 1.1
 * @see org.biojavax.bio.seq.SimpleRichFeature
 */

public class FeatureImpl {
    /**
     * Default implementation of FeatureRealizer, which wraps simple
     * implementations of Feature and StrandedFeature.  This is the
     * default FeatureRealizer used by SimpleSequence and ViewSequence,
     * and may also be used by others.  When building new FeatureRealizers,
     * you may wish to use this as a `fallback' realizer, and benefit from
     * the Feature and StrandedFeature implementations.
     */

    public final static FeatureRealizer DEFAULT;

    static {
        SimpleFeatureRealizer d  = new SimpleFeatureRealizer() {
            public Object writeReplace() {
                try {
                    return new StaticMemberPlaceHolder(FeatureImpl.class.getField("DEFAULT"));
                } catch (NoSuchFieldException ex) {
                    throw new BioError(ex);
                }
            }
        } ;

        try {
            d.addImplementation(Feature.Template.class,
                                SimpleFeature.class);
            d.addImplementation(StrandedFeature.Template.class,
                                SimpleStrandedFeature.class);
            d.addImplementation(HomologyFeature.Template.class,
                                SimpleHomologyFeature.class);
            d.addImplementation(SimilarityPairFeature.Template.class,
                                SimpleSimilarityPairFeature.class);
            d.addImplementation(RemoteFeature.Template.class,
                                SimpleRemoteFeature.class);
            d.addImplementation(FramedFeature.Template.class,
                                SimpleFramedFeature.class);
            d.addImplementation(RestrictionSite.Template.class,
                                SimpleRestrictionSite.class);
        } catch (BioException ex) {
            throw new BioError("Couldn't initialize default FeatureRealizer", ex);
        }

        DEFAULT = d;
    }
}
