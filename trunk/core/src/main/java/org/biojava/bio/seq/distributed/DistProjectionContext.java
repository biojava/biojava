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

package org.biojava.bio.seq.distributed;

import java.util.HashMap;
import java.util.Map;

import org.biojava.bio.Annotation;
import org.biojava.bio.BioError;
import org.biojava.bio.BioRuntimeException;
import org.biojava.bio.MergeAnnotation;
import org.biojava.bio.seq.ComponentFeature;
import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.FeatureHolder;
import org.biojava.bio.seq.projection.ReparentContext;
import org.biojava.utils.ChangeVetoException;

/**
 * Projection for MetaDAS.
 *
 * <p>
 * In this new version, most of the functionality is inherited from the normal
 * ProjectedFeatureHolder
 * </p>
 *
 * @author Thomas Down
 * @since 1.2
 */


class DistProjectionContext
        extends ReparentContext {
    private Annotation annotation;
    private Map componentFeatureCache = new HashMap();

    public DistProjectionContext(FeatureHolder fh,
                                 FeatureHolder parent,
                                 Annotation annotation)
    {
        super(parent, fh);
        this.annotation = annotation;
    }

    public Feature projectFeature(Feature f) {
            if (f instanceof ComponentFeature && getParent() instanceof DistributedSequence) {
            ComponentFeature pcf = (ComponentFeature) componentFeatureCache.get(f);
            if (pcf == null) {
                ComponentFeature.Template cft = (ComponentFeature.Template) ((ComponentFeature) f).makeTemplate();

                if (cft.componentSequenceName == null) {
                    cft.componentSequenceName = cft.componentSequence.getName();
                }
                if (cft.componentSequenceName == null) {
                    throw new NullPointerException("Can't get component sequence name");
                }

                cft.componentSequence = null;    // We need to go back though the DistDB for the
                                                         // proper component sequence to use here.

                try {
                    pcf = new DistComponentFeature((DistributedSequence) getParent(),
                                                                           cft);
                } catch (Exception ex) {
                    throw new BioRuntimeException("Error instantiating DistComponentFeature for: " + f, ex);
                }
                componentFeatureCache.put(f, pcf);
            }
            return pcf;
            } else {
            // Default: generate a throwaway ProjectedFeature

            return super.projectFeature(f);
            }
        }

        public Annotation getAnnotation(Feature f) {
        if (annotation != null) {
            try {
                MergeAnnotation ma = new MergeAnnotation();
                ma.addAnnotation(f.getAnnotation());
                ma.addAnnotation(annotation);
                return ma;
            } catch (ChangeVetoException cve) {
                throw new BioError(cve);
            }
        } else {
            return f.getAnnotation();
        }
        }
}
