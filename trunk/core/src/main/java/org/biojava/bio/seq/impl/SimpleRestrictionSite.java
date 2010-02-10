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

import org.biojava.bio.molbio.RestrictionEnzyme;
import org.biojava.bio.molbio.RestrictionSite;
import org.biojava.bio.seq.FeatureHolder;
import org.biojava.bio.seq.Sequence;

/**
 * <code>SimpleRestrictionSite</code> represents the recognition site
 * of a restriction enzyme.
 *
 * @author Keith James
 * @since 1.3
 */
public class SimpleRestrictionSite extends SimpleStrandedFeature
    implements RestrictionSite
{
    private RestrictionEnzyme enzyme;
    private int position;

    /**
     * Creates a new <code>SimpleRestrictionSite</code>.
     *
     * @param sourceSeq a <code>Sequence</code>.
     * @param parent a <code>FeatureHolder</code>.
     * @param template a <code>RestrictionSite.Template</code>.
     */
    public SimpleRestrictionSite(Sequence                 sourceSeq,
                                 FeatureHolder            parent,
                                 RestrictionSite.Template template)
    {
        super(sourceSeq, parent, template);

        this.enzyme = template.enzyme;

        position = template.location.getMin() +
            template.enzyme.getDownstreamCut()[0];
    }

    public int getPosition()
    {
        return position;
    }

    public RestrictionEnzyme getEnzyme()
    {
        return enzyme;
    }

    public String toString()
    {
        return super.toString() + " " + enzyme;
    }
}
