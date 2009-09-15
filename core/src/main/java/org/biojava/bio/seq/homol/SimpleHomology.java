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

package org.biojava.bio.seq.homol;
 
import java.util.Iterator;

import org.biojava.bio.BioException;
import org.biojava.bio.alignment.Alignment;
import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.FeatureHolder;
import org.biojava.bio.seq.SimpleFeatureHolder;
import org.biojava.utils.AbstractChangeable;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeVetoException;

/**
 * A no-frills implementation of Homology.
 *
 * @author Matthew Pocock
 * @author <a href="mailto:kdj@sanger.ac.uk">Keith James</a>
 * @since 1.2
 */
public class SimpleHomology extends AbstractChangeable implements Homology
{
    private SimpleFeatureHolder features;
    private Alignment           alignment;

    /**
     * Creates a new empty <code>SimpleHomology</code> containing no
     * <code>Alignment</code> and no <code>FeatureHolder</code>.
     */
    public SimpleHomology() { }

    /**
     * <code>getFeatures</code> returns the constituent
     * <code>HomologyFeature</code>s which are also used as the keys
     * in the alignment.
     *
     * @return a <code>FeatureHolder</code>.
     */
    public FeatureHolder getFeatures()
    {
        return features;
    }

    /**
     * <code>getAlignment</code> returns the alignment, which uses the
     * <code>HomologyFeature</code>s as keys.
     *
     * @return an <code>Alignment</code>.
     */
    public Alignment getAlignment()
    {
        return alignment;
    }

    /**
     * <code>setAlignment</code> sets the alignment which describes
     * the homology. The alignment, should use the
     * <code>HomologyFeature</code>s as keys. A suitable
     * <code>FeatureHolder</code> is automatically created.
     *
     * @param alignment an <code>Alignment</code>.
     *
     * @exception BioException if an error occurs.
     * @exception ChangeVetoException if the
     * <code>SimpleHomology</code> is locked.
     */
    public void setAlignment(Alignment alignment)
        throws BioException, ChangeVetoException
    {
        if (hasListeners())
        {
            ChangeSupport cs = getChangeSupport(Homology.ALIGNMENT);

            ChangeEvent ce = new ChangeEvent(this, Homology.ALIGNMENT,
                                             alignment, this.alignment);

            synchronized(cs)
            {
              cs.firePreChangeEvent(ce);
              this.alignment = alignment;
              cs.firePostChangeEvent(ce);
            }  
        }
        else
        {
            this.alignment = alignment;
        }

        features = new SimpleFeatureHolder();

        for (Iterator li = alignment.getLabels().iterator(); li.hasNext();)
        {
            Object o = li.next();
            if (! HomologyFeature.class.isInstance(o))
                throw new BioException("The labels of the Alignment used to construct a SimpleHomology should be the relevant HomologyFeatures");

            features.addFeature((Feature) o);
        }
    }

    public String toString()
    {
        return "SimpleHomology [" + alignment.getLabels().size() + " features]";
    }
}
