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

import java.util.Iterator;

import org.biojava.bio.BioException;
import org.biojava.utils.AbstractChangeable;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeListener;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;

/**
 * Wrapper implementation of FeatureHolder which calls a method
 * to create a contained FeatureHolder on demand.  This is an
 * abstract class and is normally used like:
 *
 * FeatureHolder fh = new LazyFeatureHolder() {
 *         protected FeatureHolder createFeatureHolder() {
 *             SimpleFeatureHolder features = new SimpleFeatureHolder();
 *             // Create some features in here...
 *             return features;
 *         }
 *     } ;
 * </pre>
 *
 * @author Thomas Down
 * @author Matthew Pocock
 * @since 1.2
 * @see org.biojavax.bio.seq.RichFeatureRelationshipHolder 
 */

public abstract class LazyFeatureHolder
  extends
    AbstractChangeable
  implements
    FeatureHolder
{
    private FeatureHolder featureHolder;
    private Forwarder changeForwarder;
    private FeatureFilter schema;
    
    /**
     * Construct a LazyFeatureHolder with the schema of its contained featureholder
     */
    
    protected LazyFeatureHolder() {
    }
    
    /**
     * Construct a LazyFeatureHolder with the specified schema
     */
     
    protected LazyFeatureHolder(FeatureFilter schema) {
        this.schema = schema;
    }
    
    public FeatureFilter getSchema() {
        if (schema == null) {
            return getFeatureHolder().getSchema();
        } else {
            return schema;
        }
    }

    protected abstract FeatureHolder createFeatureHolder();

    protected void flushFeatures() {
        featureHolder = null;
    }

    private FeatureHolder getFeatureHolder() {
        if (featureHolder == null) {
            featureHolder = createFeatureHolder();

            if (!hasListeners()) {
                changeForwarder = new Forwarder();
                featureHolder.addChangeListener(changeForwarder, ChangeType.UNKNOWN);
            }
        }
        return featureHolder;
    }

    public Iterator features() {
        return getFeatureHolder().features();
    }

    public int countFeatures() {
        return getFeatureHolder().countFeatures();
    }

    public FeatureHolder filter(FeatureFilter ff) {
        // FIXME: we need to optimize here
        return getFeatureHolder().filter(ff);
    }
    
    public FeatureHolder filter(FeatureFilter ff, boolean recurse) {
        if(!recurse) {
          if(FilterUtils.areDisjoint(ff, getSchema())) {
            return FeatureHolder.EMPTY_FEATURE_HOLDER;
          } else if (FilterUtils.areProperSubset(ff, getSchema())) {
            return this;
          }
        }

        // FIXME: There are other optimisations we should be making
        return getFeatureHolder().filter(ff, recurse);
    }

    public Feature createFeature(Feature.Template template) 
        throws BioException, ChangeVetoException
    {
        return getFeatureHolder().createFeature(template);
    }

    public void removeFeature(Feature f) 
        throws ChangeVetoException, BioException
    {
        getFeatureHolder().removeFeature(f);
    }

    public boolean containsFeature(Feature f) {
        return getFeatureHolder().containsFeature(f);
    }

    protected ChangeSupport getChangeSupport(ChangeType ct) {
      ChangeSupport changeSupport = super.getChangeSupport(ct);

      if (featureHolder != null) {
        changeForwarder = new Forwarder();
        featureHolder.addChangeListener(changeForwarder, ChangeType.UNKNOWN);
      }
      
      return changeSupport;
    }
	
    private class Forwarder implements ChangeListener {
	public void preChange(ChangeEvent cev)
	    throws ChangeVetoException
	{
	    getChangeSupport(cev.getType()).firePreChangeEvent(cev);
	}

	public void postChange(ChangeEvent cev) {
	    getChangeSupport(cev.getType()).firePostChangeEvent(cev);
	}
    }
}
