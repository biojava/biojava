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

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;

import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeVetoException;

/**
 * FeatureHolder which exposes all the features in a set
 * of sub-FeatureHolders.  This is provided primarily as
 * a support class for ViewSequence.  It may also be useful
 * for other applications, such as simple distributed
 * annotation systems.
 *
 * @author Thomas Down
 * @author Matthew Pocock
 * @see org.biojavax.bio.seq.RichFeatureRelationshipHolder
 */

public class MergeFeatureHolder extends AbstractFeatureHolder
  implements Serializable{

    private List featureHolders;

    /**
     * Create a new, empty, MergeFeatureHolder.
     */

    public MergeFeatureHolder() {
        featureHolders = new ArrayList();
    }

    /**
     * Create a populated MFH
     */

    private MergeFeatureHolder(List m) {
        featureHolders = m;
    }

    /**
     * Add an extra FeatureHolder to the set of FeatureHolders which
     * are merged.  This method is provided for backward compatibility,
     * and is equivalent to:
     *
     * <pre>
     *     mfh.addFeatureHolder(fh, FeatureFilter.all);
     * </pre>
     *
     * <p>
     * You should always use the two-arg version in preference if you
     * are able to define the membership of a FeatureHolder.
     * </p>
     */

    public void addFeatureHolder(FeatureHolder fh)
        throws ChangeVetoException
    {
        if(!hasListeners()) {
            featureHolders.add(fh);
        } else {
            ChangeSupport changeSupport = getChangeSupport(FeatureHolder.FEATURES);
            synchronized(changeSupport) {
                ChangeEvent ce = new ChangeEvent(this, FeatureHolder.FEATURES);
                changeSupport.firePreChangeEvent(ce);
                featureHolders.add(fh);
                changeSupport.firePostChangeEvent(ce);
            }
       }
    }

    /**
     * Remove a FeatureHolder from the set of FeatureHolders which
     * are merged.
     */

     public void removeFeatureHolder(FeatureHolder fh)
     throws ChangeVetoException
     {
       if(!hasListeners()) {
         featureHolders.remove(fh);
       } else {
         ChangeSupport changeSupport = getChangeSupport(FeatureHolder.FEATURES);
         synchronized(changeSupport) {
           ChangeEvent ce = new ChangeEvent(this, FeatureHolder.FEATURES);
           changeSupport.firePreChangeEvent(ce);
           featureHolders.remove(fh);
           changeSupport.firePostChangeEvent(ce);
         }
       }
     }

    public int countFeatures() {
        int fc = 0;
        for (Iterator i = featureHolders.iterator(); i.hasNext(); ) {
            fc += ((FeatureHolder) i.next()).countFeatures();
        }
        return fc;
    }

    public boolean containsFeature(Feature f) {
        for (Iterator i = featureHolders.iterator(); i.hasNext(); ) {
            FeatureHolder subFH = (FeatureHolder) i.next();
            FeatureFilter membership = subFH.getSchema();

            if (membership.accept(f)) {
                if(subFH.containsFeature(f)) {
                    return true;
                }
            }
        }

        return false;
    }

    /**
     * Iterate over all the features in all child FeatureHolders.
     * The Iterator may throw ConcurrantModificationException if
     * there is a change in the underlying collections during
     * iteration.
     */

    public Iterator features() {
        return new MFHIterator();
    }

    /**
     * When applied to a MergeFeatureHolder, this filters each child
     * FeatureHolder independently.
     */

    public FeatureHolder filter(FeatureFilter ff, boolean recurse) {
        List results = new ArrayList();
        for (Iterator fhi = featureHolders.iterator(); fhi.hasNext(); ) {
            FeatureHolder fh = (FeatureHolder) fhi.next();
            FeatureFilter mf = fh.getSchema();
            if (recurse && !FilterUtils.areProperSubset(mf, FeatureFilter.leaf)) {
                if (FilterUtils.areDisjoint(new FeatureFilter.Or(mf, new FeatureFilter.ByAncestor(mf)),
                                            ff))
                {
                    continue;
                }
            } else {
                if (FilterUtils.areDisjoint(mf, ff)) {
                    continue;
                }
            }

            FeatureHolder filterResult = fh.filter(ff, recurse);
            results.add(filterResult);
        }

        if (results.size() == 0) {
            return FeatureHolder.EMPTY_FEATURE_HOLDER;
        } else if (results.size() == 1) {
            return (FeatureHolder) results.get(0);
        } else {
            return new MergeFeatureHolder(results);
        }
    }

    public FeatureFilter getSchema() {
        FeatureFilter[] filters = new FeatureFilter[featureHolders.size()];
        for (int i = 0; i < filters.length; ++i) {
            filters[i] = ((FeatureHolder) featureHolders.get(i)).getSchema();
        }
        return FilterUtils.or(filters);
    }

    private class MFHIterator implements Iterator {
        private Iterator fhIterator;
        private Iterator fIterator;

        public MFHIterator() {
            fhIterator = featureHolders.iterator();
            if (fhIterator.hasNext())
                fIterator = ((FeatureHolder) fhIterator.next()).features();
            else
                fIterator = Collections.EMPTY_SET.iterator();
        }

        public boolean hasNext() {
            if (fIterator.hasNext())
                return true;
            if (fhIterator.hasNext()) {
                fIterator = ((FeatureHolder) fhIterator.next()).features();
                return hasNext();
            }
            return false;
        }

        public Object next() {
            if (fIterator.hasNext())
                return fIterator.next();
            if (fhIterator.hasNext()) {
                fIterator = ((FeatureHolder) fhIterator.next()).features();
                return next();
            }
            throw new NoSuchElementException();
        }

        public void remove() {
            throw new UnsupportedOperationException();
        }
    }
}
