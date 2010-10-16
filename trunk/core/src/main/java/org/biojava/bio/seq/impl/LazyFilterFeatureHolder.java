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

import java.util.Iterator;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.FeatureFilter;
import org.biojava.bio.seq.FeatureHolder;
import org.biojava.bio.seq.FilterUtils;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeListener;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;

/**
 * FeatureHolder which lazily applies a specified filter to another FeatureHolder.
 * This means that when you use a second filter to query the LazyFilterFeatureHolder,
 * the underlying holder receives a single <code>filter</code> request with the
 * two queries ANDed together.
 *
 * @author Thomas Down
 * @since 1.3
 */

public class LazyFilterFeatureHolder implements FeatureHolder {
	private FeatureHolder featureHolder;
	private FeatureFilter featureFilter;
	private transient ChangeSupport changeSupport;

	public LazyFilterFeatureHolder(FeatureHolder fh,
			FeatureFilter ff)
	{
		this.featureHolder = fh;
		this.featureFilter = ff;
	}

	public Iterator features() {
		return featureHolder.filter(featureFilter, false).features();
	}

	public int countFeatures() {
		return featureHolder.filter(featureFilter, false).countFeatures();
	}

	public boolean containsFeature(Feature f) {
		if (featureFilter.accept(f)) {
			return featureHolder.containsFeature(f);
		} else {
			return false;
		}
	}

	public FeatureHolder filter(FeatureFilter ff) {
		if (FilterUtils.areDisjoint(ff, featureFilter)) {
			return FeatureHolder.EMPTY_FEATURE_HOLDER;
		}
		return featureHolder.filter(new FeatureFilter.And(ff, featureFilter));
	}

	public FeatureHolder filter(FeatureFilter ff, boolean recurse) {
		if (FilterUtils.areDisjoint(ff, featureFilter)) {
			return FeatureHolder.EMPTY_FEATURE_HOLDER;
		}
		return featureHolder.filter(new FeatureFilter.And(ff, featureFilter), recurse);
	}


	public Feature createFeature(Feature.Template temp)
	throws ChangeVetoException, BioException
	{

		Feature f= null;

		synchronized (featureHolder) {
			f = featureHolder.createFeature(temp);	
		}

		return f;
	}

	public void removeFeature(Feature f)
	throws ChangeVetoException, BioException
	{
		synchronized (featureHolder) {
			featureHolder.removeFeature(f);
		}

	}


	protected boolean hasListeners() {
		return changeSupport != null;
	}

	protected ChangeSupport getChangeSupport() {
		if(changeSupport != null) {
			return changeSupport;
		}

		synchronized(this) {
			if(changeSupport == null) {
				changeSupport = new ChangeSupport();
			}
			synchronized (featureHolder) {
				featureHolder.addChangeListener(new LFFHChangeForwarder(), ChangeType.UNKNOWN);		
			}

		}

		return changeSupport;
	}

	private class LFFHChangeForwarder implements ChangeListener {
		public void preChange(ChangeEvent cev)
		throws ChangeVetoException
		{
			ChangeEvent fcev = getForwardedEvent(cev);
			if (fcev != null) {
				ChangeSupport cs = getChangeSupport();
				synchronized(cs){
					cs.firePreChangeEvent(fcev);
				}
			}
		}

		public void postChange(ChangeEvent cev)
		{
			ChangeEvent fcev = getForwardedEvent(cev);
			if (fcev != null) {
				ChangeSupport cs = getChangeSupport();
				synchronized(cs){
					cs.firePostChangeEvent(fcev);
				}
			}
		}
	}

	/**
	 * Only forward events concerning features which are either accepted by
	 * our filter, or are parents of a feature which is.
	 *
	 * In future this maybe ought to look further down the event chain
	 * to reject forwarded events from uninterested features.
	 */

	private ChangeEvent getForwardedEvent(ChangeEvent cev) {
		Object change = cev.getChange();
		if (! (change instanceof Feature)) {
			change = null;
		}
		Object previous = cev.getPrevious();
		if (! (previous instanceof Feature)) {
			previous = null;
		}
		boolean forward = false;
		if (change == null && previous == null) {
			forward = true;
		} else {
			forward = isInterestingFeature((Feature) previous) || 
			isInterestingFeature((Feature) change);
		}
		if (forward) {
			return new ChangeEvent(this,
					cev.getType(),
					cev.getChange(),
					cev.getPrevious(),
					cev);
		} else {
			return null;
		}
	}

	private boolean isInterestingFeature(Feature f) {
		if (f == null) {
			return false;
		} else if (featureFilter.accept(f)) {
			return true;
		} else {
			return f.filter(featureFilter).countFeatures() > 0;
		}
	}

	public final void addChangeListener(ChangeListener cl) {
		addChangeListener(cl, ChangeType.UNKNOWN);
	}

	public final void addChangeListener(ChangeListener cl, ChangeType ct) {
		if (!isUnchanging(ct)) {
			ChangeSupport cs = getChangeSupport();
			cs.addChangeListener(cl, ct);
		}
	}

	public final void removeChangeListener(ChangeListener cl) {
		removeChangeListener(cl, ChangeType.UNKNOWN);
	}

	public final void removeChangeListener(ChangeListener cl, ChangeType ct) {
		if(hasListeners()) {
			ChangeSupport cs = getChangeSupport();
			cs.removeChangeListener(cl, ct);
		}
	}

	public boolean isUnchanging(ChangeType ct) {
		return featureHolder.isUnchanging(ct);
	}

	public FeatureFilter getSchema() {
		return new FeatureFilter.And(featureFilter,
				new FeatureFilter.Or(featureHolder.getSchema(), 
						new FeatureFilter.ByAncestor(featureHolder.getSchema())));
	}
}
