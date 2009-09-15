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

import java.util.Collections;
import java.util.Iterator;

import org.biojava.bio.Annotatable;
import org.biojava.bio.Annotation;
import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.bio.SimpleAnnotation;
import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.FeatureFilter;
import org.biojava.bio.seq.FeatureHolder;
import org.biojava.bio.seq.FilterUtils;
import org.biojava.bio.seq.RealizingFeatureHolder;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SimpleFeatureHolder;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.ontology.OntoTools;
import org.biojava.ontology.Term;
import org.biojava.utils.AbstractChangeable;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeForwarder;
import org.biojava.utils.ChangeListener;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;

/**
 * A no-frills implementation of a feature.
 * 
 * @author Matthew Pocock
 * @author Thomas Down
 * @author Kalle N�slund
 * @author Paul Seed
 * @author Len Trigg
 * @see org.biojavax.bio.seq.SimpleRichFeature
 */

public class SimpleFeature
extends
AbstractChangeable
implements
Feature,
RealizingFeatureHolder,
java.io.Serializable
{
	private transient ChangeListener annotationForwarder;
	private transient ChangeListener featureForwarder;

	/**
	 * The FeatureHolder that we will delegate the FeatureHolder interface too.
	 * This is lazily instantiated.
	 */
	private SimpleFeatureHolder featureHolder;

	/**
	 * The location of this feature.
	 */
	private Location loc;
	/**
	 * The type of this feature - something like Exon.
	 * This is included for cheap interoperability with GFF.
	 */
	private String type;
	/**
	 * The source of this feature - the program that generated it.
	 * This is included for cheap interoperability with GFF.
	 */
	private String source;
	/**
	 * Our parent FeatureHolder.
	 */
	private FeatureHolder parent;
	/**
	 * The annotation object.
	 * This is lazily instantiated.
	 */
	private Annotation annotation;

	private Term typeTerm;
	private Term sourceTerm;

	/**
	 * A utility function to retrieve the feature holder delegate, creating it if
	 * necessary.
	 *
	 * @return  the FeatureHolder delegate
	 */
	protected SimpleFeatureHolder getFeatureHolder() {
		if(featureHolder == null) {
			featureHolder = new SimpleFeatureHolder();
		}
		return featureHolder;
	}

	/**
	 * A utility function to find out if the feature holder delegate has been
	 * instantiated yet. If it has not, we may avoid instantiating it by returning
	 * some pre-canned result.
	 *
	 * @return true if the feature holder delegate has been created and false
	 *         otherwise
	 */
	protected boolean featureHolderAllocated() {
		return featureHolder != null;
	}

	protected ChangeSupport getChangeSupport(ChangeType ct) {
		ChangeSupport cs = super.getChangeSupport(ct);

		if(
				(annotationForwarder == null) &&
				(ct.isMatchingType(Annotatable.ANNOTATION) || Annotatable.ANNOTATION.isMatchingType(ct))
		) {
			annotationForwarder =
				new ChangeForwarder.Retyper(this, cs, Annotation.PROPERTY);
			getAnnotation().addChangeListener(
					annotationForwarder,
					Annotatable.ANNOTATION
			);
		}

		if(
				(featureForwarder == null) &&
				(ct == null || ct == FeatureHolder.FEATURES)
		) {
			featureForwarder = new ChangeForwarder(
					this,
					cs
			);
			getFeatureHolder().addChangeListener(
					featureForwarder,
					FeatureHolder.FEATURES
			);
		}

		return cs;
	}

	public Location getLocation() {
		return loc;
	}

	public void setLocation(Location loc)
	throws ChangeVetoException {
		if(hasListeners()) {
			ChangeSupport cs = getChangeSupport(LOCATION);
			synchronized(cs) {
				ChangeEvent ce = new ChangeEvent(this, LOCATION, loc, this.loc);
				synchronized(ce) {
					cs.firePreChangeEvent(ce);
				}
				this.loc = loc;
				synchronized(ce) {
					cs.firePostChangeEvent(ce);
				}
			}
		} else {
			this.loc = loc;
		}
	}

	public Term getTypeTerm() {
		return typeTerm;
	}

	public String getType() {
		if(type != null) {
			return type;
		} else if (typeTerm != null) {
			return typeTerm.getName();
		} else {
			return "";
		}
	}

	public void setType(String type)
	throws ChangeVetoException {
		if(hasListeners()) {
			ChangeSupport cs = getChangeSupport(TYPE);
			synchronized(cs) {
				ChangeEvent ce = new ChangeEvent(this, TYPE, type, this.type);
				synchronized(ce) {
					cs.firePreChangeEvent(ce);
				}
				this.type = type;
				synchronized(ce) {
					cs.firePostChangeEvent(ce);
				}
			}
		} else {
			this.type = type;
		}
	}

	public void setTypeTerm(Term t)
	throws ChangeVetoException
	{
		if(hasListeners()) {
			ChangeSupport cs = getChangeSupport(TYPE);
			synchronized (cs) {
				ChangeEvent ce_term = new ChangeEvent(this, TYPETERM, t, this.getTypeTerm());
				ChangeEvent ce_name = new ChangeEvent(this, TYPE, t.getName(), this.getType());
				cs.firePreChangeEvent(ce_term);
				cs.firePreChangeEvent(ce_name);
				this.typeTerm = t;
				cs.firePostChangeEvent(ce_term);
				cs.firePostChangeEvent(ce_name);
			}
		} else {
			this.typeTerm = t;
		}
	}

	public String getSource() {
		if(source != null) {
			return source;
		} else if (sourceTerm != null) {
			return sourceTerm.getName();
		} else {
			return "";
		}
	}

	public Term getSourceTerm() {
		return sourceTerm;
	}

	public FeatureHolder getParent() {
		return parent;
	}

	public void setSource(String source)
	throws ChangeVetoException {
		if(hasListeners()) {
			ChangeSupport cs = getChangeSupport(SOURCE);
			synchronized(cs) {
				ChangeEvent ce = new ChangeEvent(this, SOURCE, this.source, source);
				cs.firePreChangeEvent(ce);
				this.source = source;
				cs.firePostChangeEvent(ce);
			}
		} else {
			this.source = source;
		}
	}


	public void setSourceTerm(Term t)
	throws ChangeVetoException
	{
		if(hasListeners()) {
			ChangeSupport cs = getChangeSupport(TYPE);
			synchronized (cs) {
				ChangeEvent ce_term = new ChangeEvent(this, SOURCETERM, t, this.getSourceTerm());
				ChangeEvent ce_name = new ChangeEvent(this, SOURCE, t.getName(), this.getSource());
				cs.firePreChangeEvent(ce_term);
				cs.firePreChangeEvent(ce_name);
				this.sourceTerm = t;
				cs.firePostChangeEvent(ce_term);
				cs.firePostChangeEvent(ce_name);
			}
		} else {
			this.sourceTerm = t;
		}
	}

	public Sequence getSequence() {
		FeatureHolder fh = this;
		while (fh instanceof Feature) {
			fh = ((Feature) fh).getParent();
		}
		try {
			return (Sequence) fh;
		} catch (ClassCastException ex) {
			throw new BioError("Feature doesn't seem to have a Sequence ancestor: " + fh);
		}
	}

	public Annotation getAnnotation() {
		if(annotation == null)
			annotation = new SimpleAnnotation();
		return annotation;
	}

	public SymbolList getSymbols() {
		return getLocation().symbols(getSequence());
	}

	public int countFeatures() {
		if(featureHolderAllocated())
			return getFeatureHolder().countFeatures();
		return 0;
	}

	public Iterator features() {
		if(featureHolderAllocated())
			return getFeatureHolder().features();
		return Collections.EMPTY_LIST.iterator();
	}

	public void removeFeature(Feature f)
	throws ChangeVetoException {
		getFeatureHolder().removeFeature(f);
	}

	public boolean containsFeature(Feature f) {
		if(featureHolderAllocated()) {
			return getFeatureHolder().containsFeature(f);
		} else {
			return false;
		}
	}


	public FeatureHolder filter(FeatureFilter ff) {
		FeatureFilter childFilter = new FeatureFilter.Not(FeatureFilter.top_level);

		if (FilterUtils.areDisjoint(ff, childFilter)) {
			return FeatureHolder.EMPTY_FEATURE_HOLDER;
		} else if (featureHolderAllocated()) {
			return getFeatureHolder().filter(ff);
		} else {
			return FeatureHolder.EMPTY_FEATURE_HOLDER;
		}
	}

	public FeatureHolder filter(FeatureFilter ff, boolean recurse) {
		if(featureHolderAllocated())
			return getFeatureHolder().filter(ff, recurse);
		return FeatureHolder.EMPTY_FEATURE_HOLDER;
	}

	public Feature.Template makeTemplate() {
		Feature.Template ft = new Feature.Template();
		fillTemplate(ft);
		return ft;
	}

	protected void fillTemplate(Feature.Template ft) {
		ft.location = getLocation();
		ft.type = getType();
		ft.source = getSource();
		ft.annotation = getAnnotation();
		ft.sourceTerm = getSourceTerm();
		ft.typeTerm = getTypeTerm();
	}

	/**
	 * Create a <code>SimpleFeature</code> on the given sequence.
	 * The feature is created underneath the parent <code>FeatureHolder</code>
	 * and populated directly from the template fields. However,
	 * if the template annotation is the <code>Annotation.EMPTY_ANNOTATION</code>,
	 * an empty <code>SimpleAnnotation</code> is attached to the feature instead.
	 * @param sourceSeq the source sequence
	 * @param parent the parent sequence or feature
	 * @param template the template for the feature
	 */
	public SimpleFeature(Sequence sourceSeq,
			FeatureHolder parent,
			Feature.Template template) {
		if (template.location == null) {
			throw new IllegalArgumentException(
					"Location can not be null. Did you mean Location.EMPTY_LOCATION? " +
					template.toString()
			);
		}
		if(!(parent instanceof Feature) && !(parent instanceof Sequence)) {
			throw new IllegalArgumentException("Parent must be sequence or feature, not: " + parent.getClass() + " " + parent);
		}

		if (template.location.getMin() < 1 || template.location.getMax() > sourceSeq.length()) {
			//throw new IllegalArgumentException("Location " + template.location.toString() + " is outside 1.." + sourceSeq.length());
		}

		this.parent = parent;
		this.loc = template.location;
		this.typeTerm = template.typeTerm != null ? template.typeTerm : OntoTools.ANY;
		this.sourceTerm = template.sourceTerm != null ? template.sourceTerm : OntoTools.ANY;
		this.type = template.type != null ? template.type : typeTerm.getName();
		this.source = template.source != null ? template.source : sourceTerm.getName();
		if (this.type == null) {
			throw new NullPointerException("Either type or typeTerm must have a non-null value");
		}
		if (this.source == null) {
			throw new NullPointerException("Either source or sourceTerm must have a non-null value");
		}
		this.annotation = template.annotation != null ? new SimpleAnnotation(template.annotation) : null;
	}

	public String toString() {
		return "Feature " + getType() + " " +
		getSource() + " " + getLocation();
	}

	public Feature realizeFeature(FeatureHolder fh, Feature.Template templ)
	throws BioException
	{
		try {
			RealizingFeatureHolder rfh = (RealizingFeatureHolder) getParent();
			return rfh.realizeFeature(fh, templ);
		} catch (ClassCastException ex) {
			throw new BioException("Couldn't propagate feature creation request.");
		}
	}

	public Feature createFeature(Feature.Template temp)
	throws BioException, ChangeVetoException
	{
		Feature f = realizeFeature(this, temp);
		getFeatureHolder().addFeature(f);
		return f;
	}

	public int hashCode() {
		return makeTemplate().hashCode();
	}

	public boolean equals(Object o) {
		if (! (o instanceof Feature)) {
			return false;
		}

		Feature fo = (Feature) o;
		if (! fo.getSequence().equals(getSequence()))
			return false;

		return makeTemplate().equals(fo.makeTemplate());
	}

	public FeatureFilter getSchema() {
		return new FeatureFilter.ByParent(new FeatureFilter.ByFeature(this));
	}
}
