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

import java.io.Serializable;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import org.biojava.bio.Annotation;
import org.biojava.bio.BioException;
import org.biojava.bio.SimpleAnnotation;
import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.FeatureFilter;
import org.biojava.bio.seq.FeatureHolder;
import org.biojava.bio.seq.FeatureRealizer;
import org.biojava.bio.seq.RealizingFeatureHolder;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SimpleFeatureHolder;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.Edit;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.utils.AbstractChangeable;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeListener;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;

/**
 * A basic implementation of the <code>Sequence</code> interface.
 * <p>
 * This class now implements all methods in the SymbolList interface by
 * delegating to another SymbolList object. This avoids unnecessary copying, but
 * means that any changes in the underlying SymbolList will be silently
 * reflected in the SimpleSequence. In general, SimpleSequences should
 * <em>only</em> be constructed from SymbolLists which are known to be
 * immutable.
 * </p>
 * 
 * <p>
 * By default, features attached to a SimpleSequence are realized as simple
 * in-memory implementations using <code>SimpleFeatureRealizer.DEFAULT</code>.
 * If you need alternative feature realization behaviour, any
 * <code>FeatureRealizer</code> implementation may be supplied at
 * construction-time.
 * </p>
 * More functionality and better persistence to biosql is offered by
 * {@link org.biojavax.bio.seq.SimpleRichSequence SimpleRichSequence}
 * 
 * @author Matthew Pocock
 * @author Thomas Down
 * @author Mark Schreiber
 * @serial WARNING serialized versions of this class may not be compatable with
 *         future versions of Biojava
 */
public class SimpleSequence extends AbstractChangeable implements Sequence,
		RealizingFeatureHolder, Serializable {
	//
	// This section is for the SymbolList implementation-by-delegation
	//

	private class FeatureForwarder implements ChangeListener {
		public void postChange(ChangeEvent cev) {
			ChangeSupport cs = getChangeSupport(cev.getType());
			synchronized (cs) {
				cs.firePostChangeEvent(cev);
			}
		}

		public void preChange(ChangeEvent cev) throws ChangeVetoException {
			ChangeSupport cs = getChangeSupport(cev.getType());
			synchronized (cs) {
				cs.firePreChangeEvent(cev);
			}
		}
	}
	/**
	 * 
	 */
	private static final long serialVersionUID = -8681680737943980721L;

	private Annotation annotation;

	private transient ChangeListener featureForwarder;

	private SimpleFeatureHolder featureHolder;

	private transient FeatureRealizer featureRealizer;

	private String name;

	/**
	 * Delegate SymbolList.
	 */

	private SymbolList symList;

	private String urn;

	/**
	 * Create a SimpleSequence with the symbols and alphabet of sym, and the
	 * sequence properties listed.
	 * 
	 * @param sym
	 *            the SymbolList to wrap as a sequence
	 * @param urn
	 *            the URN
	 * @param name
	 *            the name - should be unique if practical
	 * @param annotation
	 *            the annotation object to use or null
	 */
	public SimpleSequence(SymbolList sym, String urn, String name,
			Annotation annotation) {
		symList = sym;

		setURN(urn);
		setName(name);
		this.annotation = annotation;
		this.featureRealizer = FeatureImpl.DEFAULT;
	}

	//
	// Extra stuff which is unique to Sequences
	//

	/**
	 * Create a SimpleSequence using a specified FeatureRealizer.
	 * 
	 * @param sym
	 *            the SymbolList to wrap as a sequence
	 * @param urn
	 *            the URN
	 * @param name
	 *            the name - should be unique if practical
	 * @param annotation
	 *            the annotation object to use or null
	 * @param realizer
	 *            the FeatureRealizer implemetation to use when adding features
	 */
	public SimpleSequence(SymbolList sym, String urn, String name,
			Annotation annotation, FeatureRealizer realizer) {
		symList = sym;

		setURN(urn);
		setName(name);
		this.annotation = annotation;
		this.featureRealizer = realizer;
	}
	public boolean containsFeature(Feature f) {
		if (featureHolderAllocated()) {
			return getFeatureHolder().containsFeature(f);
		} else {
			return false;
		}
	}
	public int countFeatures() {
		if (featureHolderAllocated())
			return getFeatureHolder().countFeatures();
		return 0;
	}
	public Feature createFeature(Feature.Template template)
			throws BioException, ChangeVetoException {
		Feature f = realizeFeature(this, template);
		SimpleFeatureHolder fh = this.getFeatureHolder();
		synchronized (fh) {
			fh.addFeature(f);
		}
		return f;
	}
	/**
	 * Create a new feature in any FeatureHolder associated with this sequence.
	 * 
	 * @deprecated Please use new 1-arg createFeature instead.
	 */

	public Feature createFeature(FeatureHolder fh, Feature.Template template)
			throws BioException, ChangeVetoException {
		return fh.createFeature(template);
	}

	// private void readObject(ObjectInputStream s)throws IOException,
	// ClassNotFoundException{
	// s.defaultReadObject();
	// this.featureRealizer = FeatureImpl.DEFAULT;
	// }

	public void edit(Edit edit) throws ChangeVetoException {
		throw new ChangeVetoException("Can't edit the underlying SymbolList");
	}

	protected boolean featureHolderAllocated() {
		return featureHolder != null;
	}

	public Iterator<Feature> features() {
		if (featureHolderAllocated())
			return getFeatureHolder().features();
		return Collections.EMPTY_LIST.iterator();
	}

	public FeatureHolder filter(FeatureFilter filter) {
		return getFeatureHolder().filter(filter);
	}

	public FeatureHolder filter(FeatureFilter ff, boolean recurse) {
		if (featureHolderAllocated()) {
			return getFeatureHolder().filter(ff, recurse);
		}
		return FeatureHolder.EMPTY_FEATURE_HOLDER;
	}

	public Alphabet getAlphabet() {
		return symList.getAlphabet();
	}

	public Annotation getAnnotation() {
		if (annotation == null)
			annotation = new SimpleAnnotation();
		return annotation;
	}

	protected ChangeSupport getChangeSupport(ChangeType ct) {
		ChangeSupport changeSupport = super.getChangeSupport(ct);

		if (featureForwarder == null && featureHolder != null) {
			featureForwarder = new FeatureForwarder();
			featureHolder.addChangeListener(featureForwarder,
					ChangeType.UNKNOWN);
		}

		return changeSupport;
	}

	protected SimpleFeatureHolder getFeatureHolder() {
		if (featureHolder == null) {
			featureHolder = new SimpleFeatureHolder(FeatureFilter.top_level);
		}
		return featureHolder;
	}

	public String getName() {
		return name;
	}

	public FeatureFilter getSchema() {
		return getFeatureHolder().getSchema();
	}

	public String getURN() {
		return urn;
	}

	public Iterator iterator() {
		return symList.iterator();
	}

	public int length() {
		return symList.length();
	}

	public Feature realizeFeature(FeatureHolder parent,
			Feature.Template template) throws BioException {
		return featureRealizer.realizeFeature(this, parent, template);
	}

	/**
	 * Remove a feature attached to this sequence.
	 */

	public void removeFeature(Feature f) throws ChangeVetoException,
			BioException {
		getFeatureHolder().removeFeature(f);
	}

	public String seqString() {
		return symList.seqString();
	}

	/**
	 *Assign a name to this sequence
	 */

	public void setName(String name) {
		this.name = name;
	}

	/**
	 *Provide the URN for this sequence
	 */

	public void setURN(String urn) {
		this.urn = urn;
	}

	public SymbolList subList(int start, int end) {
		return symList.subList(start, end);
	}

	public String subStr(int start, int end) {
		return symList.subStr(start, end);
	}

	//
	// Changeable stuff
	//

	public Symbol symbolAt(int index) {
		return symList.symbolAt(index);
	}

	public List toList() {
		return symList.toList();
	}

	public String toString() {
		return super.toString() + " name: " + getName();
	}
}
