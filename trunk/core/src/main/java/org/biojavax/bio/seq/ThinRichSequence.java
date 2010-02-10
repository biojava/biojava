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

package org.biojavax.bio.seq;

import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.FeatureFilter;
import org.biojava.bio.seq.FeatureHolder;
import org.biojava.bio.seq.FilterUtils;
import org.biojava.bio.seq.ProteinTools;
import org.biojava.bio.seq.SimpleFeatureHolder;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.Edit;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.ontology.InvalidTermException;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeVetoException;
import org.biojavax.Namespace;
import org.biojavax.RichObjectFactory;
import org.biojavax.bio.SimpleBioEntry;

/**
 * A simple implementation of RichSequence. It has no sequence data, and
 * delegates to a RichSequenceHandler to do sequence handling.
 * 
 * @author Richard Holland
 * @since 1.5
 */
public class ThinRichSequence extends SimpleBioEntry implements RichSequence {

	private final static String ISCIRCULAR = "X";

	private Set<Feature> features = new TreeSet<Feature>();
	private Double symListVersion;
	private boolean circular;

	/**
	 * Creates a new instance of ThinRichSequence. Note the use of Double for
	 * seqversion, which indicates that it is nullable.
	 * 
	 * @param ns
	 *            the namespace for this sequence.
	 * @param name
	 *            the name of the sequence.
	 * @param accession
	 *            the accession of the sequence.
	 * @param version
	 *            the version of the sequence.
	 * @param seqversion
	 *            the version of the symbols for the sequence.
	 */
	public ThinRichSequence(Namespace ns, String name, String accession,
			int version, Alphabet alpha, Double seqversion) {
		super(ns, name, accession, version);
		this.symListVersion = seqversion;
		this.circular = false;
		this.alphabet = alpha;
	}

	// Hibernate requirement - not for public use.
	protected ThinRichSequence() {
	}

	/**
	 * {@inheritDoc}
	 */
	public Double getSeqVersion() {
		return this.symListVersion;
	}

	/**
	 * {@inheritDoc}
	 */
	public void setSeqVersion(Double seqVersion) throws ChangeVetoException {
		if (!this.hasListeners(RichSequence.SYMLISTVERSION)) {
			this.symListVersion = seqVersion;
		} else {
			ChangeEvent ce = new ChangeEvent(this, RichSequence.SYMLISTVERSION,
					seqVersion, this.symListVersion);
			ChangeSupport cs = this
					.getChangeSupport(RichSequence.SYMLISTVERSION);
			synchronized (cs) {
				cs.firePreChangeEvent(ce);
				this.symListVersion = seqVersion;
				cs.firePostChangeEvent(ce);
			}
		}
	}

	/**
	 * {@inheritDoc}
	 */
	public void setCircular(boolean circular) throws ChangeVetoException {
		if (!this.hasListeners(RichSequence.CIRCULAR)) {
			this.circular = circular;
		} else {
			ChangeEvent ce = new ChangeEvent(this, RichSequence.CIRCULAR,
					new Boolean(circular), new Boolean(this.circular));
			ChangeSupport cs = this.getChangeSupport(RichSequence.CIRCULAR);
			synchronized (cs) {
				cs.firePreChangeEvent(ce);
				this.circular = circular;
				cs.firePostChangeEvent(ce);
			}
		}
	}

	/**
	 * {@inheritDoc}
	 */
	public boolean getCircular() {
		return this.circular;
	}

	// Hibernate requirement - not for public use.
	String getCircularChar() {
		return getCircular() ? ISCIRCULAR : null;
	}

	// Hibernate requirement - not for public use.
	void setCircularChar(final String isHiddenChar) throws ChangeVetoException {
		setCircular(isHiddenChar != null
				|| (isHiddenChar != null && isHiddenChar.length() > 0));// any
																		// character
																		// will
																		// set
	}

	/**
	 * {@inheritDoc}
	 */
	public void edit(Edit edit) throws IndexOutOfBoundsException,
			IllegalAlphabetException, ChangeVetoException {
		RichObjectFactory.getDefaultRichSequenceHandler().edit(this, edit);
	}

	/**
	 * {@inheritDoc}
	 */
	public Symbol symbolAt(int index) throws IndexOutOfBoundsException {
		return RichObjectFactory.getDefaultRichSequenceHandler().symbolAt(this,
				index);
	}

	/**
	 * {@inheritDoc}
	 */
	public List toList() {
		return RichObjectFactory.getDefaultRichSequenceHandler().toList(this);
	}

	/**
	 * {@inheritDoc}
	 */
	public String subStr(int start, int end) throws IndexOutOfBoundsException {
		return RichObjectFactory.getDefaultRichSequenceHandler().subStr(this,
				start, end);
	}

	/**
	 * {@inheritDoc}
	 */
	public SymbolList subList(int start, int end)
			throws IndexOutOfBoundsException {
		return RichObjectFactory.getDefaultRichSequenceHandler().subList(this,
				start, end);
	}

	/**
	 * {@inheritDoc}
	 */
	public String seqString() {
		return RichObjectFactory.getDefaultRichSequenceHandler()
				.seqString(this);
	}

	/**
	 * {@inheritDoc}
	 */
	public int length() {
		return this.length;
	}

	/**
	 * {@inheritDoc}
	 */
	public Iterator iterator() {
		return RichObjectFactory.getDefaultRichSequenceHandler().iterator(this);
	}

	/**
	 * {@inheritDoc}
	 */
	public Alphabet getAlphabet() {
		return this.alphabet;
	}

	// Hibernate requirement - not for public use.
	private Alphabet alphabet;

	// Hibernate requirement - not for public use.
	protected void setAlphabetName(String alphaname)
			throws IllegalSymbolException, BioException {
		if (alphaname.equals("protein"))
			alphaname = ProteinTools.getTAlphabet().getName();
		this.alphabet = AlphabetManager.alphabetForName(alphaname);
	}

	// Hibernate requirement - not for public use.
	protected String getAlphabetName() {
		if (this.alphabet == null)
			return null;
		String name = this.alphabet.getName();
		if (name.equals(ProteinTools.getTAlphabet().getName()))
			return "protein";
		else
			return name;
	}

	// Hibernate requirement - not for public use.
	private int length = 0;

	// Hibernate requirement - not for public use.
	protected void setSequenceLength(int length) {
		this.length = length;
	}

	// Hibernate requirement - not for public use.
	protected int getSequenceLength() {
		return this.length;
	}

	/**
	 * {@inheritDoc}
	 */
	public String getURN() {
		return this.getName();
	}

	/**
	 * {@inheritDoc}
	 */
	public FeatureHolder filter(FeatureFilter fc, boolean recurse) {
		SimpleFeatureHolder fh = new SimpleFeatureHolder();
		for (Iterator<Feature> i = this.features.iterator(); i.hasNext();) {
			Feature f = (RichFeature) i.next();
			try {
				if (fc.accept(f))
					fh.addFeature(f);
			} catch (ChangeVetoException e) {
				throw new RuntimeException(
						"What? You don't like our features??");
			}
		}
		return fh;
	}

	/**
	 * {@inheritDoc}
	 */
	public Feature createFeature(Feature.Template ft) throws BioException,
			ChangeVetoException {
		Feature f;
		try {
			f = new SimpleRichFeature(this, ft);
		} catch (InvalidTermException e) {
			throw new ChangeVetoException("They don't like our term", e);
		}
		if (!this.hasListeners(RichSequence.FEATURES)) {
			this.features.add(f);
		} else {
			ChangeEvent ce = new ChangeEvent(this, RichSequence.FEATURES, f,
					null);
			ChangeSupport cs = this.getChangeSupport(RichSequence.FEATURES);
			synchronized (cs) {
				cs.firePreChangeEvent(ce);
				this.features.add(f);
				cs.firePostChangeEvent(ce);
			}
		}
		return f;
	}

	/**
	 * {@inheritDoc}
	 */
	public void removeFeature(Feature f) throws ChangeVetoException,
			BioException {
		if (!(f instanceof RichFeature))
			f = RichFeature.Tools.enrich(f);
		if (!this.hasListeners(RichSequence.FEATURES)) {
			this.features.remove(f);
		} else {
			ChangeEvent ce = new ChangeEvent(this, RichSequence.FEATURES, null,
					f);
			ChangeSupport cs = this.getChangeSupport(RichSequence.FEATURES);
			synchronized (cs) {
				cs.firePreChangeEvent(ce);
				this.features.remove(f);
				cs.firePostChangeEvent(ce);
			}
		}
	}

	/**
	 * {@inheritDoc}
	 */
	public boolean containsFeature(Feature f) {
		try {
			if (!(f instanceof RichFeature))
				f = RichFeature.Tools.enrich(f);
		} catch (ChangeVetoException e) {
			// We just can't tell!
			return false;
		}
		return this.features.contains(f);
	}

	/**
	 * {@inheritDoc}
	 */
	public FeatureHolder filter(FeatureFilter filter) {
		boolean recurse = !FilterUtils.areProperSubset(filter,
				FeatureFilter.top_level);
		return this.filter(filter, recurse);
	}

	/**
	 * {@inheritDoc} <b>Warning</b> this method gives access to the original
	 * Collection not a copy. This is required by Hibernate. If you modify the
	 * object directly the behaviour may be unpredictable.
	 */
	public Set<Feature> getFeatureSet() {
		return this.features;
	} // must be original for Hibernate

	/**
	 * {@inheritDoc} <b>Warning</b> this method gives access to the original
	 * Collection not a copy. This is required by Hibernate. If you modify the
	 * object directly the behaviour may be unpredictable.
	 */
	public void setFeatureSet(Set<Feature> features) throws ChangeVetoException {
		this.features = features;
	} // must be original for Hibernate

	/**
	 * {@inheritDoc}
	 */
	public FeatureFilter getSchema() {
		return FeatureFilter.top_level;
	}

	/**
	 * {@inheritDoc} <b>Warning</b> this method gives access to the original
	 * Collection not a copy. This is required by Hibernate. If you modify the
	 * object directly the behaviour may be unpredictable.
	 */
	public Iterator<Feature> features() {
		return this.getFeatureSet().iterator();
	}

	/**
	 * {@inheritDoc}
	 */
	public int countFeatures() {
		return this.features.size();
	}

	/**
	 * {@inheritDoc}
	 */
	public SymbolList getInternalSymbolList() {
		return SymbolList.EMPTY_LIST;
	}
}
