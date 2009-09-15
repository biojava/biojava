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

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.biojava.bio.Annotatable;
import org.biojava.bio.Annotation;
import org.biojava.bio.BioException;
import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.FeatureFilter;
import org.biojava.bio.seq.FeatureHolder;
import org.biojava.bio.seq.FilterUtils;
import org.biojava.bio.seq.RemoteFeature;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.seq.db.IllegalIDException;
import org.biojava.bio.seq.projection.ProjectedFeatureHolder;
import org.biojava.bio.seq.projection.Projection;
import org.biojava.bio.seq.projection.ReparentContext;
import org.biojava.bio.seq.projection.TranslateFlipContext;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.Edit;
import org.biojava.bio.symbol.FuzzyLocation;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.LocationTools;
import org.biojava.bio.symbol.RangeLocation;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.ontology.InvalidTermException;
import org.biojava.ontology.Term;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeForwarder;
import org.biojava.utils.ChangeListener;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;

/**
 * View a sub-section of a given sequence object, including all the
 * features intersecting that region.
 *
 * <p>
 * All features entirely contained within the region are projected by just
 * translating their locations. The features that overlap the region are
 * replaced by RemoteFeature instances with fuzzy locations that are
 * trunchated to fit inside the sub-section. All features not contained by
 * the region are not projected.
 * </p>
 *
 * @author Thomas Down
 * @author Matthew Pocock
 * @since 1.2
 */

public class SubSequence
implements Sequence, java.io.Serializable
{
	private final Sequence parent;
	private final SymbolList symbols;
	private final String name;
	private final String uri;
	private final Annotation annotation;
	private final RangeLocation parentLocation;

	private transient ProjectedFeatureHolder features;
	private transient ChangeSupport changeSupport;
	private transient ChangeListener seqListener;
	protected transient ChangeForwarder annotationForwarder;

	private void allocChangeSupport() {
		if (seqListener == null) {
			installSeqListener();
		}
		changeSupport = new ChangeSupport();
	}

	private void installSeqListener() {
		seqListener = new ChangeListener() {
			public void preChange(ChangeEvent cev)
			throws ChangeVetoException {
				if (changeSupport != null) {
					changeSupport.firePreChangeEvent(makeChainedEvent(cev));
				}
			}

			public void postChange(ChangeEvent cev) {
				if (changeSupport != null) {
					changeSupport.firePostChangeEvent(makeChainedEvent(cev));
				}
			}

			private ChangeEvent makeChainedEvent(ChangeEvent cev) {
				return new ChangeEvent(SubSequence.this,
						FeatureHolder.FEATURES,
						null, null,
						cev);
			}
		};
		parent.addChangeListener(seqListener, FeatureHolder.FEATURES);
	}

	/**
	 * Construct a new SubSequence of the specified sequence.
	 *   Generally you would use the SequenceTools.subSequence() methods
	 *  to get an instance of this class.
	 * @param seq A sequence to view
	 * @param start The start of the range to view
	 * @param end The end of the range to view
	 * @param name Name for the subsequence
	 * @throws java.lang.IndexOutOfBoundsException is the start or end position is illegal.
	 */
	public SubSequence(Sequence seq,
			final int start,
			final int end,
			final String name)
	{
		this.parent = seq;
		this.parentLocation = new RangeLocation(start, end);
		this.symbols = seq.subList(start, end);
		this.name = name;
		this.uri = seq.getURN() + "?start=" + start + ";end=" + end;
		this.annotation = seq.getAnnotation();
	}

	/**
	 * Construct a new SubSequence of the specified sequence. Generally
	 * you would use the SequenceTools.subSequence() methods to get an
	 * instance of this class.
	 *
	 * @param seq A sequence to view
	 * @param start The start of the range to view
	 * @param end The end of the range to view
	 * @throws java.lang.IndexOutOfBoundsException if the start or end position is illegal.
	 */

	public SubSequence(Sequence seq,
			final int start,
			final int end) {
		this(seq, start, end, seq.getName() + " (" + start + " - " + end + ")");
	}

	//
	// SymbolList stuff
	//

	public Symbol symbolAt(int pos) {
		return symbols.symbolAt(pos);
	}

	public Alphabet getAlphabet() {
		return symbols.getAlphabet();
	}

	public SymbolList subList(int start, int end) {
		return symbols.subList(start, end);
	}

	public String seqString() {
		return symbols.seqString();
	}

	public String subStr(int start, int end) {
		return symbols.subStr(start, end);
	}

	public List toList() {
		return symbols.toList();
	}

	public int length() {
		return symbols.length();
	}

	public Iterator iterator() {
		return symbols.iterator();
	}

	public void edit(Edit edit)
	throws ChangeVetoException {
		throw new ChangeVetoException("Can't edit SubSequences");
	}

	//
	// Implements featureholder
	//

	public int countFeatures() {
		return getFeatures().countFeatures();
	}

	public Iterator features() {
		return getFeatures().features();
	}

	public FeatureHolder filter(FeatureFilter ff, boolean recurse) {
		return getFeatures().filter(ff, recurse);
	}

	public FeatureHolder filter(FeatureFilter ff) {
		return getFeatures().filter(ff);
	}

	public boolean containsFeature(Feature f) {
		return getFeatures().containsFeature(f);
	}

	public Feature createFeature(Feature.Template templ)
	throws BioException, ChangeVetoException {

		ProjectedFeatureHolder featureHolder = getFeatures();

		Feature f = null;
		synchronized (featureHolder) {
			f = featureHolder.createFeature(templ);
		}
		return f;

	}

	public void removeFeature(Feature f)
	throws ChangeVetoException, BioException {
		ProjectedFeatureHolder featureHolder = getFeatures();

		synchronized (featureHolder){
			featureHolder.removeFeature(f);
		}
	}

	public FeatureFilter getSchema() {
		return getFeatures().getSchema();
	}

	protected ProjectedFeatureHolder getFeatures() {
		if (features == null) {
			FeatureHolder clipped = new ProjectedFeatureHolder(
					new SubProjectedFeatureContext(parent, parentLocation) );
			features = new ProjectedFeatureHolder(
					new TranslateFlipContext(this,
							clipped,
							- parentLocation.getMin() + 1) );
		}
		return features;
	}

	//
	// Identifiable
	//

	public String getName() {
		return name;
	}

	public String getURN() {
		return uri;
	}

	//
	// Annotatable
	//

	public Annotation getAnnotation() {
		return annotation;
	}

	/**
	 * Return the parent sequence of which this is a partial view
	 *
	 * @since 1.3
	 */

	public Sequence getSequence() {
		return this.parent;
	}

	public int getStart() {
		return parentLocation.getMin();
	}

	public int getEnd() {
		return parentLocation.getMax();
	}

	public void addChangeListener(ChangeListener cl, ChangeType ct) {
		if (changeSupport == null) {
			allocChangeSupport();
		}

		if (annotationForwarder == null && ct == Annotatable.ANNOTATION) {
			annotationForwarder =
				new ChangeForwarder.Retyper(this, changeSupport, Annotation.PROPERTY);
			getAnnotation().addChangeListener(annotationForwarder,
					Annotatable.ANNOTATION);
		}

		changeSupport.addChangeListener(cl, ct);
	}

	public void addChangeListener(ChangeListener cl) {
		addChangeListener(cl, ChangeType.UNKNOWN);
	}

	public void removeChangeListener(ChangeListener cl, ChangeType ct) {
		if (changeSupport != null) {
			changeSupport.removeChangeListener(cl, ct);
		}
	}

	public void removeChangeListener(ChangeListener cl) {
		removeChangeListener(cl, ChangeType.UNKNOWN);
	}

	public boolean isUnchanging(ChangeType ct) {
		return parent.isUnchanging(ct);
	}

	/**
	 * TargetContext that implements the mapping between the parent sequence and this
	 * sub-sequence.
	 *
	 * <p>
	 * This context is public because all contexts must be public, and not because
	 * it is usefull to you or part of the API. Don't use it directly.
	 * </p>
	 *
	 * <p>
	 * The context extends TranslateFlipContext so that it can translate all
	 * features within the parent location to being from index 1 in the
	 * projection. It also transforms and reverts locations so that locations
	 * falling outside the parent location are trunchated and replaced by
	 * fuzzy locations. Lastly, features with fuzzy locations are replaced by
	 * RemoteFeature instances.
	 *</p>
	 *
	 * @author Matthew Pocock
	 * @author Thomas Down
	 * @since 1.4
	 */
	public static class SubProjectedFeatureContext
	extends ReparentContext {
		private static final FeatureFilter REMOTE_FILTER =
			new FeatureFilter.ByClass(RemoteFeature.class);

		private static RemoteFeature.Resolver resolver
		= new RemoteFeature.Resolver() {
			public Feature resolve(RemoteFeature rFeat) throws IllegalIDException, BioException {
				if(!(rFeat instanceof SubRemote)) {
					throw new BioException("Unable to resolve feature: " + rFeat);
				}

				return rFeat.getRemoteFeature();
			}
		};

		private final RangeLocation parentLocation;
		private final FeatureFilter remoteLocationFilter;
		private final FeatureFilter clippingFilter;

		private SubProjectedFeatureContext(FeatureHolder wrapped,
				RangeLocation parentLocation) {
			super(FeatureHolder.EMPTY_FEATURE_HOLDER,
					new LazyFilterFeatureHolder(wrapped,
							FilterUtils.overlapsLocation(
									parentLocation)));

			this.parentLocation = parentLocation;
			this.remoteLocationFilter = new FeatureFilter.Not(
					FilterUtils.containedByLocation(parentLocation));
			this.clippingFilter = FilterUtils.overlapsLocation(parentLocation);
		}

		public Location projectLocation(Location toTransform) {
			if(LocationTools.overlaps(parentLocation, toTransform) &&
					!LocationTools.contains(parentLocation, toTransform))
			{
				toTransform = fuzzyize(toTransform);
			}

			return toTransform;
		}

		public Feature projectFeature(Feature origFeat) {
			if(remoteLocationFilter.accept(origFeat)) {
				Location loc = fuzzyize(origFeat.getLocation());
				List regions = makeRegions(origFeat.getLocation(),
						origFeat.getSequence().getName());

				origFeat = new SubRemote(origFeat, loc, regions);
			}

			return super.projectFeature(origFeat);
		}

		public Feature revertFeature(Feature projFeat) {
			Feature origFeat = super.revertFeature(projFeat);

			if(!(origFeat instanceof Projection) &&
					(origFeat instanceof SubRemote) )
			{
				origFeat = ((SubRemote) origFeat).getRemoteFeature();
			}

			return origFeat;
		}

		public FeatureHolder projectChildFeatures(Feature f, FeatureHolder parent) {
			return new LazyFilterFeatureHolder(super.projectChildFeatures(f, parent),
					clippingFilter);
		}

		private Location fuzzyize(Location loc) {
			List locList = new ArrayList();

			for(Iterator i = loc.blockIterator(); i.hasNext(); ) {
				Location l = (Location) i.next();
				if(LocationTools.overlaps(parentLocation, l)) {
					boolean fuzzyLeft = l.getMin() < parentLocation.getMin();
					boolean fuzzyRight = l.getMax() > parentLocation.getMax();

					if(fuzzyLeft || fuzzyRight) {
						locList.add(new FuzzyLocation(Math.min(l.getMin(), parentLocation.getMin()),
								Math.max(l.getMax(), parentLocation.getMax()),
								Math.max(l.getMin(), parentLocation.getMin()),
								Math.min(l.getMax(), parentLocation.getMax()),
								fuzzyLeft,
								fuzzyRight,
								FuzzyLocation.RESOLVE_INNER
						));
					} else {
						locList.add(l);
					}
				}
			}

			return LocationTools.union(locList);
		}

		private List makeRegions(Location oldLoc, String seqID) {
			List regions = new ArrayList();

			for(Iterator i = oldLoc.blockIterator(); i.hasNext(); ) {
				Location loc = (Location) i.next();

				if(!LocationTools.overlaps(loc, parentLocation)) {
					regions.add(new RemoteFeature.Region(loc, seqID, true));
				} else if(LocationTools.contains(parentLocation, loc)) {
					regions.add(new RemoteFeature.Region(loc, null, false));
				} else {
					// straddles boundary
					Location remote = LocationTools.subtract(loc, parentLocation);
					Location local = LocationTools.subtract(parentLocation, loc);
					RemoteFeature.Region remoteR = new RemoteFeature.Region(remote, seqID, true);
					RemoteFeature.Region localR = new RemoteFeature.Region(local, null, false);

					if(remote.getMin() < local.getMin()) {
						regions.add(remoteR);
						regions.add(localR);
					} else {
						regions.add(localR);
						regions.add(remoteR);
					}
				}
			}

			return regions;
		}

		protected FilterUtils.FilterTransformer getReverter() {
			final FilterUtils.FilterTransformer delegate = super.getReverter();

			return new FilterUtils.FilterTransformer() {
				public FeatureFilter transform(FeatureFilter filt) {
					// I may have got this logic wrong
					if (filt.equals(REMOTE_FILTER)) {
						filt = remoteLocationFilter;
					}

					filt = new FeatureFilter.And(
							filt,
							FilterUtils.overlapsLocation(parentLocation)
					);

					filt = delegate.transform(filt);

					return filt;
				}
			};
		}

		protected FilterUtils.FilterTransformer getTransformer() {
			return super.getTransformer();
		}

		// fixme: this isn't wired to forward events
		// however, the context should be skipping these anyhow during
		// projectFeature and restoreFeature.
		private class SubRemote
		implements RemoteFeature {
			private final Feature wrapped;
			private final Location loc;
			private final List regionList;

			public SubRemote(Feature wrapped, Location loc, List regionList) {
				this.wrapped = wrapped;
				this.loc = loc;
				this.regionList = regionList;
			}

			public Location getLocation() {
				return loc;
			}

			public void setLocation(Location loc)
			throws ChangeVetoException {
				throw new ChangeVetoException("Can't set location: " + loc +
						" on " + this);
			}

			public List getRegions() {
				return regionList;
			}

			public RemoteFeature.Resolver getResolver() {
				return resolver;
			}

			public Feature getRemoteFeature() {
				return wrapped;
			}

			public StrandedFeature.Strand getStrand() {
				if(wrapped instanceof StrandedFeature) {
					return ((StrandedFeature) wrapped).getStrand();
				} else {
					return StrandedFeature.UNKNOWN;
				}
			}

			public void setStrand(StrandedFeature.Strand strand)
			throws ChangeVetoException
			{
				if(wrapped instanceof StrandedFeature) {
					((StrandedFeature) wrapped).setStrand(strand);
				} else {
					throw new ChangeVetoException("Can't set strand. The underlying feature is not stranded: " + wrapped);
				}
			}

			public String getType() {
				return wrapped.getType();
			}

			public void setType(String type)
			throws ChangeVetoException {
				wrapped.setType(type);
			}

			public Term getTypeTerm() {
				return wrapped.getTypeTerm();
			}

			public void setTypeTerm(Term t) throws ChangeVetoException, InvalidTermException {
				wrapped.setTypeTerm(t);
			}

			public Term getSourceTerm() {
				return wrapped.getSourceTerm();
			}

			public void setSourceTerm(Term t) throws ChangeVetoException, InvalidTermException {
				wrapped.setSourceTerm(t);
			}

			public String getSource() {
				return wrapped.getSource();
			}

			public void setSource(String source)
			throws ChangeVetoException {
				wrapped.setSource(source);
			}

			public SymbolList getSymbols() {
				return wrapped.getSymbols();
			}

			public FeatureHolder getParent() {
				return wrapped.getParent();
			}

			public Sequence getSequence() {
				return wrapped.getSequence();
			}

			public Feature.Template makeTemplate() {
				return wrapped.makeTemplate();
			}

			public int countFeatures() {
				return wrapped.countFeatures();
			}

			public Iterator features() {
				return wrapped.features();
			}

			public FeatureHolder filter(FeatureFilter fc, boolean recurse) {
				return wrapped.filter(fc, recurse);
			}

			public FeatureHolder filter(FeatureFilter filter) {
				return wrapped.filter(filter);
			}

			public Feature createFeature(Feature.Template ft)
			throws BioException, ChangeVetoException {
				return wrapped.createFeature(ft);
			}

			public void removeFeature(Feature f)
			throws ChangeVetoException, BioException {
				wrapped.removeFeature(f);
			}

			public boolean containsFeature(Feature f) {
				return wrapped.containsFeature(f);
			}

			public FeatureFilter getSchema() {
				return wrapped.getSchema();
			}

			public void addChangeListener(ChangeListener cl) {
				wrapped.addChangeListener(cl);
			}

			public void addChangeListener(ChangeListener cl, ChangeType ct) {
				wrapped.addChangeListener(cl, ct);
			}

			public void removeChangeListener(ChangeListener cl) {
				wrapped.removeChangeListener(cl);
			}

			public void removeChangeListener(ChangeListener cl, ChangeType ct) {
				wrapped.removeChangeListener(cl, ct);
			}

			public boolean isUnchanging(ChangeType ct) {
				return wrapped.isUnchanging(ct);
			}

			public Annotation getAnnotation() {
				return wrapped.getAnnotation();
			}
		}
	}
}
