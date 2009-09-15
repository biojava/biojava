package org.biojava.bio.seq.projection;



import java.io.Serializable;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import org.biojava.bio.BioException;
import org.biojava.bio.BioRuntimeException;
import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.FeatureFilter;
import org.biojava.bio.seq.FeatureHolder;
import org.biojava.bio.seq.FilterUtils;
import org.biojava.bio.seq.Sequence;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeListener;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.Unchangeable;



/**
 * A good base class to implement ProjectionContext from.
 *

 * <p>

 * ReparentContext takes care of much of the ugliness of implementing

 * ProjectionContext, such as handling listeners and grafting features onto

 * a new parent. It also sets up a framework for mutating feature filters.

 * Think carefully before overriding methods in this class.

 * </p>

 *

 * @author Matthew Pocock

 * @author Thomas Down

 */

public class ReparentContext implements ProjectionContext, Serializable {
	private transient final Map forwardersByFeature = new HashMap();
	private final FeatureHolder parent;
	private final FeatureHolder wrapped;



	public ReparentContext(FeatureHolder parent,
			FeatureHolder wrapped) 
	{
		this.parent = parent;
		this.wrapped = wrapped;
	}

	public final FeatureHolder getParent() {
		return parent;
	}

	public final FeatureHolder getUnprojectedFeatures() {
		return wrapped;
	}

	/**
	 * Create a single projected feature using the rules of this <code>ProjectedFeatureHolder</code>.
	 */

	public Feature projectFeature(Feature feat) {
		return ProjectionEngine.DEFAULT.projectFeature(feat, this);
	}

	public Feature revertFeature(Feature feat) {
		return ((Projection) feat).getViewedFeature();
	}

	public final FeatureFilter projectFilter(FeatureFilter ff) {
		return FilterUtils.transformFilter(ff, getTransformer());
	}

	public final FeatureFilter revertFilter(FeatureFilter ff) {
		return FilterUtils.transformFilter(ff, getReverter());
	}

	protected FilterUtils.FilterTransformer getTransformer() {
		return new FilterUtils.FilterTransformer() {
			public FeatureFilter transform(FeatureFilter ff) {
				// fixme: should we be mapping filters on feature instance or sequence
				// instance?

				return ff;
			}
		};
	}

	protected FilterUtils.FilterTransformer getReverter() {
		return new FilterUtils.FilterTransformer() {
			public FeatureFilter transform(FeatureFilter ff) {
				// fixme: should we be mapping filters on feature instance or sequence
				// instance?

				return ff;
			}
		};
	}

	public final FeatureHolder getParent(Feature f) {
		FeatureHolder oldP = f.getParent();
		if (oldP instanceof Feature) {
			if (wrapped.containsFeature(f)) {
				return parent;
			} else {
				return projectFeature((Feature) oldP);
			}
		} else {
			return parent;
		}
	}

	public final Sequence getSequence(Feature f) {
		FeatureHolder fh = getParent();
		while (fh instanceof Feature) {
			fh = ((Feature) fh).getParent();
		}

		if (! (fh instanceof Sequence)) {
			throw new BioRuntimeException("Chasing up parents to get sequence: actually hit: " + fh.toString());
		}

		return (Sequence) fh;
	}

	public FeatureHolder projectChildFeatures(Feature f, FeatureHolder parent) {
		return new ProjectionSet(f);
	}

	public final Feature createFeature(Feature.Template projTempl)
	throws BioException, ChangeVetoException
	{

		Feature f = null;

		synchronized(wrapped){
			f = wrapped.createFeature(
					ProjectionEngine.DEFAULT.revertTemplate(projTempl, this));
		}

		f = projectFeature(f);





		return f;


	}
	public final void removeFeature(Feature dyingChild)
	throws BioException, ChangeVetoException 
	{
		wrapped.removeFeature(revertFeature(dyingChild));
	}

	public final Feature createFeature(Feature f, Feature.Template projTempl)
	throws BioException, ChangeVetoException
	{



		Feature f1 = null;
		synchronized (f) {
			f1 = revertFeature(f);
		}

		Feature f2 = null;
		synchronized(f1){
			f2 = f1.createFeature(ProjectionEngine.DEFAULT.revertTemplate(projTempl, this));
		}

		Feature newF = null;

		synchronized (f2) {
			newF = projectFeature(f2);
		}

		return newF;


	}

	public final void removeFeature(Feature f, Feature f2)
	throws ChangeVetoException, BioException 
	{
		revertFeature(f).removeFeature(revertFeature(f2));
	}

	public final FeatureFilter getSchema(Feature f) {
		return projectFilter(f.getSchema());
	}



	public final void addChangeListener(Feature f, ChangeListener cl, ChangeType ct) {
		if (!f.isUnchanging(ct)) {
			PFChangeForwarder forwarder = (PFChangeForwarder) forwardersByFeature.get(f);
			if (forwarder == null) {
				forwarder = new PFChangeForwarder(f);
				forwardersByFeature.put(f, forwarder);
				f.addChangeListener(forwarder, ChangeType.UNKNOWN);
			}

			forwarder.addChangeListener(cl, ct);
		}
	}

	public final void removeChangeListener(Feature f, ChangeListener cl, ChangeType ct) {
		PFChangeForwarder forwarder = (PFChangeForwarder) forwardersByFeature.get(f);
		if (forwarder != null) {
			forwarder.removeChangeListener(cl, ct);
			if (!forwarder.hasListeners()) {
				forwardersByFeature.remove(f);
				f.removeChangeListener(forwarder, ChangeType.UNKNOWN);
			}
		}
	}

	public final FeatureHolder projectFeatures(FeatureHolder fh) {
		return new ProjectionSet(fh);
	}

	//
	// Dumb set of features to which we delegate everything except the
	// ChangeEvent stuff.
	//

	private class ProjectionSet
	extends Unchangeable
	implements FeatureHolder, Serializable
	{
		private final FeatureHolder baseSet;

		ProjectionSet(FeatureHolder baseSet) {
			this.baseSet = baseSet;
		}

		public int countFeatures() {
			return baseSet.countFeatures();
		}

		public Iterator features() {
			final Iterator wrappedIterator = baseSet.features();

			return new Iterator() {
				public boolean hasNext() {
					return wrappedIterator.hasNext();
				}

				public Object next() {
					return projectFeature((Feature) wrappedIterator.next());
				}

				public void remove() {
					throw new UnsupportedOperationException();
				}
			};
		}

		public boolean containsFeature(Feature f) { 
			if (! (f instanceof Projection)) {
				return false;
			} else {
				Projection p = (Projection) f;
				return p.getProjectionContext() == ReparentContext.this && baseSet.containsFeature(p.getViewedFeature());
			}
		}

		public FeatureHolder filter(FeatureFilter ff) {
			return filter(ff, true); // bit of a hack for now.
		}

		public FeatureHolder filter(FeatureFilter ff, boolean recurse) {
			ff = revertFilter(ff);
			FeatureHolder toProject = baseSet.filter(ff, recurse);
			return new ProjectionSet(toProject);
		}

		public Feature createFeature(Feature.Template templ)
		throws ChangeVetoException, BioException 
		{
			throw new ChangeVetoException("Can't create features in this projection");
		}

		public void removeFeature(Feature toDie)
		throws ChangeVetoException, BioException 
		{
			throw new ChangeVetoException("Can't remove features in this projection");
		}

		public FeatureFilter getSchema() {
			return projectFilter(baseSet.getSchema());
		}
	}

	private class PFChangeForwarder
	extends ChangeSupport
	implements ChangeListener
	{
		private Feature master;

		public PFChangeForwarder(Feature master) {
			super(1);
			this.master = master;
		}

		public void preChange(ChangeEvent cev)
		throws ChangeVetoException 
		{
			ChangeEvent cev2 = forwardFeatureChangeEvent(master, cev);
			if (cev2 != null) {
				synchronized (cev2) {
					firePreChangeEvent(cev2);
				}

			}
		}

		public void postChange(ChangeEvent cev) {
			ChangeEvent cev2 = forwardFeatureChangeEvent(master, cev);
			if (cev2 != null) {
				synchronized (cev2) {
					firePostChangeEvent(cev2);
				}
			}
		}

		/**
		 * Called internally to generate a forwarded version of a ChangeEvent from a ProjectedFeature
		 *
		 * @param f the feature who's projection is due to receive an event.
		 * @return a tranformed event, or <code>null</code> to cancel the event.
		 */

		protected ChangeEvent forwardFeatureChangeEvent(Feature f, ChangeEvent cev) {
			return new ChangeEvent(projectFeature(f),
					cev.getType(),
					cev.getChange(),
					cev.getPrevious(),
					cev);

		}
	}
}

