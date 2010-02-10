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

package org.biojava.bio.gui.sequence;

import java.awt.Graphics2D;
import java.awt.event.MouseEvent;
import java.io.Serializable;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.biojava.bio.BioError;
import org.biojava.bio.seq.FeatureFilter;
import org.biojava.bio.seq.FeatureHolder;
import org.biojava.utils.AbstractChangeable;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeForwarder;
import org.biojava.utils.ChangeListener;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.Changeable;
import org.biojava.utils.cache.CacheMap;
import org.biojava.utils.cache.FixedSizeMap;

/**
 * <p><code>PairwiseFilteringRenderer</code> wraps a
 * <code>PairwiseSequenceRenderer</code> and filters the
 * <code>PairwiseRenderContext</code>s passed to it. The renderer
 * receives a new <code>PairwiseRenderContext</code> which has had
 * both of its <code>FeatureHolder</code>s filtered with the
 * <code>FeatureFilter</code>.<p>
 *
 * <p>Instances of this class cache up to 5 of the derived
 * <code>PairwiseRenderContext</code>s. Therefore cycling through up
 * to 5 different <code>FeatureFilter</code>s will only be hitting the
 * cache rather than recalculating everthing. Should the
 * <code>FeatureHolder</code>s themselves change, the cache will be
 * flushed. As only the features overlapping the context's range are
 * filtered, changing the range will also result in re-filtering.</p>
 *
 * @author Keith James
 * @author Matthew Pocock
 * @since 1.2
 */
public class PairwiseFilteringRenderer extends AbstractChangeable
    implements PairwiseSequenceRenderer, Serializable
{
    /**
     * Constant <code>FILTER</code> indicating a change to the
     * renderer's filter.
     */
    public static final ChangeType FILTER =
        new ChangeType("The filter has changed",
                       "org.biojava.bio.gui.sequence.PairwiseFilteringRenderer",
                       "FILTER",
                       SequenceRenderContext.LAYOUT);

    /**
     * Constant <code>RECURSE</code> indicating a change to the
     * renderer's filter recursion flag.
     */
    public static final ChangeType RECURSE =
        new ChangeType("The recurse flag has changed",
                       "org.biojava.bio.gui.sequence.PairwiseFilteringRenderer",
                       "RECURSE",
                       SequenceRenderContext.LAYOUT);

    /**
     * Constant <code>RENDERER</code> indicating a change to the
     * renderer.
     */
    public static final ChangeType RENDERER =
        new ChangeType("The renderer has changed",
                       "org.biojava.bio.gui.sequence.PairwiseFilteringRenderer",
                       "RENDERER",
                       SequenceRenderContext.REPAINT);

    /**
     * <code>filter</code> is the filter applied to both
     * <code>FeatureHolder</code>s.
     */
    protected FeatureFilter filter;
    /**
     * <code>recurse</code> indicates whether the filter should
     * recurse through any subfeatures.
     */
    protected boolean recurse;

    // Cache of previously created subcontexts
    private CacheMap subContextCache = new FixedSizeMap(5);
    // Set of listeners to cache keys
    private Set        cacheFlushers = new HashSet();

    // Delegate renderer
    private PairwiseSequenceRenderer  renderer;
    private transient ChangeForwarder rendererForwarder;

    /**
     * Creates a new <code>PairwiseFilteringRenderer</code> which uses
     * a filter which accepts all features.
     *
     * @param renderer a <code>PairwiseSequenceRenderer</code>.
     */
    public PairwiseFilteringRenderer(PairwiseSequenceRenderer renderer)
    {
        this(renderer, FeatureFilter.all, false);
    }

    /**
     * Creates a new <code>PairwiseFilteringRenderer</code>.
     *
     * @param renderer a <code>PairwiseSequenceRenderer</code>.
     * @param filter a <code>FeatureFilter</code>.
     * @param recurse a <code>boolean</code>.
     */
    public PairwiseFilteringRenderer(PairwiseSequenceRenderer renderer,
                                     FeatureFilter            filter,
                                     boolean                  recurse)
    {
        try
        {
            this.renderer = renderer;
            setFilter(filter);
            setRecurse(recurse);
        }
        catch (ChangeVetoException cve)
        {
            throw new BioError("Assertion failed: should have no listeners", cve);
        }
    }

    protected ChangeSupport getChangeSupport(ChangeType ct)
    {
        ChangeSupport cs = super.getChangeSupport(ct);

        if (rendererForwarder == null)
        {
            rendererForwarder =
                new PairwiseSequenceRenderer.PairwiseRendererForwarder(this, cs);

            if (renderer instanceof Changeable)
            {
                Changeable c = (Changeable) renderer;
                c.addChangeListener(rendererForwarder,
                                    SequenceRenderContext.REPAINT);
            }
        }

        return cs;
    }

    /**
     * <code>getRenderer</code> return the current renderer.
     *
     * @return a <code>PairwiseSequenceRenderer</code>.
     */
    public PairwiseSequenceRenderer getRenderer()
    {
        return renderer;
    }

    /**
     * <code>setRenderer</code> sets the renderer.
     *
     * @param renderer a <code>PairwiseSequenceRenderer</code>.
     *
     * @exception ChangeVetoException if the change is vetoed.
     */
    public void setRenderer(PairwiseSequenceRenderer renderer)
        throws ChangeVetoException
    {
        if (hasListeners())
        {
            ChangeEvent ce = new ChangeEvent(this, RENDERER,
                                             renderer, this.renderer);

            ChangeSupport cs = getChangeSupport(RENDERER);
            synchronized(cs)
            {
                cs.firePreChangeEvent(ce);

                if (this.renderer instanceof Changeable)
                {
                    Changeable c = (Changeable) this.renderer;
                    c.removeChangeListener(rendererForwarder);
                }

                this.renderer = renderer;

                if (renderer instanceof Changeable)
                {
                    Changeable c = (Changeable) renderer;
                    c.addChangeListener(rendererForwarder);
                }
                cs.firePostChangeEvent(ce);
            }
        }
        else
        {
            this.renderer = renderer;
        }
    }

    /**
     * <code>getFilter</code> returns the current filter.
     *
     * @return a <code>FeatureFilter</code>.
     */
    public FeatureFilter getFilter()
    {
        return filter;
    }

    /**
     * <code>setFilter</code> sets the filter.
     *
     * @param filter a <code>FeatureFilter</code>.
     *
     * @exception ChangeVetoException if the change is vetoed.
     */
    public void setFilter(FeatureFilter filter)
        throws ChangeVetoException
    {
        if (hasListeners())
        {
            ChangeSupport cs = getChangeSupport(FILTER);
            synchronized(cs)
            {
                ChangeEvent ce = new ChangeEvent(this, FILTER,
                                                 this.filter, filter);
                cs.firePreChangeEvent(ce);
                this.filter = filter;
                cs.firePostChangeEvent(ce);
            }
        }
        else
        {
            this.filter = filter;
        }
    }

    /**
     * <code>getRecurse</code> returns the recursion flag of the
     * filter.
     *
     * @return a <code>boolean</code>.
     */
    public boolean getRecurse()
    {
        return recurse;
    }

    /**
     * <code>setRecurse</code> sets the recursion flag on the filter.
     *
     * @param recurse a <code>boolean</code>.
     *
     * @exception ChangeVetoException if the change is vetoed.
     */
    public void setRecurse(boolean recurse) throws ChangeVetoException
    {
        if (hasListeners())
        {
            ChangeSupport cs = getChangeSupport(RECURSE);
            synchronized(cs)
            {
                ChangeEvent ce = new ChangeEvent(this, RECURSE,
                                                 new Boolean(recurse),
                                                 new Boolean(this.recurse));
                cs.firePreChangeEvent(ce);
                this.recurse = recurse;
                cs.firePostChangeEvent(ce);
            }
        }
        else
        {
            this.recurse = recurse;
        }
    }

    public void paint(Graphics2D g2, PairwiseRenderContext context)
    {
        renderer.paint(g2, getSubContext(context));
    }

    public SequenceViewerEvent processMouseEvent(PairwiseRenderContext context,
                                                 MouseEvent            me,
                                                 List                  path)
    {
        path.add(this);
        return renderer.processMouseEvent(getSubContext(context), me, path);
    }

    /**
     * <code>getSubContext</code> creates a new context which has
     * <code>FeatureHolder</code>s filtered using the current filter.
     *
     * @param context a <code>PairwiseRenderContext</code>.
     *
     * @return a <code>PairwiseRenderContext</code>.
     */
    protected PairwiseRenderContext getSubContext(PairwiseRenderContext context)
    {
        // Filter the sequence features
        FeatureFilter ff =
            new FeatureFilter.And(filter,
                                  new FeatureFilter.OverlapsLocation(context.getRange()));

        // Filter the secondary sequence features
        FeatureFilter ffSec =
            new FeatureFilter.And(filter,
                                  new FeatureFilter.OverlapsLocation(context.getSecondaryRange()));

        // Create a cache key
        FilteredSubContext cacheKey = new FilteredSubContext(context,
                                                             ff,
                                                             ffSec,
                                                             recurse);
        // Try the cache for a context first
        PairwiseRenderContext filtered =
            (PairwiseRenderContext) subContextCache.get(cacheKey);

        if (filtered == null)
        {
            // None in cache, so make a new one
            filtered =
                new SubPairwiseRenderContext(context, // context delegate
                                             null,    // symbols
                                             null,    // secondary symbols
                                             context.getFeatures().filter(ff,    recurse),
                                             context.getFeatures().filter(ffSec, recurse),
                                             null,    // range
                                             null);   // secondary range

            // Add to cache
            subContextCache.put(cacheKey, filtered);
            // Create a listener for changes in the feature holder
            CacheFlusher cf = new CacheFlusher(cacheKey);
            // Add the listener for changes in features
            ((Changeable) context.getSymbols())
                .addChangeListener(cf, FeatureHolder.FEATURES);
            cacheFlushers.add(cf);
        }

        return filtered;
    }

    /**
     * <code>FilteredSubContext</code> is a cache key whose equality
     * with another such key is determined by context, filter and
     * recurse values.
     */
    private class FilteredSubContext
    {
        private PairwiseRenderContext context;
        private FeatureFilter         filter;
        private FeatureFilter         secondaryFilter;
        private boolean               recurse;

        public FilteredSubContext(PairwiseRenderContext context,
                                  FeatureFilter         filter,
                                  FeatureFilter         secondaryFilter,
                                  boolean               recurse)
        {
            this.context         = context;
            this.filter          = filter;
            this.secondaryFilter = secondaryFilter;
            this.recurse         = recurse;
        }

        public boolean equals(Object o)
        {
            if (! (o instanceof FilteredSubContext))
            {
                return false;
            }

            FilteredSubContext that = (FilteredSubContext) o;

            return
                context.equals(that.context) &&
                filter.equals(that.filter) &&
                secondaryFilter.equals(that.secondaryFilter) &&
                (recurse == that.recurse);
        }

        public int hashCode()
        {
            return context.hashCode() ^ filter.hashCode();
        }
    }

    /**
     * <code>CacheFlusher</code> is a listener for changes to the
     * original context's feature holders. A change forces removal of
     * its associated cache key from the cache.
     */
    private class CacheFlusher implements ChangeListener
    {
        private FilteredSubContext fsc;

        public CacheFlusher(FilteredSubContext fsc)
        {
            this.fsc = fsc;
        }

        public void preChange(ChangeEvent ce) { }

        public void postChange(ChangeEvent ce)
        {
            subContextCache.remove(fsc);
            cacheFlushers.remove(this);

            if (hasListeners())
            {
                ChangeSupport cs = getChangeSupport(SequenceRenderContext.LAYOUT);
                synchronized(cs)
                {
                    ChangeEvent ce2 =
                        new ChangeEvent(PairwiseFilteringRenderer.this,
                                        SequenceRenderContext.LAYOUT);
                    cs.firePostChangeEvent(ce2);
                }
            }
        }
    }
}
