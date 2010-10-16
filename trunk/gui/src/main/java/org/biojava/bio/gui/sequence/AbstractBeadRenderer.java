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

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Paint;
import java.awt.Stroke;
import java.awt.event.MouseEvent;
import java.io.Serializable;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.FeatureFilter;
import org.biojava.bio.seq.FeatureHolder;
import org.biojava.bio.seq.OptimizableFilter;
import org.biojava.utils.AbstractChangeable;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;

/**
 * <p><code>AbstractBeadRenderer</code> is a an abstract base class
 * for the creation of <code>FeatureRenderer</code>s which use a
 * 'string of beads' metaphor for displaying features. Each subclass
 * of <code>AbstractBeadRenderer</code> should override the abstract
 * method <code>renderBead()</code> and provide the drawing routine
 * for its particular bead type.</p>
 *
 * <p>A concrete <code>BeadFeatureRenderer</code> may render a series
 * of features in more than one style by delegating to other
 * <code>BeadFeatureRenderer</code>s for the additional style(s). This
 * is achieved using the <code>setDelegateRenderer()</code> method
 * which associates an <code>OptimizableFilter</code> with another
 * <code>BeadFeatureRenderer</code>. Any feature accepted by the
 * filter is rendered with that renderer, while the remainder are
 * rendered by the current renderer.</p>
 *
 * @author Keith James
 * @author Paul Seed
 * @since 1.2
 */
public abstract class AbstractBeadRenderer extends AbstractChangeable
    implements BeadFeatureRenderer, Serializable
{
    /**
     * Constant <code>DISPLACEMENT</code> indicating a change to the
     * Y-axis displacement of the features.
     */
    public static final ChangeType DISPLACEMENT =
        new ChangeType("The displacement of the features has changed",
                       "org.biojava.bio.gui.sequence.AbstractBeadRenderer",
                       "DISPLACEMENT", SequenceRenderContext.LAYOUT);
    
    /**
     * Constant <code>DEPTH</code> indicating a change to the depth of
     * the renderer.
     */
    public static final ChangeType DEPTH =
        new ChangeType("The depth of the renderer has changed",
                       "org.biojava.bio.gui.sequence.AbstractBeadRenderer",
                       "DEPTH", SequenceRenderContext.LAYOUT);

    /**
     * Constant <code>OUTLINE</code> indicating a change to the
     * outline paint of the features.
     */
    public static final ChangeType OUTLINE =
        new ChangeType("The outline of the features has changed",
                       "org.biojava.bio.gui.sequence.AbstractBeadRenderer",
                       "OUTLINE", SequenceRenderContext.REPAINT);

    /**
     * Constant <code>STROKE</code> indicating a change to the outline
     * stroke of the features.
     */
    public static final ChangeType STROKE =
        new ChangeType("The stroke of the features has changed",
                       "org.biojava.bio.gui.sequence.AbstractBeadRenderer",
                       "STROKE", SequenceRenderContext.REPAINT);

    /**
     * Constant <code>FILL</code> indicating a change to the fill of
     * the features.
     */
    public static final ChangeType FILL =
        new ChangeType("The fill of the features has changed",
                       "org.biojava.bio.gui.sequence.AbstractBeadRenderer",
                       "FILL", SequenceRenderContext.REPAINT);

    protected double        beadDepth;
    protected double beadDisplacement;
    protected Paint       beadOutline;
    protected Paint          beadFill;
    protected Stroke       beadStroke;

    // Map of FeatureFilter -> FeatureRenderer
    protected Map delegates;
    // Map of Feature -> FeatureRenderer
    protected Map delegationCache;

    /**
     * Creates a new <code>AbstractBeadRenderer</code> with no
     * delegates. It will render all features itself, using its own
     * style settings.
     */
    public AbstractBeadRenderer()
    {
        this(10.0f, 0.0f, Color.black, Color.black, new BasicStroke());
    }

    /**
     * Creates a new <code>AbstractBeadRenderer</code> object.
     *
     * @param beadDepth a <code>double</code>.
     * @param beadDisplacement a <code>double</code>.
     * @param beadOutline a <code>Paint</code>.
     * @param beadFill a <code>Paint</code>.
     * @param beadStroke a <code>Stroke</code>.
     */
    public AbstractBeadRenderer(double beadDepth,
                                double beadDisplacement,
                                Paint  beadOutline,
                                Paint  beadFill,
                                Stroke beadStroke)
    {
        this.beadDepth        = beadDepth;
        this.beadDisplacement = beadDisplacement;
        this.beadOutline      = beadOutline;
        this.beadFill         = beadFill;
        this.beadStroke       = beadStroke;

        delegates       = new HashMap();
        delegationCache = new HashMap();
    }

    /**
     * <code>processMouseEvent</code> defines the behaviour on
     * revieving a mouse event.
     *
     * @param holder a <code>FeatureHolder</code>.
     * @param context a <code>SequenceRenderContext</code>.
     * @param mEvent a <code>MouseEvent</code>.
     *
     * @return a <code>FeatureHolder</code>.
     */
    public FeatureHolder processMouseEvent(FeatureHolder         holder,
                                           SequenceRenderContext context,
                                           MouseEvent            mEvent)
    {
        return holder;
    }

    /**
     * <code>renderFeature</code> draws a feature using the supplied
     * graphics context. The rendering may be delegated to another
     * <code>FeatureRenderer</code> instance.
     *
     * @param g2 a <code>Graphics2D</code> context.
     * @param f a <code>Feature</code> to render.
     * @param context a <code>SequenceRenderContext</code> context.
     */
    public void renderFeature(Graphics2D            g2,
                              Feature               f,
                              SequenceRenderContext context)
    {
        // Check the cache first
        if (delegationCache.containsKey(f))
        {
            // System.err.println("Used cache for: " + f);

            BeadFeatureRenderer cachedRenderer =
                (BeadFeatureRenderer) delegationCache.get(f);

            cachedRenderer.renderBead(g2, f, context);
            return;
        }

        for (Iterator di = delegates.keySet().iterator(); di.hasNext();)
        {
            FeatureFilter filter = (FeatureFilter) di.next();

            if (filter.accept(f))
            {
                // System.err.println(filter + " accepted " + f);

                FeatureRenderer delegate =
                    (FeatureRenderer) delegates.get(filter);

                delegate.renderFeature(g2, f, context);
                return;
            }
        }

        delegationCache.put(f, this);
        // System.err.println("Rendering: " + f);
        renderBead(g2, f, context);
    }

    /**
     * <code>setDelegateRenderer</code> associates an
     * <code>OptimizableFilter</code> with a
     * <code>BeadFeatureRenderer</code>. Any feature accepted by the
     * filter will be passed to the associated renderer for
     * drawing. The <code>OptimizableFilter</code>s should be disjoint
     * with respect to each other (a feature may not be rendered more
     * than once).
     *
     * @param filter an <code>OptimizableFilter</code>.
     * @param renderer a <code>BeadFeatureRenderer</code>.
     *
     * @exception IllegalArgumentException if the filter is not
     * disjoint with existing delegate filters.
     */
    public void setDelegateRenderer(OptimizableFilter   filter,
                                    BeadFeatureRenderer renderer)
        throws IllegalArgumentException
    {
        // Ensure the cache doesn't hide the new delegate
        delegationCache.clear();
	
        Set delegateFilters = delegates.keySet();

        if (delegateFilters.size() == 0)
        {
            delegates.put(filter, renderer);
        }
        else
        {
            for (Iterator fi = delegateFilters.iterator(); fi.hasNext();)
            {
                OptimizableFilter thisFilter = (OptimizableFilter) fi.next();

                if (! thisFilter.isDisjoint(filter))
                {
                    throw new IllegalArgumentException("Failed to apply filter as it clashes with existing filter "
                                                       + thisFilter
                                                       + " (filters must be disjoint)");
                }
                else
                {
                    delegates.put(filter, renderer);
                    break;
                }
            }
        }
    }

    /**
     * <code>removeDelegateRenderer</code> removes any association
     * of the given <code>OptimizableFilter</code> with a
     * <code>BeadFeatureRenderer</code>.
     *
     * @param filter an <code>OptimizableFilter</code>.
     */
    public void removeDelegateRenderer(OptimizableFilter   filter)
    {
        // Ensure the cache doesn't hide the change of delegate
        delegationCache.clear();
        delegates.remove(filter);
    }

    /**
     * <code>getDepth</code> calculates the depth required by this
     * renderer to display its beads. It recurses through its delegate
     * renderers and returns the highest value. Concrete renderers
     * should override this method and supply code to calculate their
     * own depth. If a subclass needs to know the depth of its
     * delegates (as is likely if it has any) they can call this
     * method using <code>super.getDepth()</code>.
     *
     * @param context a <code>SequenceRenderContext</code>.
     *
     * @return a <code>double</code>.
     */
    public double getDepth(SequenceRenderContext context)
    {
        Collection delegateRenderers = delegates.values();
        double maxDepth = 0.0;

        if (delegateRenderers.size() == 0)
        {
            return maxDepth + 1.0;
        }
        else
        {
            for (Iterator ri = delegateRenderers.iterator(); ri.hasNext();)
            {
                maxDepth = Math.max(maxDepth, ((FeatureRenderer) ri.next()).getDepth(context));
            }

            return maxDepth + 1.0;
        }
    }

    /**
     * <code>getBeadDepth</code> returns the depth of a single bead
     * produced by this renderer.
     *
     * @return a <code>double</code>.
     */
    public double getBeadDepth()
    {
        return beadDepth;
    }

    /**
     * <code>setBeadDepth</code> sets the depth of a single bead
     * produced by this renderer.
     *
     * @param depth a <code>double</code>.
     *
     * @exception ChangeVetoException if an error occurs.
     */
    public void setBeadDepth(double depth) throws ChangeVetoException
    {
        if (hasListeners())
        {
            ChangeSupport cs = getChangeSupport(SequenceRenderContext.LAYOUT);
            synchronized(cs)
            {
                ChangeEvent ce = new ChangeEvent(this, SequenceRenderContext.LAYOUT,
                                                 null, null,
                                                 new ChangeEvent(this, DEPTH,
                                                                 new Double(beadDepth),
                                                                 new Double(depth)));
                cs.firePreChangeEvent(ce);
                beadDepth = depth;
                cs.firePostChangeEvent(ce);
            }
        }
        else
        {
            beadDepth = depth;
        }
    }

    /**
     * <code>getBeadDisplacement</code> returns the displacement of
     * beads from the centre line of the renderer. A positive value
     * indicates displacment downwards (for horizontal renderers) or
     * to the right (for vertical renderers).
     *
     * @return a <code>double</code>.
     */
    public double getBeadDisplacement()
    {
        return beadDisplacement;
    }

    /**
     * <code>setBeadDisplacement</code> sets the displacement of
     * beads from the centre line of the renderer. A positive value
     * indicates displacment downwards (for horizontal renderers) or
     * to the right (for vertical renderers).
     *
     * @param displacement a <code>double</code>.
     *
     * @exception ChangeVetoException if an error occurs.
     */
    public void setBeadDisplacement(double displacement)
        throws ChangeVetoException
    {
        if (hasListeners())
        {
            ChangeSupport cs = getChangeSupport(SequenceRenderContext.LAYOUT);
            synchronized(cs)
            {
                ChangeEvent ce = new ChangeEvent(this, SequenceRenderContext.LAYOUT,
                                                 null, null,
                                                 new ChangeEvent(this, DISPLACEMENT,
                                                                 new Double(beadDisplacement),
                                                                 new Double(displacement)));
                cs.firePreChangeEvent(ce);
                beadDisplacement = displacement;
                cs.firePostChangeEvent(ce);
            }
        }
        else
        {
            beadDisplacement = displacement;
        }
    }

    /**
     * <code>getBeadOutline</code> returns the bead outline paint.
     *
     * @return a <code>Paint</code>.
     */
    public Paint getBeadOutline()
    {
        return beadOutline;
    }

    /**
     * <code>setBeadOutline</code> sets the bead outline paint.
     *
     * @param outline a <code>Paint</code>.
     *
     * @exception ChangeVetoException if an error occurs.
     */
    public void setBeadOutline(Paint outline) throws ChangeVetoException
    {
        if (hasListeners())
        {
            ChangeSupport cs = getChangeSupport(SequenceRenderContext.LAYOUT);
            synchronized(cs)
            {
                ChangeEvent ce = new ChangeEvent(this, SequenceRenderContext.LAYOUT,
                                                 null, null,
                                                 new ChangeEvent(this, OUTLINE,
                                                                 outline,
                                                                 beadOutline));
                cs.firePreChangeEvent(ce);
                beadOutline = outline;
                cs.firePostChangeEvent(ce);
            }
        }
        else
        {
            beadOutline = outline;
        }
    }

    /**
     * <code>getBeadStroke</code> returns the bead outline stroke.
     *
     * @return a <code>Stroke</code>.
     */
    public Stroke getBeadStroke()
    {
        return beadStroke;
    }

    /**
     * <code>setBeadStroke</code> sets the bead outline stroke.
     *
     * @param stroke a <code>Stroke</code>.
     *
     * @exception ChangeVetoException if an error occurs.
     */
    public void setBeadStroke(Stroke stroke) throws ChangeVetoException
    {
        if (hasListeners())
        {
            ChangeSupport cs = getChangeSupport(SequenceRenderContext.LAYOUT);
            synchronized(cs)
            {
                ChangeEvent ce = new ChangeEvent(this, SequenceRenderContext.LAYOUT,
                                                 null, null,
                                                 new ChangeEvent(this, STROKE,
                                                                 stroke,
                                                                 beadStroke));
                cs.firePreChangeEvent(ce);
                beadStroke = stroke;
                cs.firePostChangeEvent(ce);
            }
        }
        else
        {
            beadStroke = stroke;
        }
    } 

    /**
     * <code>getBeadFill</code> returns the bead fill paint.
     *
     * @return a <code>Paint</code>.
     */
    public Paint getBeadFill()
    {
        return beadFill;
    }

    /**
     * <code>setBeadFill</code> sets the bead fill paint.
     *
     * @param fill a <code>Paint</code>.
     *
     * @exception ChangeVetoException if an error occurs.
     */
    public void setBeadFill(Paint fill) throws ChangeVetoException
    {
        if (hasListeners())
        {
            ChangeSupport cs = getChangeSupport(SequenceRenderContext.LAYOUT);
            synchronized(cs)
            {
                ChangeEvent ce = new ChangeEvent(this, SequenceRenderContext.LAYOUT,
                                                 null, null,
                                                 new ChangeEvent(this, FILL,
                                                                 fill,
                                                                 beadFill));
                cs.firePreChangeEvent(ce);
                beadFill = fill;
                cs.firePostChangeEvent(ce);
            }
        }
        else
        {
            beadFill = fill;
        }
    }

    /**
     * <code>renderBead</code> should be overridden by the concrete
     * <code>BeadRenderer</code>.
     *
     * @param g2 a <code>Graphics2D</code>.
     * @param f a <code>Feature</code> to render.
     * @param context a <code>SequenceRenderContext</code> context.
     */
    public abstract void renderBead(Graphics2D            g2,
                                    Feature               f,
                                    SequenceRenderContext context);
}
