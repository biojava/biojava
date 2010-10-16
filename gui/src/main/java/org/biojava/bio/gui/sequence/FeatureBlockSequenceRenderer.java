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
import java.awt.Shape;
import java.awt.event.MouseEvent;
import java.awt.geom.AffineTransform;
import java.util.Iterator;
import java.util.List;

import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.FeatureFilter;
import org.biojava.bio.seq.FeatureHolder;
import org.biojava.bio.seq.FilterUtils;
import org.biojava.bio.symbol.RangeLocation;
import org.biojava.utils.AbstractChangeable;
import org.biojava.utils.AssertionFailure;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeForwarder;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.Changeable;

/**
 * <code>FeatureBlockSequenceRenderer</code> forms a bridge between
 * <code>Sequence</code> rendering and <code>Feature</code>
 * rendering. It is a <code>SequenceRenderer</code> which iterates
 * through a <code>Sequence</code>'s <code>Feature</code>s and makes
 * method calls on a <code>FeatureRenderer</code>.
 *
 * @author Matthew Pocock
 * @author Keith James
 * @author Thomas Down
 */
public class FeatureBlockSequenceRenderer extends AbstractChangeable
    implements SequenceRenderer {
    public static ChangeType FEATURE_RENDERER =
        new ChangeType("The associated FeatureRenderer has changed",
                       "org.biojava.bio.gui.sequence.FeatureBlockSequenceRenderer",
                       "FEATURE_RENDERER",
                       SequenceRenderContext.LAYOUT);
    public static ChangeType FEATURE_COLLAPSING =
        new ChangeType("Changed whether the render collapses when no features are visible",
                       "org.biojava.bio.gui.sequence.FeatureBlockSequenceRenderer",
                       "FEATURE_COLLAPSING",
                       SequenceRenderContext.LAYOUT);

    private FeatureRenderer renderer;
    private boolean isCollapsing = true;;
    private transient ChangeForwarder rendForwarder;

    protected ChangeSupport getChangeSupport(ChangeType ct) {
        ChangeSupport cs = super.getChangeSupport(ct);

        if (rendForwarder == null) {
            rendForwarder = new SequenceRenderer.RendererForwarder(this, cs);
            if ((renderer != null) && (renderer instanceof Changeable)) {
                Changeable c = (Changeable) this.renderer;
                c.addChangeListener(rendForwarder,
                                    SequenceRenderContext.REPAINT);
            }
        }

        return cs;
    }

    /**
     * Creates a new <code>FeatureBlockSequenceRenderer</code> which
     * uses a <code>BasicFeatureRenderer</code> as its renderer.
     */
    public FeatureBlockSequenceRenderer() {
        try {
            setFeatureRenderer(new BasicFeatureRenderer());
        } catch (ChangeVetoException cve) {
            throw new AssertionFailure("Assertion Failure: Should have no listeners", cve);
        }
    }

    /**
     * Creates a new <code>FeatureBlockSequenceRenderer</code> which
     * uses the specified <code>FeatureRenderer</code>.
     *
     * @param fRend a <code>FeatureRenderer</code>.
     */
    public FeatureBlockSequenceRenderer(FeatureRenderer fRend) {
        try {
            setFeatureRenderer(fRend);
        } catch (ChangeVetoException cve) {
            throw new AssertionFailure("Assertion Failure: Should have no listeners", cve);
        }
    }

    /**
     * <code>getFeatureRenderer</code> returns the currently active
     * renderer.
     *
     * @return a <code>FeatureRenderer</code>.
     */
    public FeatureRenderer getFeatureRenderer() {
        return renderer;
    }

    /**
     * <code>setFeatureRenderer</code> sets the renderer to be used.
     *
     * @param renderer a <code>FeatureRenderer</code>.
     * @exception ChangeVetoException if the renderer can not be
     * changed.
     */
    public void setFeatureRenderer(FeatureRenderer renderer)
        throws ChangeVetoException {
        if (hasListeners()) {
            ChangeSupport cs = getChangeSupport(FEATURE_RENDERER);
            synchronized(cs) {
                ChangeEvent ce = new ChangeEvent(this,
                                                 FEATURE_RENDERER,
                                                 this.renderer,
                                                 renderer);
                cs.firePreChangeEvent(ce);
                if ((this.renderer != null) &&
                    (this.renderer instanceof Changeable)) {
                    Changeable c = (Changeable) this.renderer;
                    c.removeChangeListener(rendForwarder, SequenceRenderContext.REPAINT);
                }
                this.renderer = renderer;
                if (renderer instanceof Changeable) {
                    Changeable c = (Changeable) renderer;
                    c.removeChangeListener(rendForwarder, SequenceRenderContext.REPAINT);
                }
                cs.firePostChangeEvent(ce);
            }
        } else {
            this.renderer = renderer;
        }
    }

    /**
     * Specifies if the renderer should collapse to zero depth when no
     * features are visible (default <code>true</code>).
     *
     * @since 1.3
     */

    public void setCollapsing(boolean b)
        throws ChangeVetoException
    {
        if (hasListeners()) {
            ChangeSupport cs = getChangeSupport(FEATURE_COLLAPSING);
            synchronized (cs) {
                ChangeEvent ce = new ChangeEvent(this,
                                                 FEATURE_COLLAPSING,
                                                 new Boolean(this.isCollapsing),
                                                 new Boolean(b));
                cs.firePreChangeEvent(ce);
                this.isCollapsing = b;
                cs.firePostChangeEvent(ce);
            }
        } else {
            this.isCollapsing = b;
        }
    }

    /**
     * Returns <code>true</code> if this class collapses to zero depth when there are
     * no visible features.
     *
     * @since 1.3
     */

    public boolean getCollapsing() {
        return isCollapsing;
    }

    public double getDepth(SequenceRenderContext src) {
        FeatureHolder features = src.getFeatures();
        FeatureFilter filter =
            FilterUtils.overlapsLocation(src.getRange());
        FeatureHolder fh = features.filter(filter, false);
        if (!isCollapsing || fh.countFeatures() > 0) {
            return renderer.getDepth(src);
        } else {
            return 0.0;
        }
    }

    public double getMinimumLeader(SequenceRenderContext src) {
        return 0.0;
    }

    public double getMinimumTrailer(SequenceRenderContext src) {
        return 0.0;
    }

  public void paint(Graphics2D g, SequenceRenderContext src) {
    FeatureFilter filt = FilterUtils.shadowOverlapsLocation(
            GUITools.getVisibleRange(src, g) );
    FeatureHolder feats = src.getFeatures().filter(filt, false);

    for (Iterator i =feats.features(); i.hasNext();) {
      Shape clip = g.getClip();
      AffineTransform at = g.getTransform();

      Feature f = (Feature) i.next();
      renderer.renderFeature(g, f, src);

      g.setTransform(at);
      g.setClip(clip);
    }
  }

    public SequenceViewerEvent processMouseEvent(SequenceRenderContext src,
                                                 MouseEvent me,
                                                 List path) {
        double pos;
        if (src.getDirection() == SequenceRenderContext.HORIZONTAL) {
            pos = me.getPoint().getX();
        } else {
            pos = me.getPoint().getY();
        }

        int sMin = src.graphicsToSequence(pos);
        int sMax = src.graphicsToSequence(pos + 1);

        FeatureHolder hits =
            src.getFeatures().filter(new FeatureFilter.OverlapsLocation(new RangeLocation(sMin, sMax)),
                                     false);

        hits = renderer.processMouseEvent(hits, src, me);

        return new SequenceViewerEvent(this, hits, sMin, me, path);
    }
}

