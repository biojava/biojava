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
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import org.biojava.utils.AbstractChangeable;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeForwarder;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.Changeable;

/**
 * <code>MultiLineRenderer</code> is a <code>SequenceRenderer</code>
 * which collects a number of other <code>SequenceRenderer</code>s
 * each of which render their own view of a <code>Sequence</code>.
 *
 * @author Matthew Pocock
 * @author Keith James
 */
public class MultiLineRenderer extends AbstractChangeable
    implements SequenceRenderer, Serializable {
    public static final ChangeType RENDERERS =
        new ChangeType("A SequenceRenderer has been added or removed.",
                       "org.biojava.bio.gui.sequence.MultiLineRenderer",
                       "RENDERERS",
                       SequenceRenderContext.LAYOUT);

    protected List renderers = new ArrayList();
    private transient ChangeForwarder rendererForwarder = null;

    protected ChangeSupport getChangeSupport(ChangeType ct) {
        ChangeSupport cs = super.getChangeSupport(ct);
    
        if (rendererForwarder == null) {
            rendererForwarder = new SequenceRenderer.RendererForwarder(this, cs);
            for (Iterator i = renderers.iterator(); i.hasNext(); ) {
                SequenceRenderer sRend = (SequenceRenderer) i.next();
                if (sRend instanceof Changeable) {
                    Changeable c = (Changeable) sRend;
                    c.addChangeListener(rendererForwarder,
                                        SequenceRenderContext.REPAINT);
                }
            }
        }

        return cs;
    }

    /**
     * <code>addRenderer</code> adds a renderer as a new track.
     *
     * @param renderer a <code>SequenceRenderer</code> to add.
     *
     * @exception ChangeVetoException if the renderer cannot be added.
     */
    public void addRenderer(SequenceRenderer renderer)
        throws ChangeVetoException {
        if (hasListeners()) {
            ChangeEvent ce = new ChangeEvent(this, RENDERERS, renderer, null);
            ChangeSupport cs = getChangeSupport(RENDERERS);

            synchronized(cs) {
                cs.firePreChangeEvent(ce);
                _addRenderer(renderer);
                if (renderer instanceof Changeable) {
                    Changeable c = (Changeable) renderer;
                    c.addChangeListener(rendererForwarder,
                                        SequenceRenderContext.REPAINT);
                }
                cs.firePostChangeEvent(ce);
            }
        } else {
            _addRenderer(renderer);
        }
    }

    protected void _addRenderer(SequenceRenderer renderer) {
        renderers.add(renderer);
    }

    /**
     * <code>removeRenderer</code> removes a renderer.
     *
     * @param renderer a <code>SequenceRenderer</code> to remove.
     *
     * @exception ChangeVetoException if the renderer can not be
     * removed.
     */
    public void removeRenderer(SequenceRenderer renderer)
        throws ChangeVetoException {
        if (hasListeners()) {
            ChangeEvent ce = new ChangeEvent(this, RENDERERS, null, renderer);
            ChangeSupport cs = getChangeSupport(RENDERERS);

            synchronized(cs) {
                cs.firePreChangeEvent(ce);
                _removeRenderer(renderer);
                if (renderer instanceof Changeable) {
                    Changeable c = (Changeable) renderer;
                    c.removeChangeListener(rendererForwarder,
                                           SequenceRenderContext.REPAINT);
                }
                cs.firePostChangeEvent(ce);
            }
        } else {
            _removeRenderer(renderer);
        }
    }

    protected void _removeRenderer(SequenceRenderer renderer) {
        renderers.remove(renderer);
    }

    /**
     * <code>clearRenderers</code> removes all renderers from this
     * renderer.
     *
     * @exception ChangeVetoException if the renderers can not be
     * cleared.
     */
    public void clearRenderers()
        throws ChangeVetoException {
        if (hasListeners()) {
            ChangeEvent ce = new ChangeEvent(this, RENDERERS);
            ChangeSupport cs = getChangeSupport(RENDERERS);

            synchronized(cs) {
                cs.firePreChangeEvent(ce);
                for (Iterator i = renderers.iterator(); i.hasNext(); ) {
                    Object r = i.next();
                    if (r instanceof Changeable) {
                        Changeable c = (Changeable) r;
                        c.removeChangeListener(rendererForwarder,
                                               SequenceRenderContext.REPAINT);
                    }
                }
                renderers.clear();
                cs.firePostChangeEvent(ce);
            }
        } else {
            renderers.clear();
        }
    }
  
    public double getDepth(SequenceRenderContext src) {
        return LayeredRenderer.INSTANCE.getDepth(Collections.nCopies(renderers.size(), src),
                                                 renderers);
    }

    public double getMinimumLeader(SequenceRenderContext src) {
        return LayeredRenderer.INSTANCE.getMinimumLeader(Collections.nCopies(renderers.size(), src),
                                                         renderers);
    }

    public double getMinimumTrailer(SequenceRenderContext src) {
        return LayeredRenderer.INSTANCE.getMinimumTrailer(Collections.nCopies(renderers.size(), src),
                                                          renderers);
    }
  
    public void paint(Graphics2D g, SequenceRenderContext src) {
        LayeredRenderer.INSTANCE.paint(g,
                                       Collections.nCopies(renderers.size(), src),
                                       renderers);
    }
  
    public SequenceViewerEvent processMouseEvent(SequenceRenderContext src,
                                                 MouseEvent me,
                                                 List path) {
        path.add(this);
        SequenceViewerEvent sve =
            LayeredRenderer.INSTANCE.processMouseEvent(Collections.nCopies(renderers.size(), src),
                                                       me,
                                                       path,
                                                       renderers);

        if (sve == null) {
            sve = new SequenceViewerEvent(this,
                                          null,
                                          src.graphicsToSequence(me.getPoint()),
                                          me,
                                          path);
        }

        return sve;
    }
}
