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
 * <code>PairwiseOverlayRenderer</code> allows a list of other
 * <code>PairwiseSequenceRenderer</code>s to superimpose their
 * output. Operations on the child renderers (rendering, event
 * handling) are carried out in the order in which they were
 * added.
 *
 * @author Keith James
 * @author Matthew Pocock
 * @since 1.2
 */
public class PairwiseOverlayRenderer extends AbstractChangeable
    implements PairwiseSequenceRenderer, Serializable
{
    /**
     * Constant <code>RENDERERS</code> indicating a change to the
     * renderers handled by the overlay.
     */
    public static final ChangeType RENDERERS =
        new ChangeType("A PairwiseSequenceRenderer has been added or removed",
                       "org.biojava.bio.gui.sequence.PairwiseOverlayRenderer",
                       "RENDERERS",
                       SequenceRenderContext.LAYOUT);

    private List renderers;
    private transient ChangeForwarder rendererForwarder = null;

    /**
     * Creates a new, empty <code>PairwiseOverlayRenderer</code>.
     */
    public PairwiseOverlayRenderer()
    {
        renderers = new ArrayList();
    }
  
    /**
     * <code>addRenderer</code> adds a renderer.
     *
     * @param renderer a <code>PairwiseSequenceRenderer</code>.
     *
     * @exception ChangeVetoException if an error occurs.
     */
    public void addRenderer(PairwiseSequenceRenderer renderer)
        throws ChangeVetoException
    {
        if (hasListeners())
        {
            ChangeSupport cs = getChangeSupport(RENDERERS);

            ChangeEvent ce = new ChangeEvent(this, RENDERERS, renderer, null);

            synchronized(cs)
            {
                cs.firePreChangeEvent(ce);
                _addRenderer(renderer);

                if (renderer instanceof Changeable)
                {
                    Changeable c = (Changeable) renderer;
                    c.addChangeListener(rendererForwarder,
                                        SequenceRenderContext.REPAINT);
                }
                cs.firePostChangeEvent(ce);
            }
        }
        else
        {
            _addRenderer(renderer);
        }
    }
  
    /**
     * <code>removeRenderer</code> removes a renderer.
     *
     * @param renderer a <code>PairwiseSequenceRenderer</code>.
     *
     * @exception ChangeVetoException if an error occurs.
     */
    public void removeRenderer(PairwiseSequenceRenderer renderer)
        throws ChangeVetoException
    {
        if (hasListeners())
        {
            ChangeSupport cs = getChangeSupport(RENDERERS);

            ChangeEvent ce = new ChangeEvent(this, RENDERERS, null, renderer);

            synchronized(cs)
            {
                cs.firePreChangeEvent(ce);
                _removeRenderer(renderer);

                if (renderer instanceof Changeable)
                {
                    Changeable c = (Changeable) renderer;
                    c.removeChangeListener(rendererForwarder,
                                           SequenceRenderContext.REPAINT);
                }
                cs.firePostChangeEvent(ce);
            }
        }
        else
        {
            _removeRenderer(renderer);
        }
    }

    /**
     * <code>clearRenderers</code> removes all the renderers.
     *
     * @exception ChangeVetoException if an error occurs.
     */
    public void clearRenderers() throws ChangeVetoException
    {
        if (hasListeners())
        {
            ChangeSupport cs = getChangeSupport(RENDERERS);

            ChangeEvent ce = new ChangeEvent(this, RENDERERS);

            synchronized(cs)
            {
                cs.firePreChangeEvent(ce);
                for (Iterator i = renderers.iterator(); i.hasNext();)
                {
                    Object renderer = i.next();
                    if (renderer instanceof Changeable)
                    {
                        Changeable c = (Changeable) renderer;
                        c.removeChangeListener(rendererForwarder,
                                               SequenceRenderContext.REPAINT);
                    }
                }
                renderers.clear();
                cs.firePostChangeEvent(ce);
            }
        }
        else
        {
            renderers.clear();
        }
    }

    /**
     * <code>paint</code> applies all renderers in the order in which
     * they were added.
     *
     * @param g2 a <code>Graphics2D</code>.
     * @param context a <code>PairwiseRenderContext</code>.
     */
    public void paint(Graphics2D g2, PairwiseRenderContext context)
    {
        for (Iterator ri = renderers.iterator(); ri.hasNext();)
        {
            PairwiseSequenceRenderer renderer = (PairwiseSequenceRenderer) ri.next();
            renderer.paint(g2, context);
        }
    }

    protected ChangeSupport getChangeSupport(ChangeType ct)
    {
        ChangeSupport cs = super.getChangeSupport(ct);
    
        if (rendererForwarder == null)
        {
            rendererForwarder =
                new PairwiseSequenceRenderer.PairwiseRendererForwarder(this, cs);

            for (Iterator i = renderers.iterator(); i.hasNext();)
            {
                PairwiseSequenceRenderer renderer = (PairwiseSequenceRenderer) i.next();
                if (renderer instanceof Changeable)
                {
                    Changeable c = (Changeable) renderer;
                    c.addChangeListener(rendererForwarder,
                                        SequenceRenderContext.REPAINT);
                }
            }
        }
    
        return cs;
    }

    /**
     * <code>processMouseEvent</code> produces a
     * <code>SequenceViewerEvent</code> in response to a mouse
     * gesture. The list of renderers are probed in the order in which
     * they were added and the first renderer to accept the event will
     * return.
     *
     * @param context a <code>PairwiseRenderContext</code>.
     * @param me a <code>MouseEvent</code> that caused the request.
     * @param path a <code>List</code> of
     * <code>PairwiseSequenceRenderer</code> instances passed through
     * so far.
     *
     * @return a <code>SequenceViewerEvent</code> encapsulating the
     * mouse gesture.
     */
    public SequenceViewerEvent processMouseEvent(PairwiseRenderContext context,
                                                 MouseEvent            me,
                                                 List                  path)
    {
        path.add(this);

        SequenceViewerEvent event = null;

        List contextCopies = Collections.nCopies(renderers.size(), context);
        Iterator ci = contextCopies.iterator();
        Iterator ri = renderers.iterator();

        while (ci.hasNext() && ri.hasNext())
        {
            PairwiseRenderContext contextCopy = (PairwiseRenderContext) ci.next();
            PairwiseSequenceRenderer renderer = (PairwiseSequenceRenderer) ri.next();
            event = renderer.processMouseEvent(contextCopy, me, path);
        }

        if (event == null)
            event = new SequenceViewerEvent(this, null,
                                            context.graphicsToSequence(me.getPoint()),
                                            me, path);

        return event;
    }

    protected void _addRenderer(PairwiseSequenceRenderer renderer)
    {
        renderers.add(renderer);
    }
 
    protected void _removeRenderer(PairwiseSequenceRenderer renderer)
    {
        renderers.remove(renderer);
    }
}
