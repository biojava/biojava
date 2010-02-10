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

import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Insets;
import java.awt.Point;
import java.awt.RenderingHints;
import java.awt.Shape;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.geom.AffineTransform;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.Serializable;
import java.util.ArrayList;

import javax.swing.JComponent;

import org.biojava.bio.seq.FeatureHolder;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.symbol.RangeLocation;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.utils.ChangeAdapter;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeListener;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.Changeable;

/**
 * <p>A <code>PairwiseSequencePanel</code> is a panel that displays a
 * pair of sequences; one sequence (the primary) may be either
 * left-to-right (HORIZONTAL) or from top-to-bottom (VERTICAL). The
 * other sequence (secondary) will then occupy the remaining
 * direction. It has an associated scale which is the number of pixels
 * per symbol and applies to both sequences. The leading and trailing
 * borders apply to the primary sequence only.</p>
 *
 * <p>The primary purpose of this component is to provide a means for
 * representing graphically comparisons between two sequences. This
 * could be anything from traditional dotplots (or variants created
 * with lines) to a more complex layered plot involving superimposed
 * renderers.</p>
 *
 * <p>Each sequence has a translation which is the number of
 * <code>Symbol</code>s to skip before rendering starts. In order to
 * produce a scrolling effect, the <code>setSymbolTranslation</code>
 * or <code>setSecondarySymbolTranslation</code> method may be hooked
 * up to an <code>Adjustable</code> such as <code>JScrollBar</code> or
 * to an event listener.</p>
 *
 * <p>The exact number of <code>Symbol</code>s rendered in each
 * sequence depends on the dimensions of the panel and the
 * scale. Resizing the panel will cause the number of
 * <code>Symbol</code>s rendered to change accordingly.</p>
 *
 * <p>The panel will fill its background to the <code>Color</code>
 * defined by the <code>setBackground()</code> method provided that it
 * has been defined as opaque using <code>setOpaque()</code>.</p>
 *
 * <p>The change event handling code is based on the original panel
 * and other BioJava components by Matthew and Thomas.</p>
 *
 * @author Keith James
 * @author Matthew Pocock
 * @since 1.2
 */
public class PairwiseSequencePanel extends JComponent
    implements PairwiseRenderContext, Changeable, Serializable
{
    /**
     * Constant <code>RENDERER</code> is a <code>ChangeType</code>
     * which indicates a change to the renderer, requiring a layout
     * update.
     */
    public static final ChangeType RENDERER =
        new ChangeType("The renderer for this PairwiseSequencePanel has changed",
                       "org.biojava.bio.gui.sequence.PairwiseSequencePanel",
                       "RENDERER", SequenceRenderContext.LAYOUT);

    /**
     * Constant <code>TRANSLATION</code> is a <code>ChangeType</code>
     * which indicates a change to the translation, requiring a paint
     * update.
     */
    public static final ChangeType TRANSLATION =
        new ChangeType("The translation for this PairwiseSequencePanel has changed",
                       "org.biojava.bio.gui.sequence.PairwiseSequencePanel",
                       "TRANSLATION", SequenceRenderContext.REPAINT);

    // The query sequence to be rendered
    private Sequence sequence;
    // The number of residues to skip before starting to render
    private int translation;
    // The rendering direction (HORIZONTAL or VERTICAL)
    private int direction;

    // The subject sequence to be rendered
    private Sequence secSequence;
    // The number of residues to skip before starting to render
    private int secTranslation;
    // The rendering direction (HORIZONTAL or VERTICAL)
    private int secDirection;

    // The rendering scale in pixels per residue
    private double scale;
    // The homology renderer
    private PairwiseSequenceRenderer renderer;

    // RenderingHints to be used by renderers
    private RenderingHints hints;

    // The rendering context borders
    private SequenceRenderContext.Border leadingBorder;
    private SequenceRenderContext.Border trailingBorder;

    // Listens for bound property changes which require a repaint
    // afterwards
    private PropertyChangeListener propertyListener =
    new PropertyChangeListener()
    {
        public void propertyChange(PropertyChangeEvent pce)
        {
            repaint();
        }
    };

    // ChangeSupport helper for BioJava ChangeListeners
    private transient ChangeSupport changeSupport = null;

    // Listens for BioJava changes which will require a repaint
    // afterwards
    private ChangeListener repaintListener = new ChangeAdapter()
    {
        public void postChange(ChangeEvent ce)
        {
            repaint();
        }
    };

    // Listens for BioJava changes which will require a revalidation
    // afterwards
    private ChangeListener layoutListener = new ChangeAdapter()
    {
        public void postChange(ChangeEvent ce)
        {
            revalidate();
        }
    };

    // SequenceViewerSupport helper for BioJava SequenceViewerListeners
    private SequenceViewerSupport svSupport = new SequenceViewerSupport();

    // Listens for mouse click events
    private MouseListener mouseListener = new MouseAdapter()
    {
        public void mouseClicked(MouseEvent me)
        {
            if (! isActive())
                return;

            Insets insets = getInsets();
            me.translatePoint(-insets.left, -insets.top);

            SequenceViewerEvent sve =
                renderer.processMouseEvent(PairwiseSequencePanel.this, me,
                                           new ArrayList());

            me.translatePoint(insets.left, insets.top);
            svSupport.fireMouseClicked(sve);
        }

        public void mousePressed(MouseEvent me)
        {
            if (! isActive())
                return;

            Insets insets = getInsets();
            me.translatePoint(-insets.left, -insets.top);

            SequenceViewerEvent sve =
                renderer.processMouseEvent(PairwiseSequencePanel.this, me,
                                           new ArrayList());

            me.translatePoint(insets.left, insets.top);
            svSupport.fireMousePressed(sve);
        }

        public void mouseReleased(MouseEvent me)
        {
            if (! isActive())
                return;

            Insets insets = getInsets();
            me.translatePoint(-insets.left, -insets.top);

            SequenceViewerEvent sve =
                renderer.processMouseEvent(PairwiseSequencePanel.this, me,
                                           new ArrayList());

            me.translatePoint(insets.left, insets.top);
            svSupport.fireMouseReleased(sve);
        }
    };

    // SequenceViewerMotionSupport helper for BioJava
    // SequenceViewerMotionListeners
    private SequenceViewerMotionSupport svmSupport =
        new SequenceViewerMotionSupport();

    // Listens for mouse movement events
    private MouseMotionListener mouseMotionListener = new MouseMotionListener()
    {
        public void mouseDragged(MouseEvent me)
        {
            if (! isActive())
                return;

            Insets insets = getInsets();
            me.translatePoint(-insets.left, -insets.top);

            SequenceViewerEvent sve =
                renderer.processMouseEvent(PairwiseSequencePanel.this, me,
                                           new ArrayList());

            me.translatePoint(insets.left, insets.top);
            svmSupport.fireMouseDragged(sve);
        }

        public void mouseMoved(MouseEvent me)
        {
            if (! isActive())
                return;

            Insets insets = getInsets();
            me.translatePoint(-insets.left, -insets.top);

            SequenceViewerEvent sve =
                renderer.processMouseEvent(PairwiseSequencePanel.this, me,
                                           new ArrayList());

            me.translatePoint(insets.left, insets.top);
            svmSupport.fireMouseMoved(sve);
        }
    };

    /**
     * Creates a new <code>PairwiseSequencePanel</code> with the
     * default settings (primary sequence direction HORIZONTAL, scale
     * 10.0 pixels per symbol, symbol translation 0, secondary symbol
     * translation 0, leading border 0.0, trailing border 0.0, 12
     * point sanserif font).
     */
    public PairwiseSequencePanel()
    {
        super();

        // Direction of the primary sequence
        direction      = HORIZONTAL;
        scale          = 10.0;
        translation    = 0;
        secTranslation = 0;

        leadingBorder  = new SequenceRenderContext.Border();
        trailingBorder = new SequenceRenderContext.Border();

        leadingBorder.setSize(0.0);
        trailingBorder.setSize(0.0);

        hints = new RenderingHints(null);

        this.addPropertyChangeListener(propertyListener);
        this.addMouseListener(mouseListener);
        this.addMouseMotionListener(mouseMotionListener);
    }

    /**
     * <code>getSequence</code> returns the entire
     * <code>Sequence</code> currently being rendered.
     *
     * @return a <code>Sequence</code>.
     */
    public Sequence getSequence()
    {
        return sequence;
    }

    /**
     * <code>setSequence</code> sets the <code>Sequence</code>
     * to be rendered.
     *
     * @param sequence a <code>Sequence</code>.
     */
    public void setSequence(Sequence sequence)
    {
        Sequence prevSequence = sequence;

        // Remove out listener from the sequence, if necessary
        if (prevSequence != null)
            prevSequence.removeChangeListener(layoutListener);

        this.sequence = sequence;

        // Add our listener to the sequence, if necessary
        if (sequence != null)
            sequence.addChangeListener(layoutListener);

        firePropertyChange("sequence", prevSequence, sequence);
        resizeAndValidate();
    }

    /**
     * <code>getSecondarySequence</code> returns the entire secondary
     * <code>Sequence</code> currently being rendered.
     *
     * @return a <code>Sequence</code>.
     */
    public Sequence getSecondarySequence()
    {
        return secSequence;
    }

    /**
     * <code>setSecondarySequence</code> sets the secondary
     * <code>Sequence</code> to be rendered.
     *
     * @param sequence a <code>Sequence</code>.
     */
    public void setSecondarySequence(Sequence sequence)
    {
        Sequence prevSecSequence = secSequence;

        // Remove out listener from the sequence, if necessary
        if (prevSecSequence != null)
            prevSecSequence.removeChangeListener(layoutListener);

        secSequence = sequence;

        // Add our listener to the sequence, if necessary
        if (sequence != null)
            sequence.addChangeListener(layoutListener);

        firePropertyChange("secSequence", prevSecSequence, sequence);
        resizeAndValidate();
    }

    /**
     * <code>getSymbols</code> returns all of the <code>Symbol</code>s
     * belonging to the currently rendered <code>Sequence</code>.
     *
     * @return a <code>SymbolList</code>.
     */
    public SymbolList getSymbols()
    {
        return sequence;
    }

    /**
     * <code>getSecondarySymbols</code> returns all of the
     * <code>Symbol</code>s belonging to the currently rendered
     * secondary <code>Sequence</code>.
     *
     * @return a <code>SymbolList</code>.
     */
    public SymbolList getSecondarySymbols()
    {
        return secSequence;
    }

    /**
     * <code>getFeatures</code> returns all of the
     * <code>Feature</code>s belonging to the currently rendered
     * <code>Sequence</code>.
     *
     * @return a <code>FeatureHolder</code>.
     */
    public FeatureHolder getFeatures()
    {
        return sequence;
    }

    /**
     * <code>getSecondaryFeatures</code> returns all of the
     * <code>Feature</code>s belonging to the currently rendered
     * secondary <code>Sequence</code>.
     *
     * @return a <code>FeatureHolder</code>.
     */
    public FeatureHolder getSecondaryFeatures()
    {
        return secSequence;
    }

    /**
     * <code>getRange</code> returns a <code>RangeLocation</code>
     * representing the region of the sequence currently being
     * rendered. This is calculated from the size of the
     * <code>PairwiseSequencePanel</code>, the current rendering
     * translation and the current scale. The value will therefore
     * change when the <code>PairwiseSequencePanel</code> is resized
     * or "scrolled" by changing the translation.
     *
     * @return a <code>RangeLocation</code>.
     */
    public RangeLocation getRange()
    {
        int visibleSymbols = getVisibleSymbolCount();

        // This is a fudge as we have to return a RangeLocation, which
        // can not have start == end
        if (visibleSymbols == 0)
            return new RangeLocation(translation + 1,
                                     translation + 2);
        else
            return new RangeLocation(translation + 1, visibleSymbols);
    }

    /**
     * <code>getSecondaryRange</code> returns a
     * <code>RangeLocation</code> representing the region of the
     * secondary sequence currently being rendered. This is calculated
     * from the size of the <code>PairwiseSequencePanel</code>, the
     * current rendering translation and the current scale. The value
     * will therefore change when the
     * <code>PairwiseSequencePanel</code> is resized or "scrolled" by
     * changing the translation.
     *
     * @return a <code>RangeLocation</code>.
     */
    public RangeLocation getSecondaryRange()
    {
        int visibleSecSymbols = getVisibleSecondarySymbolCount();

        // This is a fudge as we have to return a RangeLocation, which
        // can not have start == end
        if (visibleSecSymbols == 0)
            return new RangeLocation(secTranslation + 1,
                                     secTranslation + 2);
        else
            return new RangeLocation(secTranslation + 1, visibleSecSymbols);
    }

    /**
     * <code>getDirection</code> returns the direction in which this
     * context expects the sequence to be rendered - HORIZONTAL or
     * VERTICAL.
     *
     * @return an <code>int</code>.
     */
    public int getDirection()
    {
        return direction;
    }

    /**
     * <code>setDirection</code> sets the direction in which this
     * context will render the sequence - HORIZONTAL or VERTICAL.
     *
     * @param direction an <code>int</code>.
     *
     * @exception IllegalArgumentException if an invalid direction is
     * used.
     */
    public void setDirection(int direction)
        throws IllegalArgumentException
    {
        int    prevDirection = direction;
        int prevSecDirection = secDirection;

         if (direction == HORIZONTAL)
             secDirection = VERTICAL;
         else if (direction == VERTICAL)
             secDirection = HORIZONTAL;
         else
             throw new IllegalArgumentException("Direction must be either HORIZONTAL or VERTICAL");

         this.direction = direction;

         firePropertyChange("direction",    prevDirection,    direction);
         firePropertyChange("secDirection", prevSecDirection, secDirection);
         resizeAndValidate();
    }

    /**
     * <code>getSecondaryDirection</code> returns the direction in
     * which this context expects the secondary sequence to be
     * rendered - HORIZONTAL or VERTICAL.
     *
     * @return an <code>int</code>.
     */
    public int getSecondaryDirection()
    {
        return secDirection;
    }

    /**
     * <code>getScale</code> returns the scale in pixels per
     * <code>Symbol</code>.
     *
     * @return a <code>double</code>.
     */
    public double getScale()
    {
        return scale;
    }

    /**
     * <code>setScale</code> sets the scale in pixels per
     * <code>Symbol</code>.
     *
     * @param scale a <code>double</code>.
     */
    public void setScale(double scale)
    {
        double prevScale = this.scale;
        this.scale = scale;

        firePropertyChange("scale", prevScale, scale);
        resizeAndValidate();
    }

    /**
     * <code>getSymbolTranslation</code> returns the current
     * translation in <code>Symbol</code>s which will be applied when
     * rendering. The sequence will be rendered starting at this
     * translation. Values may be from 0 to the length of the rendered
     * sequence.
     *
     * @return an <code>int</code>.
     */
    public int getSymbolTranslation()
    {
        return translation;
    }

    /**
     * <code>setSymbolTranslation</code> sets the translation in
     * <code>Symbol</code>s which will be applied when rendering. The
     * sequence will be rendered starting at that translation. Values
     * may be from 0 to the length of the rendered sequence.
     *
     * @param translation an <code>int</code>.
     *
     * @exception IndexOutOfBoundsException if the translation is
     * greater than the sequence length.
     */
    public void setSymbolTranslation(int translation)
        throws IndexOutOfBoundsException
    {
        if (translation >= sequence.length())
            throw new IndexOutOfBoundsException("Tried to set symbol translation offset equal to or greater than SymbolList length");

        int prevTranslation = this.translation;

        if (hasListeners())
        {
            ChangeSupport cs = getChangeSupport(TRANSLATION);

            ChangeEvent ce = new ChangeEvent(this, TRANSLATION);

            cs.firePostChangeEvent(ce);
            this.translation = translation;
            cs.firePostChangeEvent(ce);
        }
        else
        {
            this.translation = translation;
        }

        firePropertyChange("translation", prevTranslation, translation);
        resizeAndValidate();
    }

    /**
     * <code>getSecondarySymbolTranslation</code> returns the current
     * translation in <code>Symbol</code>s which will be applied when
     * rendering. The secondary sequence will be rendered starting at
     * this translation. Values may be from 0 to the length of the
     * rendered sequence.
     *
     * @return an <code>int</code>.
     */
    public int getSecondarySymbolTranslation()
    {
        return secTranslation;
    }

    /**
     * <code>setSecondarySymbolTranslation</code> sets the translation
     * in <code>Symbol</code>s which will be applied when
     * rendering. The secondary sequence will be rendered starting at
     * that translation. Values may be from 0 to the length of the
     * rendered sequence.
     *
     * @param translation an <code>int</code>.
     *
     * @exception IndexOutOfBoundsException if the translation is
     * greater than the sequence length.
     */
    public void setSecondarySymbolTranslation(int translation)
        throws IndexOutOfBoundsException
    {
        if (translation >= secSequence.length())
            throw new IndexOutOfBoundsException("Tried to set secondary symbol translation offset equal to or greater than SymbolList length");

        int prevSecTranslation = secTranslation;

        if (hasListeners())
        {
            ChangeSupport cs = getChangeSupport(TRANSLATION);

            ChangeEvent ce = new ChangeEvent(this, TRANSLATION);

            cs.firePostChangeEvent(ce);
            secTranslation = translation;
            cs.firePostChangeEvent(ce);
        }
        else
        {
            secTranslation = translation;
        }

        firePropertyChange("secTranslation", prevSecTranslation, translation);
        resizeAndValidate();
    }

    /**
     * <code>getLeadingBorder</code> returns the leading border of the
     * primary sequence.
     *
     * @return a <code>SequenceRenderContext.Border</code>.
     */
    public SequenceRenderContext.Border getLeadingBorder()
    {
        return leadingBorder;
    }

    /**
     * <code>getTrailingBorder</code> returns the trailing border of
     * the primary sequence.
     *
     * @return a <code>SequenceRenderContext.Border</code>.
     */
    public SequenceRenderContext.Border getTrailingBorder()
    {
        return trailingBorder;
    }

    /**
     * <code>getRenderer</code> returns the current
     * <code>PairwiseSequenceRenderer</code>.
     *
     * @return a <code>PairwiseSequenceRenderer</code>.
     */
    public PairwiseSequenceRenderer getRenderer()
    {
        return renderer;
    }

    /**
     * <code>setRenderer</code> sets the current
     * <code>PairwiseSequenceRenderer</code>.
     */
    public void setRenderer(PairwiseSequenceRenderer renderer)
        throws ChangeVetoException
    {
        if (hasListeners())
        {
            ChangeSupport cs = getChangeSupport(RENDERER);

            // We are originating a change event so use this
            // constructor
            ChangeEvent ce = new ChangeEvent(this, RENDERER,
                                             renderer, this.renderer);

            synchronized(cs)
            {
                cs.firePreChangeEvent(ce);
                _setRenderer(renderer);
                cs.firePostChangeEvent(ce);
            }
        }
        else
        {
            _setRenderer(renderer);
        }
        resizeAndValidate();
    }

    /**
     * <code>getRenderingHints</code> returns the
     * <code>RenderingHints</code> currently being used by the
     * <code>Graphics2D</code> instances of delegate renderers. If
     * none is set, the constructor creates one with a null
     * <code>Map</code>.
     *
     * @return a <code>RenderingHints</code>.
     */
    public RenderingHints getRenderingHints()
    {
        return hints;
    }

    /**
     * <code>setRenderingHints</code> sets the
     * <code>RenderingHints</code> which will be used by the
     * <code>Graphics2D</code> instances of delegate renderers.
     *
     * @param hints a <code>RenderingHints</code>.
     */
    public void setRenderingHints(RenderingHints hints)
    {
        RenderingHints prevHints = this.hints;
        this.hints = hints;

        firePropertyChange("hints", prevHints, hints);
    }

    /**
     * <code>sequenceToGraphics</code> converts a sequence index
     * to a graphical position.
     *
     * @param sequencePos an <code>int</code>.
     *
     * @return a <code>double</code>.
     */
    public double sequenceToGraphics(int sequencePos)
    {
        return (sequencePos - translation - 1) * scale;
    }

    /**
     * <code>secondarySequenceToGraphics</code> converts a sequence
     * index to a graphical position.
     *
     * @param sequencePos an <code>int</code>.
     *
     * @return a <code>double</code>.
     */
    public double secondarySequenceToGraphics(int sequencePos)
    {
        return (sequencePos - secTranslation - 1) * scale;
    }

    /**
     * <code>graphicsToSequence</code> converts a graphical position
     * to a sequence index.
     *
     * @param graphicsPos a <code>double</code>.
     *
     * @return an <code>int</code>.
     */
    public int graphicsToSequence(double graphicsPos)
    {
        return ((int) (graphicsPos / scale)) + translation + 1;
    }

    /**
     * <code>graphicsToSequence</code> converts a graphical position
     * to a sequence index.
     *
     * @param point a graphic position.
     *
     * @return an <code>int</code>.
     */
    public int graphicsToSequence(Point2D point)
    {
        if (direction == HORIZONTAL)
            return graphicsToSequence(point.getX());
        else
            return graphicsToSequence(point.getY());
    }

    /**
     * <code>graphicsToSecondarySequence</code> converts a graphical
     * position to a secondary sequence index.
     *
     * @param graphicsPos a <code>double</code>.
     *
     * @return an <code>int</code>.
     */
    public int graphicsToSecondarySequence(double graphicsPos)
    {
        return ((int) (graphicsPos / scale)) + secTranslation + 1;
    }

    /**
     * <code>graphicsToSecondarySequence</code> converts a graphical
     * position to a secondary sequence index.
     *
     * @param point a <code>Point</code>.
     *
     * @return an <code>int</code>.
     */
    public int graphicsToSecondarySequence(Point point)
    {
        if (secDirection == HORIZONTAL)
            return graphicsToSecondarySequence(point.getX());
        else
            return graphicsToSecondarySequence(point.getY());
    }

    /**
     * <code>getVisibleSymbolCount</code> returns the
     * <strong>maximum</strong> number of <code>Symbol</code>s which
     * can be rendered in the visible area (excluding all borders) of
     * the <code>PairwiseSequencePanel</code> at the current
     * scale. Note that if the translation is greater than 0, the
     * actual number of <code>Symbol</code>s rendered will be less.
     *
     * @return an <code>int</code>.
     */
    public int getVisibleSymbolCount()
    {
        // The Insets
        Insets insets = getInsets();

        int visible;

        if (direction == HORIZONTAL)
            visible = getWidth() - insets.left - insets.right;
        else
            visible = getHeight() - insets.top - insets.bottom;

        return Math.min(graphicsToSequence(visible), sequence.length());
    }

    /**
     * <code>getVisibleSecondarySymbolCount</code> returns the
     * <strong>maximum</strong> number of secondary
     * <code>Symbol</code>s which can be rendered in the visible area
     * (excluding all borders) of the
     * <code>PairwiseSequencePanel</code> at the current scale. Note
     * that if the translation is greater than 0, the actual number of
     * <code>Symbol</code>s rendered will be less.
     *
     * @return an <code>int</code>.
     */
    public int getVisibleSecondarySymbolCount()
    {
        // The Insets
        Insets insets = getInsets();

        int visible;

        if (secDirection == HORIZONTAL)
            visible = getWidth() - insets.left - insets.right;
        else
            visible = getHeight() - insets.top - insets.bottom;

        return Math.min(graphicsToSecondarySequence(visible),
                        secSequence.length());
    }

    public void paintComponent(Graphics g)
    {
        if (! isActive())
            return;

        super.paintComponent(g);

        // Set hints
        Graphics2D g2 = (Graphics2D) g;
        g2.addRenderingHints(hints);

        // As we subclass JComponent we have to paint our own
        // background, but only if we are opaque
        if (isOpaque())
        {
            g2.setPaint(getBackground());
            g2.fillRect(0, 0, getWidth(), getHeight());
        }

        // Save current transform and clip
        AffineTransform prevTransform = g2.getTransform();
        Shape                prevClip = g2.getClip();
        Insets                 insets = getInsets();

        Rectangle2D.Double clip = new Rectangle2D.Double();

        clip.x = 0.0;
        clip.y = 0.0;

        if (direction == HORIZONTAL)
        {
            clip.width  = sequenceToGraphics(getVisibleSymbolCount() + 1);
            clip.height = secondarySequenceToGraphics(getVisibleSecondarySymbolCount() + 1);
            g2.translate(leadingBorder.getSize() + insets.left, insets.top);
        }
        else
        {
            clip.width  = secondarySequenceToGraphics(getVisibleSecondarySymbolCount() + 1);
            clip.height = sequenceToGraphics(getVisibleSymbolCount() + 1);
            g2.translate(insets.left, leadingBorder.getSize() + insets.top);
        }

        // Clip and paint
        g2.clip(clip);
        renderer.paint(g2, this);

        // Restore
        g2.setTransform(prevTransform);
        g2.setClip(prevClip);
    }

    /**
     * <code>resizeAndValidate</code> sets the minimum, preferred and
     * maximum sizes of the component according to the current visible
     * symbol count.
     */
    public void resizeAndValidate()
    {
        Dimension d = null;

        if (! isActive())
        {
            d = new Dimension(0, 0);
        }
        else
        {
            double width;
            double height;

            if (direction == HORIZONTAL)
            {
                width  = sequenceToGraphics(getVisibleSymbolCount());
                height = secondarySequenceToGraphics(getVisibleSecondarySymbolCount());
            }
            else
            {
                width  = secondarySequenceToGraphics(getVisibleSecondarySymbolCount());
                height = sequenceToGraphics(getVisibleSymbolCount());
            }

            d = new Dimension((int) Math.ceil(width),
                              (int) Math.ceil(height));
        }

        setMinimumSize(d);
        setPreferredSize(d);
        setMaximumSize(d);
        revalidate();
    }

    /**
     * <code>addChangeListener</code> adds a listener for all types of
     * change.
     *
     * @param cl a <code>ChangeListener</code>.
     */
    public void addChangeListener(ChangeListener cl)
    {
        addChangeListener(cl, ChangeType.UNKNOWN);
    }

    /**
     * <code>addChangeListener</code> adds a listener for specific
     * types of change.
     *
     * @param cl a <code>ChangeListener</code>.
     * @param ct a <code>ChangeType</code>.
     */
    public void addChangeListener(ChangeListener cl, ChangeType ct)
    {
      ChangeSupport cs = getChangeSupport(ct);
      cs.addChangeListener(cl, ct);
    }

    /**
     * <code>removeChangeListener</code> removes a listener.
     *
     * @param cl a <code>ChangeListener</code>.
     */
    public void removeChangeListener(ChangeListener cl)
    {
        removeChangeListener(cl, ChangeType.UNKNOWN);
    }

    /**
     * <code>removeChangeListener</code> removes a listener.
     *
     * @param cl a <code>ChangeListener</code>.
     * @param ct a <code>ChangeType</code>.
     */
    public void removeChangeListener(ChangeListener cl, ChangeType ct)
    {
      if(hasListeners()) {
        ChangeSupport cs = getChangeSupport(ct);
        cs.removeChangeListener(cl);
      }
    }

    public boolean isUnchanging(ChangeType ct) {
      ChangeSupport cs = getChangeSupport(ct);
      return cs.isUnchanging(ct);
    }

    /**
     * <code>addSequenceViewerListener</code> adds a listener for
     * mouse click <code>SequenceViewerEvent</code>s.
     *
     * @param svl a <code>SequenceViewerListener</code>.
     */
    public void addSequenceViewerListener(SequenceViewerListener svl)
    {
        svSupport.addSequenceViewerListener(svl);
    }

    /**
     * <code>removeSequenceViewerListener</code> removes a listener
     * for mouse click <code>SequenceViewerEvent</code>s.
     *
     * @param svl a <code>SequenceViewerListener</code>.
     */
    public void removeSequenceViewerListener(SequenceViewerListener svl)
    {
        svSupport.removeSequenceViewerListener(svl);
    }

    /**
     * <code>addSequenceViewerMotionListener</code> adds a listener for
     * mouse motion <code>SequenceViewerEvent</code>s.
     *
     * @param svml a <code>SequenceViewerMotionListener</code>.
     */
    public void addSequenceViewerMotionListener(SequenceViewerMotionListener svml)
    {
        svmSupport.addSequenceViewerMotionListener(svml);
    }

    /**
     * <code>addSequenceViewerMotionListener</code> removes a listener for
     * mouse motion <code>SequenceViewerEvent</code>s.
     *
     * @param svml a <code>SequenceViewerMotionListener</code>.
     */
    public void removeSequenceViewerMotionListener(SequenceViewerMotionListener svml)
    {
        svmSupport.removeSequenceViewerMotionListener(svml);
    }

    /**
     * <code>getChangeSupport</code> lazily instantiates a helper for
     * change listeners.
     *
     * @param ct a <code>ChangeType</code>.
     *
     * @return a <code>ChangeSupport</code> object.
     */
    protected ChangeSupport getChangeSupport(ChangeType ct)
    {
      if(changeSupport != null) {
        return changeSupport;
      }

      synchronized(this) {
        if (changeSupport == null) {
          changeSupport = new ChangeSupport();
        }

        return changeSupport;
      }
    }

    /**
     * <code>hasListeners</code> returns true if there are active
     * listeners for BioJava events.
     *
     * @return a <code>boolean</code> value.
     */
    protected boolean hasListeners()
    {
        return changeSupport != null;
    }

    /**
     * <code>isActive</code> returns true if both the
     * <code>Sequence</code>s to be rendered and the
     * <code>PairwiseHomologyRenderer</code> are not null.
     *
     * @return a <code>boolean</code> value.
     */
    protected boolean isActive()
    {
        return (sequence != null) && (secSequence != null) &&
            (renderer != null);
    }

    private void _setRenderer(PairwiseSequenceRenderer renderer)
    {
        // Remove our listeners from the old renderer, if necessary
        if (this.renderer != null &&
            Changeable.class.isInstance(this.renderer))
        {
            Changeable c = (Changeable) this.renderer;

            c.removeChangeListener(layoutListener,  SequenceRenderContext.LAYOUT);
            c.removeChangeListener(repaintListener, SequenceRenderContext.REPAINT);
        }

        this.renderer = renderer;

        // Add our listeners to the new renderer, if necessary
        if (renderer != null &&
            Changeable.class.isInstance(renderer))
        {
            Changeable c = (Changeable) renderer;
            c.addChangeListener(layoutListener,  SequenceRenderContext.LAYOUT);
            c.addChangeListener(repaintListener, SequenceRenderContext.REPAINT);
        }
    }
}
