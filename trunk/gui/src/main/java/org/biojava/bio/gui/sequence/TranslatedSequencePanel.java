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
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Insets;
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
import java.util.ArrayList;

import javax.swing.JComponent;

import org.biojava.bio.seq.FeatureHolder;
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
 * <p><code>TranslatedSequencePanel</code> is a panel that displays a
 * Sequence. Its features are that it will always draw at low pixel
 * coordinates when using Java2D to render very long sequences and it
 * is quite fast (approximately 8x faster than
 * <code>SequencePanel</code></p>.
 *
 * <p>A <code>TranslatedSequencePanel</code> can either display the
 * sequence from left-to-right (HORIZONTAL) or from top-to-bottom
 * (VERTICAL). It has an associated scale which is the number of
 * pixels per symbol and a translation which is the number of
 * <code>Symbol</code>s to skip before rendering starts. In order to
 * produce a scrolling effect, the <code>setSymbolTranslation</code>
 * method may be hooked up to an <code>Adjustable</code> such as
 * <code>JScrollBar</code> or to an event listener.</p>
 *
 * <p>The exact number of <code>Symbol</code>s rendered depends on the
 * width of the panel and the scale. Resizing the panel will cause the
 * number of <code>Symbol</code>s rendered to change accordingly.</p>
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
 * @author Thomas Down
 * @author Jolyon Holdstock
 * @since 1.2
 */
public class TranslatedSequencePanel extends JComponent
    implements SequenceRenderContext, Changeable
{
    /**
	 * Generated Serial Version UID
	 */
	private static final long serialVersionUID = 3269477379497205817L;

	/**
     * Constant <code>RENDERER</code> is a <code>ChangeType</code>
     * which indicates a change to the renderer, requiring a layout
     * update.
     */
    public static final ChangeType RENDERER =
        new ChangeType("The renderer for this TranslatedSequencePanel has changed",
                       "org.biojava.bio.gui.sequence.TranslatedSequencePanel",
                       "RENDERER", SequenceRenderContext.LAYOUT);

    /**
     * Constant <code>TRANSLATION</code> is a <code>ChangeType</code>
     * which indicates a change to the translation, requiring a paint
     * update.
     */
    public static final ChangeType TRANSLATION =
        new ChangeType("The translation for this TranslatedSequencePanel has changed",
                       "org.biojava.bio.gui.sequence.TranslatedSequencePanel",
                       "TRANSLATION", SequenceRenderContext.REPAINT);

    // The sequence to be rendered
    private SymbolList sequence;
    // The number of residues to skip before starting to render
    private int translation;
    // The rendering direction (HORIZONTAL or VERTICAL)
    private int direction;
    // The rendering scale in pixels per residue
    private double scale;
    // The sequence renderer
    private SequenceRenderer renderer;
    // The total border size of the renderer. Cached to avoid
    // recursive method calls between us and the renderer.
    private double rendererBorders;
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
                renderer.processMouseEvent(TranslatedSequencePanel.this, me,
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
                renderer.processMouseEvent(TranslatedSequencePanel.this, me,
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
                renderer.processMouseEvent(TranslatedSequencePanel.this, me,
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
                renderer.processMouseEvent(TranslatedSequencePanel.this, me,
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
                renderer.processMouseEvent(TranslatedSequencePanel.this, me,
                                           new ArrayList());

            me.translatePoint(insets.left, insets.top);
            svmSupport.fireMouseMoved(sve);
        }
    };

    /**
     * Creates a new <code>TranslatedSequencePanel</code> with the
     * default settings (direction HORIZONTAL, scale 10.0 pixels per
     * symbol, symbol translation 0, leading border 0.0, trailing
     * border 0.0, 12 point sanserif font).
     */
    public TranslatedSequencePanel()
    {
        super();

        if (getFont() == null)
            setFont(new Font("sanserif", Font.PLAIN, 12));

        direction   = HORIZONTAL;
        scale       = 10.0;
        translation = 0;

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
    public SymbolList getSequence()
    {
        return sequence;
    }

    /**
     * <code>setSequence</code> sets the <code>Sequence</code> to be
     * rendered.
     *
     * @param sequence a <code>Sequence</code>.
     */
    public void setSequence(SymbolList sequence)
    {
        SymbolList prevSequence = this.sequence;

        // Remove out listener from the sequence, if necessary
        if (prevSequence != null)
            prevSequence.removeChangeListener(layoutListener);

        this.sequence = sequence;

        // Add our listener to the sequence, if necessary. Also update
        // rendererBorders cache which can be affected by changes in
        // sequence
        if (sequence != null)
        {
            if (renderer != null)
                rendererBorders = renderer.getMinimumLeader(this)
                + renderer.getMinimumTrailer(this);

            sequence.addChangeListener(layoutListener);
        }

        firePropertyChange("sequence", prevSequence, sequence);
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
     * <code>getFeatures</code> returns all of the
     * <code>Feature</code>s belonging to the currently rendered
     * <code>Sequence</code>.
     *
     * @return a <code>FeatureHolder</code>.
     */
    public FeatureHolder getFeatures()
    {
      if(sequence instanceof FeatureHolder) {
        return (FeatureHolder) sequence;
      } else {
        return FeatureHolder.EMPTY_FEATURE_HOLDER;
      }
    }

    /**
     * <code>getRange</code> returns a <code>RangeLocation</code>
     * representing the region of the sequence currently being
     * rendered. This is calculated from the size of the
     * <code>TranslatedSequencePanel</code>, minus its
     * <code>SequenceRenderContext.Border</code>s and its delegate
     * renderer borders (if any), the current rendering translation
     * and the current scale. The value will therefore change when the
     * <code>TranslatedSequencePanel</code> is resized or "scrolled"
     * by changing the translation.
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
     * <code>getDirection</code> returns the direction in which this
     * context expects sequences to be rendered - HORIZONTAL or
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
     * context will render sequences - HORIZONTAL or VERTICAL.
     *
     * @param direction an <code>int</code>.
     *
     * @exception IllegalArgumentException if an error occurs.
     */
    public void setDirection(int direction) throws IllegalArgumentException
    {
        if (direction != HORIZONTAL && direction != VERTICAL)
            throw new IllegalArgumentException("Direction must be either HORIZONTAL or VERTICAL");

        int prevDirection = this.direction;
        this.direction = direction;

        firePropertyChange("direction", prevDirection, direction);
        resizeAndValidate();
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
     * rendering. The sequence will be rendered, immediately after any
     * borders, starting at this translation. Values may be from 0 to
     * the length of the rendered sequence.
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
     * sequence will be rendered, immediately after any borders,
     * starting at that translation. Values may be from 0 to the
     * length of the rendered sequence.
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
     * <code>getLeadingBorder</code> returns the leading border.
     *
     * @return a <code>SequenceRenderContext.Border</code>.
     */
    public SequenceRenderContext.Border getLeadingBorder()
    {
        return leadingBorder;
    }

    /**
     * <code>getTrailingBorder</code> returns the trailing border.
     *
     * @return a <code>SequenceRenderContext.Border</code>.
     */
    public SequenceRenderContext.Border getTrailingBorder()
    {
        return trailingBorder;
    }

    /**
     * <code>getRenderer</code> returns the current
     * <code>SequenceRenderer</code>.
     *
     * @return a <code>SequenceRenderer</code>.
     */
    public SequenceRenderer getRenderer()
    {
        return renderer;
    }

    /**
     * <code>setRenderer</code> sets the current
     * <code>SequenceRenderer</code>.
     *
     * @param renderer  set the <code>SequenceRenderer</code> used
     */
    public void setRenderer(SequenceRenderer renderer)
        throws ChangeVetoException
    {
        if (! isActive())
            _setRenderer(renderer);

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
     * <code>sequenceToGraphics</code> converts a sequence index to a
     * graphical position.
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
     * @param point the <code>Point2D</code> to transform
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
     * <code>getVisibleSymbolCount</code> returns the
     * <strong>maximum</strong> number of <code>Symbol</code>s which
     * can be rendered in the visible area (excluding all borders) of
     * the <code>TranslatedSequencePanel</code> at the current
     * scale. Note that if the translation is greater than 0, the
     * actual number of <code>Symbol</code>s rendered will be less.
     *
     * @return an <code>int</code>.
     */
    public int getVisibleSymbolCount()
    {
        // The BioJava borders
        double totalBorders = leadingBorder.getSize()
            + trailingBorder.getSize() + rendererBorders;

        // The Insets
        Insets insets = getInsets();

        int visible;

        if (direction == HORIZONTAL)
        {
            int width = getWidth() - insets.left - insets.right;

            if (width <= totalBorders) 
                return 0;
            else
                visible = width - (int) totalBorders;
        }
        else
        {
            int height = getHeight() - insets.top - insets.bottom;

            if (height <= totalBorders)
                return 0;
            else
                visible = height - (int) totalBorders;
        }

        return Math.min(graphicsToSequence(visible), sequence.length());
    }

    /**
     * <code>paintComponent</code> paints this component.
     *
     * @param g a <code>Graphics</code> object.
     */
    public void paintComponent(Graphics g)
    {
        if (! isActive())
            return;

        super.paintComponent(g);

        // Set hints
        Graphics2D g2 = (Graphics2D) g;
        g2.addRenderingHints(hints);

        // Save current transform and clip
        AffineTransform prevTransform = g2.getTransform();
        Shape                prevClip = g2.getClip();
        Insets                 insets = getInsets();

        // As we subclass JComponent we have to paint our own
        // background, but only if we are opaque
        if (isOpaque())
        {
            g2.setPaint(getBackground());
            g2.fillRect(0, 0, getWidth(), getHeight());
        }

         Rectangle2D.Double clip = new Rectangle2D.Double();
         if (direction == HORIZONTAL) {
        	  // Clip x to edge of delegate renderer's leader
        	  //clip.x = renderer.getMinimumLeader(this);
        	  clip.x = 0 - renderer.getMinimumLeader(this);
        	  clip.y = 0.0;
        	  // Set the width to visible symbols + the delegate
        	  // renderer's minimum trailer (which may have something in
        	  // it to render).
        	  clip.width = sequenceToGraphics(getVisibleSymbolCount() + 1) +
        	renderer.getMinimumLeader(this) + renderer.getMinimumTrailer(this);
        	  clip.height = renderer.getDepth(this);
        	  g2.translate(leadingBorder.getSize() - clip.x + insets.left, insets.top); 
         }
         else
         {
             clip.x = 0.0;
             // Clip y to edge of delegate renderer's leader
             clip.y = renderer.getMinimumLeader(this);
             clip.width = renderer.getDepth(this);
             // Set the height to visible symbols + the delegate
             // renderer's minimum trailer (which may have something in
             // it to render).
             clip.height = sequenceToGraphics(getVisibleSymbolCount() + 1)
                 + renderer.getMinimumTrailer(this);

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
     * maximum sizes of the component according to the current leading
     * and trailing borders, renderer depth and visible symbol count.
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
            double width = sequenceToGraphics(getVisibleSymbolCount())
                + rendererBorders;
            double depth = renderer.getDepth(this);

            if (direction == HORIZONTAL)
            {
                d = new Dimension((int) Math.ceil(width),
                                  (int) Math.ceil(depth));
            }
            else
            {
                d = new Dimension((int) Math.ceil(depth),
                                  (int) Math.ceil(width));
            }
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
        if (hasListeners())
        {
            ChangeSupport cs = getChangeSupport(ct);
            cs.removeChangeListener(cl, ct);
        }
    }
    
    public boolean isUnchanging(ChangeType ct)
    {
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
        if (changeSupport != null)
        {
            return changeSupport;
        }
      
        synchronized(this)
        {
            if (changeSupport == null)
            {
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
     * <code>Sequence</code> to be rendered and the
     * <code>SequenceRenderer</code> are not null.
     *
     * @return a <code>boolean</code> value.
     */
    protected boolean isActive()
    {
        return (sequence != null) && (renderer != null);
    }

    /**
     * <code>_setRenderer</code> handles the details of listeners
     * during changes.
     *
     * @param renderer a <code>SequenceRenderer</code>.
     */
    private void _setRenderer(SequenceRenderer renderer)
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

        // Update our cache of the renderer's total border size, but
        // only is the sequence is not null. If the sequence was null
        // this value is no longer correct, so we have to update
        // rendererBorders in setSequence too.
        if (sequence != null)
            rendererBorders = renderer.getMinimumLeader(this)
                + renderer.getMinimumTrailer(this);

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

