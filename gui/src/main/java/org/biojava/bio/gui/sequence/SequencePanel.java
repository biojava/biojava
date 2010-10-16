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
import java.beans.PropertyChangeSupport;
import java.io.Serializable;
import java.util.ArrayList;

import javax.swing.JComponent;
import javax.swing.SwingConstants;

import org.biojava.bio.seq.FeatureHolder;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.LocationTools;
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
 * A panel that displays a Sequence.
 * <p>
 * A SequencePanel can either display the sequence from left-to-right
 * (HORIZONTAL) or from top-to-bottom (VERTICAL). It has an associated scale
 * which is the number of pixels per symbol. It also has a lines property that
 * controls how to wrap the sequence off one end and onto the other.
 * <p>
 * Each line in the SequencePanel is broken down into a list of strips,
 * each rendered by an individual SequenceRenderer object.
 * You could add a SequenceRenderer that draws on genes, another that
 * draws repeats and another that prints out the DNA sequence. They are
 * responsible for rendering their view of the sequence in the place that the
 * SequencePanel positions them.
 *
 * @author Thomas Down
 * @author Matthew Pocock
 * @author David Huen
 */
public class SequencePanel
  extends
    JComponent
  implements
    SwingConstants,
    SequenceRenderContext,
    Changeable
{
  public static final ChangeType RENDERER = new ChangeType(
    "The renderer for this SequencePanel has changed",
    "org.biojava.bio.gui.sequence.SequencePanel",
    "RENDERER",
    SequenceRenderContext.LAYOUT
  );

  private SymbolList sequence;
  private RangeLocation range;
  private int direction;
  private double scale;
  private double pixelOffset;

  private SequenceRenderContext.Border leadingBorder;
  private SequenceRenderContext.Border trailingBorder;

  private SequenceRenderer renderer;
  private RendererMonitor theMonitor;

  private RenderingHints hints = null;

  private transient ChangeSupport changeSupport = null;

  private SequenceViewerSupport svSupport = new SequenceViewerSupport();

  /**
   * Use this to switch on effects like Anti-aliasing etc
   * @param hints the desired rendering properties
   */
  public void setRenderingHints(RenderingHints hints){
    this.hints = hints;
  }

  /**
   * @return the current rendering properties
   */
  public RenderingHints getRenderingHints(){
    return hints;
  }

  private MouseListener mouseListener = new MouseAdapter() {
    public void mouseClicked(MouseEvent me) {
      if(!isActive()) {
        return;
      }

      int [] dist = calcDist();
      me.translatePoint(+dist[0], +dist[1]);
      SequenceViewerEvent sve = renderer.processMouseEvent(
        SequencePanel.this,
        me,
        new ArrayList()
      );
      me.translatePoint(-dist[0], -dist[1]);
      svSupport.fireMouseClicked(sve);
    }

    public void mousePressed(MouseEvent me) {
      if(!isActive()) {
        return;
      }

      int [] dist = calcDist();
      me.translatePoint(+dist[0], +dist[1]);
      SequenceViewerEvent sve = renderer.processMouseEvent(
        SequencePanel.this,
        me,
        new ArrayList()
      );
      me.translatePoint(-dist[0], -dist[1]);
      svSupport.fireMousePressed(sve);
    }

    public void mouseReleased(MouseEvent me) {
      if(!isActive()) {
        return;
      }

      int [] dist = calcDist();
      me.translatePoint(+dist[0], +dist[1]);
      SequenceViewerEvent sve = renderer.processMouseEvent(
        SequencePanel.this,
        me,
        new ArrayList()
      );
      me.translatePoint(-dist[0], -dist[1]);
      svSupport.fireMouseReleased(sve);
    }
  };
  public void addSequenceViewerListener(SequenceViewerListener svl) {
    svSupport.addSequenceViewerListener(svl);
  }
  public void removeSequenceViewerListener(SequenceViewerListener svl) {
    svSupport.removeSequenceViewerListener(svl);
  }

  private SequenceViewerMotionSupport svmSupport = new SequenceViewerMotionSupport();
  private MouseMotionListener mouseMotionListener = new MouseMotionListener() {
    public void mouseDragged(MouseEvent me) {
      if(!isActive()) {
        return;
      }

      int [] dist = calcDist();
      me.translatePoint(+dist[0], +dist[1]);
      SequenceViewerEvent sve = renderer.processMouseEvent(
        SequencePanel.this,
        me,
        new ArrayList()
      );
      me.translatePoint(-dist[0], -dist[1]);
      svmSupport.fireMouseDragged(sve);
    }

    public void mouseMoved(MouseEvent me) {
      if(!isActive()) {
        return;
      }

      int [] dist = calcDist();
      me.translatePoint(+dist[0], +dist[1]);
      SequenceViewerEvent sve = renderer.processMouseEvent(
        SequencePanel.this,
        me,
        new ArrayList()
      );
      me.translatePoint(-dist[0], -dist[1]);
      svmSupport.fireMouseMoved(sve);
    }
  };
  public void addSequenceViewerMotionListener(SequenceViewerMotionListener svml) {
    svmSupport.addSequenceViewerMotionListener(svml);
  }
  public void removeSequenceViewerMotionListener(SequenceViewerMotionListener svml) {
    svmSupport.removeSequenceViewerMotionListener(svml);
  }

  protected boolean hasChangeListeners() {
    return changeSupport != null;
  }

  protected ChangeSupport getChangeSupport(ChangeType ct) {
    if(changeSupport != null) {
      return changeSupport;
    }

    synchronized(this) {
      if(changeSupport == null) {
        changeSupport = new ChangeSupport();
      }

      return changeSupport;
    }
  }

  protected boolean hasListeners() {
    return changeSupport != null;
  }

  public void addChangeListener(ChangeListener cl) {
    addChangeListener(cl, ChangeType.UNKNOWN);
  }

  public void addChangeListener(ChangeListener cl, ChangeType ct) {
    ChangeSupport cs = getChangeSupport(ct);
    cs.addChangeListener(cl, ct);
  }

  public void removeChangeListener(ChangeListener cl) {
    removeChangeListener(cl, ChangeType.UNKNOWN);
  }

  public void removeChangeListener(ChangeListener cl, ChangeType ct) {
    if(hasListeners()) {
      ChangeSupport cs = getChangeSupport(ct);
      cs.removeChangeListener(cl, ct);
    }
  }

  public boolean isUnchanging(ChangeType ct) {
    ChangeSupport cs = getChangeSupport(ct);
    return cs.isUnchanging(ct);
  }

  private ChangeListener layoutListener = new ChangeAdapter() {
    public void postChange(ChangeEvent ce) {
        System.err.println("Layout event");
      resizeAndValidate();
    }
  };
  private ChangeListener repaintListener = new ChangeAdapter() {
    public void postChange(ChangeEvent ce) {
        System.err.println("Repaint event for " + hashCode());
      repaint();
    }
  };

  /**
   * Initializer.
   */

  {
    direction = HORIZONTAL;
    scale = 12.0;
    pixelOffset = 0.0;

    theMonitor = new RendererMonitor();
    leadingBorder = new SequenceRenderContext.Border();
    trailingBorder = new SequenceRenderContext.Border();
  }

  /**
   * Create a new SequencePanel.
   */
  public SequencePanel() {
    super();
    if(getFont() == null) {
      setFont(new Font("serif", Font.PLAIN, 12));
    }
    this.addPropertyChangeListener(theMonitor);
    this.addMouseListener(mouseListener);
    this.addMouseMotionListener(mouseMotionListener);
  }

  /**
   * Set the SymboList to be rendered. This symbol list will be passed onto the
   * SequenceRenderer instances registered with this SequencePanel.
   *
   * @param s  the SymboList to render
   */
  public void setSequence(SymbolList s) {
    SymbolList oldSequence = sequence;
    if(oldSequence != null) {
      oldSequence.removeChangeListener(layoutListener, ChangeType.UNKNOWN);
    }
    this.sequence = s;
    if(s != null) {
      sequence.addChangeListener(layoutListener, ChangeType.UNKNOWN);
    }

    resizeAndValidate();
    firePropertyChange("sequence", oldSequence, s);
  }

  public SymbolList getSequence() {
    return sequence;
  }

  /**
   * Retrieve the currently rendered SymbolList
   *
   * @return  the current SymbolList
   */
  public SymbolList getSymbols() {
    return sequence;
  }

  public FeatureHolder getFeatures() {
    if(sequence instanceof FeatureHolder) {
      return (FeatureHolder) sequence;
    } else {
      return FeatureHolder.EMPTY_FEATURE_HOLDER;
    }
  }

  public void setRange(RangeLocation range) {
    RangeLocation oldRange = this.range;
    this.range = range;
    resizeAndValidate();
    firePropertyChange("range", oldRange, range);
  }

  public RangeLocation getRange() {
    return this.range;
  }

  /**
   * Set the direction that this SequencePanel renders in. The direction can be
   * one of HORIZONTAL or VERTICAL. Once the direction is set, the display will
   * redraw. HORIZONTAL represents left-to-right rendering. VERTICAL represents
   * AceDB-style vertical rendering.
   *
   * @param dir  the new rendering direction
   */
  public void setDirection(int dir)
  throws IllegalArgumentException {
    if(dir != HORIZONTAL && dir != VERTICAL) {
      throw new IllegalArgumentException(
        "Direction must be either HORIZONTAL or VERTICAL"
      );
    }
    int oldDirection = direction;
    direction = dir;
    resizeAndValidate();
    firePropertyChange("direction", oldDirection, direction);
  }

  /**
   * Retrieve the current rendering direction.
   *
   * @return the rendering direction (one of HORIZONTAL and VERTICAL)
   */
  public int getDirection() {
    return direction;
  }

  /**
   * Set the scale.
   * <p>
   * The scale parameter is interpreted as the number of pixels per symbol. This
   * may take on a wide range of values - for example, to render the symbols as
   * text, you will need a scale of > 8, where as to render chromosome 1 you
   * will want a scale &lt; 0.00000001
   *
   * @param scale the new pixels-per-symbol ratio
   */
  public void setScale(double scale) {
    double oldScale = this.scale;
    this.scale = scale;
    resizeAndValidate();
    firePropertyChange("scale", oldScale, scale);
  }

  /**
   * Retrieve the current scale.
   *
   * @return the number of pixels used to render one symbol
   */
  public double getScale() {
    return scale;
  }

  /**
   * Retrieve the object that encapsulates the leading border area - the space
   * before sequence information is rendered.
   *
   * @return a SequenceRenderContext.Border instance
   */
  public SequenceRenderContext.Border getLeadingBorder() {
    return leadingBorder;
  }

  /**
   * Retrieve the object that encapsulates the trailing border area - the space
   * after sequence information is rendered.
   *
   * @return a SequenceRenderContext.Border instance
   */
  public SequenceRenderContext.Border getTrailingBorder() {
    return trailingBorder;
  }

  /**
   * Paint this component.
   * <p>
   * This calls the paint method of the currently registered SequenceRenderer
   * after setting up the graphics appropriately.
   */
  public synchronized void paintComponent(Graphics g) {
          if(!isActive()) {
                  return;
          }

          Graphics2D g2 = (Graphics2D) g;
          if(hints != null){
            g2.setRenderingHints(hints);
          }
          super.paintComponent(g);


          AffineTransform oldTransform = g2.getTransform();
          //Rectangle2D currentClip = g2.getClip().getBounds2D();

          Insets insets = getInsets();

          if (isOpaque())
          {
                  g2.setPaint(getBackground());
                  g2.fillRect(0, 0, getWidth(), getHeight());
          }

          // do a transform to offset drawing to the neighbourhood of zero.
          adjustOffset(sequenceToGraphics(range.getMin()));

          double minAcross = sequenceToGraphics(range.getMin()) -
                  renderer.getMinimumLeader(this);
          double maxAcross = sequenceToGraphics(range.getMax()) + 1 +
                  renderer.getMinimumTrailer(this);
          double alongDim = maxAcross - minAcross;
          double depth = renderer.getDepth(this);
          Rectangle2D.Double clip = new Rectangle2D.Double();
          if (direction == HORIZONTAL) {
                  clip.x = minAcross;
                  clip.y = 0.0;
                  clip.width = alongDim;
                  clip.height = depth;
                  g2.translate(leadingBorder.getSize() - minAcross + insets.left, insets.top);
          } else {
                  clip.x = 0.0;
                  clip.y = minAcross;
                  clip.width = depth;
                  clip.height = alongDim;
                  g2.translate(insets.left, leadingBorder.getSize() - minAcross + insets.top);
          }

          Shape oldClip = g2.getClip();
          g2.clip(clip);
          renderer.paint(g2, new PaintContext());
          g2.setClip(oldClip);
          g2.setTransform(oldTransform);
  }

  public void setRenderer(SequenceRenderer r)
  throws ChangeVetoException {
    if(hasChangeListeners()) {
      ChangeEvent ce = new ChangeEvent(
        this,
        RENDERER,
        r,
        this.renderer
      );
      ChangeSupport cs = getChangeSupport(RENDERER);
      synchronized(cs) {
        cs.firePreChangeEvent(ce);
        _setRenderer(r);
        cs.firePostChangeEvent(ce);
      }
    } else {
      _setRenderer(r);
    }
    resizeAndValidate();
  }

  protected void _setRenderer(SequenceRenderer r) {
    if( (this.renderer != null) && (this.renderer instanceof Changeable) ) {
      Changeable c = (Changeable) this.renderer;
      c.removeChangeListener(layoutListener, SequenceRenderContext.LAYOUT);
      c.removeChangeListener(repaintListener, SequenceRenderContext.REPAINT);
    }

    this.renderer = r;

    if( (r != null) && (r instanceof Changeable) ) {
      Changeable c = (Changeable) r;
      c.addChangeListener(layoutListener, SequenceRenderContext.LAYOUT);
      c.addChangeListener(repaintListener, SequenceRenderContext.REPAINT);
    }
  }

  private void adjustOffset(double newOrigin) {
    pixelOffset -= newOrigin;
  }

  public double sequenceToGraphics(int seqPos) {
    return ((double) (seqPos-1)) * scale + pixelOffset;
  }

  public int graphicsToSequence(double gPos) {
    return ((int) ((gPos - pixelOffset) / scale)) + 1;
  }

  public int graphicsToSequence(Point2D point) {
    if(direction == HORIZONTAL) {
      return graphicsToSequence(point.getX());
    } else {
      return graphicsToSequence(point.getY());
    }
  }

  public void resizeAndValidate() {
    //System.out.println("resizeAndValidate starting");
    Dimension mind = null;
    Dimension maxd = null;

    if(!isActive()) {
      // System.out.println("No sequence");
      // no sequence - collapse down to no size at all
      leadingBorder.setSize(0.0);
      trailingBorder.setSize(0.0);
      mind = maxd = new Dimension(0, 0);
    } else {
      double minAcross = sequenceToGraphics(range.getMin());
      double maxAcross = sequenceToGraphics(range.getMax());
      double maxDropAcross = sequenceToGraphics(range.getMax() - 1);
      double lb = renderer.getMinimumLeader(this);
      double tb = renderer.getMinimumTrailer(this) + trailingBorder.getSize();
      double alongDim =
        (maxAcross - minAcross) +
        lb + tb;
      double alongDropDim =
    (maxDropAcross - minAcross) +
    lb + tb;
      double depth = renderer.getDepth(this);
      if(direction == HORIZONTAL) {
      mind = new Dimension((int) Math.ceil(alongDropDim), (int) Math.ceil(depth));
      maxd = new Dimension((int) Math.ceil(alongDim), (int) Math.ceil(depth));
      } else {
      mind = new Dimension((int) Math.ceil(depth), (int) Math.ceil(alongDropDim));
      maxd = new Dimension((int) Math.ceil(depth), (int) Math.ceil(alongDim));
      }
    }

    setMinimumSize(mind);
    setPreferredSize(maxd);
    setMaximumSize(maxd);
    revalidate();
    // System.out.println("resizeAndValidate ending");
  }

  private class RendererMonitor implements PropertyChangeListener {
    public void propertyChange(PropertyChangeEvent ev) {
      repaint();
    }
  }

    protected int [] calcDist() {
        double minAcross = sequenceToGraphics(range.getMin()) -
            renderer.getMinimumLeader(this);
        Insets insets = getInsets();

        int [] dist = new int[2];
        if(direction == HORIZONTAL) {
            dist[0] = (int) minAcross - insets.left;
            dist[1] = -insets.top;
        } else {
            dist[0] = -insets.left;
            dist[1] = (int) minAcross - insets.top;
        }

        return dist;
    }

  protected boolean isActive() {
    return
      (sequence != null) &&
      (renderer != null) &&
      (range != null);
  }

  public class Border
  implements Serializable, SwingConstants {
    protected final PropertyChangeSupport pcs;
    private double size = 0.0;
    private int alignment = CENTER;

    public double getSize() {
      return size;
    }

    public int getAlignment() {
      return alignment;
    }

    public void setAlignment(int alignment)
        throws IllegalArgumentException
    {
    if (alignment == LEADING || alignment == TRAILING || alignment == CENTER) {
        int old = this.alignment;
        this.alignment = alignment;
        pcs.firePropertyChange("alignment", old, alignment);
    } else {
        throw new IllegalArgumentException(
          "Alignment must be one of the constants LEADING, TRAILING or CENTER"
            );
    }
    }

    private Border() {
      alignment = CENTER;
      pcs = new PropertyChangeSupport(this);
    }

    public void addPropertyChangeListener(PropertyChangeListener listener) {
      pcs.addPropertyChangeListener(listener);
    }

    public void removePropertyChangeListener(PropertyChangeListener listener) {
      pcs.removePropertyChangeListener(listener);
    }
  }

    private boolean eq(Object a, Object b) {
    if (a == null || b == null) {
        return a == b;
    } else {
        return a.equals(b);
    }
    }


    public boolean equals(Object o) {
    if (! (o instanceof SequencePanel)) {
        return false;
    }

    SequencePanel osp = (SequencePanel) o;
    return (eq(getSymbols(), osp.getSymbols()) && eq(getRange(), osp.getRange()));
    }

    public int hashCode() {
    int hc = 653;
    SymbolList sl = getSymbols();
    if (sl != null) {
        hc = hc ^ sl.hashCode();
    }

    Location l = getRange();
    if (l != null) {
        hc = hc ^ l.hashCode();
    }

    return hc;
    }

  private class PaintContext
          implements SequenceRenderContext
  {
    private final RangeLocation range;

    public PaintContext() {
      this.range = (RangeLocation) LocationTools.intersection(
              SequencePanel.this.getRange(),
              new RangeLocation(1, SequencePanel.this.getSequence().length()));
    }

    public RangeLocation getRange()
    {
      return range;
    }

    public int getDirection()
    {
      return SequencePanel.this.getDirection();
    }

    public double getScale()
    {
      return SequencePanel.this.getScale();
    }

    public double sequenceToGraphics(int i)
    {
      return SequencePanel.this.sequenceToGraphics(i);
    }

    public int graphicsToSequence(double d)
    {
      return SequencePanel.this.graphicsToSequence(d);
    }

    public int graphicsToSequence(Point2D point)
    {
      return SequencePanel.this.graphicsToSequence(point);
    }

    public SymbolList getSymbols()
    {
      return SequencePanel.this.getSymbols();
    }

    public FeatureHolder getFeatures()
    {
      return SequencePanel.this.getFeatures();
    }

    public SequenceRenderContext.Border getLeadingBorder()
    {
      return SequencePanel.this.getLeadingBorder();
    }

    public SequenceRenderContext.Border getTrailingBorder()
    {
      return SequencePanel.this.getTrailingBorder();
    }

    public Font getFont()
    {
      return SequencePanel.this.getFont();
    }
  }
}
