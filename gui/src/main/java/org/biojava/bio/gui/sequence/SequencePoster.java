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
import java.awt.RenderingHints;
import java.awt.Shape;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.beans.PropertyChangeSupport;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;

import javax.swing.JComponent;
import javax.swing.SwingConstants;

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
 * A panel that displays a Sequence.
 * <p>
 * A SequencePoster can either display the sequence from left-to-right
 * (HORIZONTAL) or from top-to-bottom (VERTICAL). It has an associated scale
 * which is the number of pixels per symbol. It also has a lines property that
 * controls how to wrap the sequence off one end and onto the other.
 * <p>
 * Each line in the SequencePoster is broken down into a list of strips,
 * each rendered by an individual SequenceRenderer object.
 * You could add a SequenceRenderer that draws on genes, another that
 * draws repeats and another that prints out the DNA sequence. They are
 * responsible for rendering their view of the sequence in the place that the
 * SequencePoster positions them.
 *
 * @author Thomas Down
 * @author Matthew Pocock
 * @deprecated This doesn't handle loads of stuff. Use SequencePoster.
 */
public class SequencePoster
  extends
    JComponent
  implements
    SwingConstants,
    SequenceRenderContext,
    Changeable
{
  public static final ChangeType RENDERER = new ChangeType(
    "The renderer for this SequencePoster has changed",
    "org.biojava.bio.gui.sequence.SequencePoster",
    "RENDERER",
    SequenceRenderContext.LAYOUT
  );

  private Sequence sequence;
  private int direction;
  private double scale;
  private int lines;
  private int spacer;

  private SequenceRenderContext.Border leadingBorder;
  private SequenceRenderContext.Border trailingBorder;

  private SequenceRenderer renderer;
  private double[] offsets;
  private int realLines;
  private double alongDim = 0.0;
  private double acrossDim = 0.0;
  private int symbolsPerLine = 0;

  private RendererMonitor theMonitor;
  private RenderingHints renderingHints = null;

  private transient ChangeSupport changeSupport = null;

  private SequenceViewerSupport svSupport = new SequenceViewerSupport();
  private MouseListener mouseListener = new MouseAdapter() {
    public void mouseClicked(MouseEvent me) {
      if(!isActive()) {
        return;
      }
      int[] lineExtent = calcLineExtent(me);
      me.translatePoint(-lineExtent[2], -lineExtent[3]);
      SequenceViewerEvent sve = renderer.processMouseEvent(
        new SubSequenceRenderContext(
          SequencePoster.this, null, null,
          new RangeLocation(lineExtent[0], lineExtent[1])
        ),
        me,
        new ArrayList()
      );
      me.translatePoint(+lineExtent[2], +lineExtent[3]);
      svSupport.fireMouseClicked(sve);
    }

    public void mousePressed(MouseEvent me) {
      if(!isActive()) {
        return;
      }
      int[] lineExtent = calcLineExtent(me);
      me.translatePoint(-lineExtent[2], -lineExtent[3]);
      SequenceViewerEvent sve = renderer.processMouseEvent(
        new SubSequenceRenderContext(
          SequencePoster.this, null, null,
          new RangeLocation(lineExtent[0], lineExtent[1])
        ),
        me,
        new ArrayList()
      );
      me.translatePoint(+lineExtent[2], +lineExtent[3]);
      svSupport.fireMousePressed(sve);
    }

    public void mouseReleased(MouseEvent me) {
      if(!isActive()) {
        return;
      }
      int[] lineExtent = calcLineExtent(me);
      me.translatePoint(-lineExtent[2], -lineExtent[3]);
      SequenceViewerEvent sve = renderer.processMouseEvent(
        new SubSequenceRenderContext(
          SequencePoster.this, null, null,
          new RangeLocation(lineExtent[0], lineExtent[1])
        ),
        me,
        new ArrayList()
      );
      me.translatePoint(+lineExtent[2], +lineExtent[3]);
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
      int[] lineExtent = calcLineExtent(me);
      me.translatePoint(-lineExtent[2], -lineExtent[3]);
      SequenceViewerEvent sve = renderer.processMouseEvent(
        new SubSequenceRenderContext(
          SequencePoster.this, null, null,
          new RangeLocation(lineExtent[0], lineExtent[1])
        ),
        me,
        new ArrayList()
      );
      me.translatePoint(+lineExtent[2], +lineExtent[3]);
      svmSupport.fireMouseDragged(sve);
    }

    public void mouseMoved(MouseEvent me) {
      if(!isActive()) {
        return;
      }
      int[] lineExtent = calcLineExtent(me);
      me.translatePoint(-lineExtent[2], -lineExtent[3]);
      SequenceViewerEvent sve = renderer.processMouseEvent(
        new SubSequenceRenderContext(
          SequencePoster.this, null, null,
          new RangeLocation(lineExtent[0], lineExtent[1])
        ),
        me,
        new ArrayList()
      );
      me.translatePoint(+lineExtent[2], +lineExtent[3]);
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
    if(hasChangeListeners()) {
      getChangeSupport(ct).removeChangeListener(cl, ct);
    }
  }

  public boolean isUnchanging(ChangeType ct) {
    return getChangeSupport(ct).isUnchanging(ct);
  }

  private ChangeListener layoutListener = new ChangeAdapter() {
    public void postChange(ChangeEvent ce) {
      resizeAndValidate();
    }
  };
  private ChangeListener repaintListener = new ChangeAdapter() {
    public void postChange(ChangeEvent ce) {
      repaint();
    }
  };

  /**
   * Initializer.
   */

  {
    direction = HORIZONTAL;
    scale = 12.0;
    lines = 1;
    spacer = 0;

    theMonitor = new RendererMonitor();
    leadingBorder = new SequenceRenderContext.Border();
    trailingBorder = new SequenceRenderContext.Border();
  }

  /**
   * Create a new SeqeuncePanel.
   */
  public SequencePoster() {
    super();
    if(getFont() == null) {
      setFont(new Font("Times New Roman", Font.PLAIN, 12));
    }
    this.addPropertyChangeListener(theMonitor);
    this.addMouseListener(mouseListener);
    this.addMouseMotionListener(mouseMotionListener);
  }

  /**
   * Set the SymboList to be rendered. This symbol list will be passed onto the
   * SequenceRenderer instances registered with this SequencePoster.
   *
   * @param s  the SymboList to render
   */
  public void setSequence(Sequence s) {
    Sequence oldSequence = sequence;
    if(oldSequence != null) {
      oldSequence.removeChangeListener(layoutListener);
    }
    this.sequence = s;
    sequence.addChangeListener(layoutListener);

    resizeAndValidate();
    firePropertyChange("sequence", oldSequence, s);
  }

  public Sequence getSequence() {
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
    return sequence;
  }

  public RangeLocation getRange() {
    return new RangeLocation(1, sequence.length());
  }

  public RangeLocation getVisibleRange() {
    return getRange();
  }

  /**
   * Set the direction that this SequencePoster renders in. The direction can be
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
   * Set the number of pixels to leave blank between each block of sequence
   * information.
   * <p>
   * If the SeqeuncePanel chooses to display the sequence information split
   * across multiple lines, then the spacer parameter indicates how many pixles
   * will seperate each line.
   *
   * @param spacer  the number of pixels seperating each line of sequence
   * information
   */
  public void setSpacer(int spacer) {
    int oldSpacer = this.spacer;
    this.spacer = spacer;
    resizeAndValidate();
    firePropertyChange("spacer", oldSpacer, spacer);
  }

  /**
   * Retrieve the current spacer value
   *
   * @return the number of pixels between each line of sequence information
   */
  public int getSpacer() {
    return spacer;
  }

  /**
   * Set the scale.
   * <p>
   * The scale parameter is interpreted as the number of pixels per symbol. This
   * may take on a wide range of values - for example, to render the symbols as
   * text, you will need a scale of > 8, where as to render chromosome 1 you
   * will want a scale &lt; 0.00000001
   *
   * @param scale the new pixles-per-symbol ratio
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
   * @return the number of pixles used to render one symbol
   */
  public double getScale() {
    return scale;
  }

  /**
   * Set the absolute number of lines that the sequence will be rendered on. If
   * this is set to 0, then the number of lines will be calculated according to
   * how many lines will be needed to render the sequence in the currently
   * available space. If it is set to any positive non-zero value, the sequence
   * will be rendered using that many lines, and the SequencePoster will request
   * enough space to accomplish this.
   *
   * @param lines  the number of lines to split the sequence information over
   */
  public void setLines(int lines) {
    int oldLines = this.lines;
    this.lines = lines;
    resizeAndValidate();
    firePropertyChange("lines", oldLines, lines);
  }

  /**
   * Retrieve the number of lines that the sequence will be rendered over.
   *
   * @return the current number of lines (0 if autocalculated)
   */
  public int getLines() {
    return lines;
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
  public void paintComponent(Graphics g) {
    if(!isActive()) {
      return;
    }

    Graphics2D g2 = (Graphics2D) g;
    if(renderingHints != null){
      g2.setRenderingHints(renderingHints);
    }

    Rectangle2D currentClip = g2.getClip().getBounds2D();
    double minPos;
    double maxPos;
    if(direction == HORIZONTAL) {
      minPos = currentClip.getMinY();
      maxPos = currentClip.getMaxY();
    } else {
      minPos = currentClip.getMinX();
      maxPos = currentClip.getMaxX();
    }

    //System.out.println("minPos: " + minPos);
    int minOffset = Arrays.binarySearch(offsets, minPos);
    if(minOffset < 0) {
      minOffset = -minOffset - 1;
    }
    //System.out.println("minOffset: " + minOffset);
    double minCoord = (minOffset == 0) ? 0.0 : offsets[minOffset-1];
    //System.out.println("minCoord: " + minCoord);
    int minP = 1 + (int) ((double) minOffset * symbolsPerLine);
    //System.out.println("minP: " + minP);

    Rectangle2D.Double clip = new Rectangle2D.Double();
    if (direction == HORIZONTAL) {
        clip.width = alongDim;
        clip.height = acrossDim;
        g2.translate(leadingBorder.getSize() - alongDim * minOffset, minCoord);
    } else {
        clip.width = acrossDim;
        clip.height = alongDim;
        g2.translate(minCoord, leadingBorder.getSize() - alongDim * minOffset);
    }

    int min = minP;
    for(int l = minOffset; l < realLines; l++) {

      if (direction == HORIZONTAL) {
          clip.x = l * alongDim;
          clip.y = 0.0;
      } else {
          clip.x = 0.0;
          clip.y = l * alongDim;
      }

      double depth = offsets[l] - spacer;
      if(l != 0) {
        depth -= offsets[l-1];
      }

      if (direction == HORIZONTAL) {
          clip.height = depth;
      } else {
          clip.width = depth;
      }

      Shape oldClip = g2.getClip();
      g2.clip(clip);
      renderer.paint(g2, this);
      g2.setClip(oldClip);

      if (direction == HORIZONTAL) {
          g2.translate(-alongDim, spacer + depth);
      } else {
          g2.translate(spacer + depth, -alongDim);
      }

      min += symbolsPerLine;

      if( (min > sequence.length()) || (offsets[l] > maxPos)) {
        break;
      }
    }
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

  public double sequenceToGraphics(int seqPos) {
    return ((double) (seqPos-1) * scale);
  }

  public int graphicsToSequence(double gPos) {
    return (int) (gPos / scale) + 1;
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
    Dimension d = null;
    double acrossDim;

    if(!isActive()) {
      System.out.println("No sequence");
      // no sequence - collapse down to no size at all
      alongDim = 0.0;
      acrossDim = 0.0;
      realLines = 0;
      leadingBorder.setSize(0.0);
      trailingBorder.setSize(0.0);
      d = new Dimension(0, 0);
    } else {
      System.out.println("Fitting to sequence");

      int width;
      Dimension parentSize = (getParent() != null)
                ? getParent().getSize()
                : new Dimension(500, 400);
      if (direction == HORIZONTAL) {
        width = parentSize.width;
      } else {
        width = parentSize.height;
      }

      System.out.println("Initial width: " + width);
      // got a sequence - fit the size according to sequence length & preferred
      // number of lines.
      alongDim = scale * sequence.length();
      System.out.println("alongDim (pixles needed for sequence only): "
      + alongDim);
      acrossDim = 0.0;

      double insetBefore = renderer.getMinimumLeader(this);
      double insetAfter = renderer.getMinimumTrailer(this);

      leadingBorder.setSize(insetBefore);
      trailingBorder.setSize(insetAfter);
      double insets = insetBefore + insetAfter;
      //System.out.println("insetBefore: " + insetBefore);
      //System.out.println("insetAfter: " + insetAfter);

      if(lines > 0) {
        // Fixed number of lines. Calculate width needed to lay out rectangle.
        realLines = lines;
        width = (int) Math.ceil(
          insets +
          alongDim / (double) lines
        );
      } else {
        // Calculated number of lines for a fixed width
        double dWidth = (double) width;
        dWidth -= insets; // leave space for insets
        realLines = (int) Math.ceil(alongDim / (double) width);
        width = (int) Math.ceil(
          insets +
          alongDim / (double) realLines
        );
      }

      acrossDim = 0.0;
      symbolsPerLine = (int) Math.ceil((double) width / (double) scale);
      //System.out.println("symbolsPerLine: " + symbolsPerLine);
      //System.out.println("width: " + width);
      //System.out.println("lines: " + lines);
      //System.out.println("realLines: " + realLines);
      if(symbolsPerLine < 1) {
        throw new Error("Pants");
      }
      int min = 1;
      this.offsets = new double[realLines];
      int li = 0;
      while(min <= sequence.length()) {
        int max = min + symbolsPerLine - 1;
        double depth = renderer.getDepth(this);
        acrossDim += depth + spacer;
        offsets[li] = acrossDim;
        min = max + 1;
        li++;
      }

      acrossDim += spacer * (realLines - 1);
      alongDim = /* Math.ceil((double) width); */ symbolsPerLine * scale;
      if (direction == HORIZONTAL) {
        d = new Dimension(
          (int) Math.ceil(alongDim + insetBefore + insetAfter),
          (int) acrossDim
        );
      } else {
        d = new Dimension(
          (int) acrossDim,
          (int) Math.ceil(alongDim + insetBefore + insetAfter)
        );
      }
    }

    setMinimumSize(d);
    setPreferredSize(d);
    revalidate();
    //System.out.println("resizeAndValidate ending");
  }

  private class RendererMonitor implements PropertyChangeListener {
    public void propertyChange(PropertyChangeEvent ev) {
      repaint();
    }
  }

  protected int[] calcLineExtent(MouseEvent me) {
    int pos;
    if(direction == HORIZONTAL) {
      pos = me.getY();
    } else {
      pos = me.getX();
    }

    int minOffset = Arrays.binarySearch(offsets, pos);
    if(minOffset < 0) {
      minOffset = -minOffset - 1;
    }
    int min = 1 + (int) ((double) minOffset * symbolsPerLine);
    int max = min + symbolsPerLine - 1;
    double minPos;
    if(minOffset > 0) {
      minPos = offsets[minOffset - 1];
    } else {
      minPos = 0.0;
    }

    double ad = alongDim * minOffset;

    int xdiff;
    int ydiff;
    if(direction == HORIZONTAL) {
      xdiff = (int) -ad;
      ydiff = (int) minPos;
    } else {
      xdiff = (int) minPos;
      ydiff = (int) -ad;
    }

    return new int[] { min, max, xdiff, ydiff };
  }

  protected boolean isActive() {
    return
      (sequence != null) &&
      (renderer != null);
  }

  /**
   * @return the current rendering properties
   */
  public RenderingHints getRenderingHints() {
    return renderingHints;
  }
  /**
   * Use this to switch on effects like Anti-aliasing etc
   * @param renderingHints the desired rendering properties
   */
  public void setRenderingHints(RenderingHints renderingHints) {
    this.renderingHints = renderingHints;
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
}

