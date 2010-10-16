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
import java.awt.geom.AffineTransform;
import java.util.List;

import org.biojava.utils.AssertionFailure;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;

/**
 * A renderer that adds padding before and after a delegate renderer.
 *
 * @author Matthew Pocock
 * @since 1.2
 */ 
public class PaddingRenderer
extends SequenceRendererWrapper {
  /**
   * Event type for when the size of the padding changes.
   */
  public static ChangeType PADDING = new ChangeType(
    "The padding has changed",
    "org.biojava.bio.gui.sequence.PaddingRenderer",
    "PADDING",
    SequenceRenderContext.LAYOUT
  );

  private double padding;

  protected boolean hasListeners() {
    return super.hasListeners();
  }

  protected ChangeSupport getChangeSupport(ChangeType ct) {
    return super.getChangeSupport(ct);
  }
  
  /**
   * Build a new PaddingRenderer with zero padding.
   * <p>
   * This will cause a rendering effect equivalent to missing out the padding
   * renderer all together.
   */
  public PaddingRenderer() {
    padding = 0.0;
  }
  
  /**
   * Build a new PaddingRenderer that wraps <code>renderer</code> and has
   * padding depth <code>padding</code>.
   *
   * @param renderer  the SequenceRenderer that will actually do the rendering
   * @param padding  the number of pixels to leave both before and after
   *        rendering the child renderer
   */
  public PaddingRenderer(
    SequenceRenderer renderer,
    double padding
  ) {
    super(renderer);
    try {
      setPadding(padding);
    } catch (ChangeVetoException cve) {
      throw new AssertionFailure("Assertion Failure: Should have no listeners", cve);
    }
  }
  
  /**
   * Set the padding.
   * <p>
   * The padding will be added to the area before and after that required to
   * render the delegate renderer.
   *
   * @param padding  the new padding size
   * @throws ChangeVetoException if padding is negative or if any listener
   *         objected to the change
   */
  public void setPadding(double padding)
  throws ChangeVetoException {
    if(padding < 0.0) {
      ChangeEvent ce = new ChangeEvent(
        this, PADDING, new Double(this.padding), new Double(padding)
      );
      throw new ChangeVetoException(
        ce,
        "Can't set padding to a negative value"
      );
    }
    if(hasListeners()) {
      ChangeSupport cs = getChangeSupport(PADDING);
      synchronized(cs) {
        ChangeEvent ce = new ChangeEvent(
          this, PADDING, new Double(this.padding), new Double(padding)
        );
        cs.firePreChangeEvent(ce);
        this.padding = padding;
        cs.firePostChangeEvent(ce);
      }
    } else {
      this.padding = padding;
    }
  }
  
  /**
   * Retrieve the current padding.
   *
   * @return the current padding
   */
  public double getPadding() {
    return this.padding;
  }
  
  public double getDepth(SequenceRenderContext src) {
    return super.getDepth(src) + padding * 2.0;
  }    
  
  public double getMinimumLeader(SequenceRenderContext src) {
    return super.getMinimumLeader(src);
  }
  
  public double getMinimumTrailer(SequenceRenderContext src) {
    return super.getMinimumTrailer(src);
  }
  
  public void paint(
    Graphics2D g,
    SequenceRenderContext src
  ) {
    AffineTransform old = g.getTransform();
    if(src.getDirection() == SequenceRenderContext.HORIZONTAL) {
      g.translate(0.0, getPadding());
    } else {
      g.translate(getPadding(), 0.0);
    }
    super.paint(g, src);
    g.setTransform(old);
  }
  
  public SequenceViewerEvent processMouseEvent(
    SequenceRenderContext src,
    MouseEvent me,
    List path
  ) {
    int padding = (int) getPadding();
    if(src.getDirection() == SequenceRenderContext.HORIZONTAL) {
      me.translatePoint(0, -padding);
    } else {
      me.translatePoint(-padding, 0);
    }
    SequenceViewerEvent sve = super.processMouseEvent(
      src,
      me,
      path
    );
    if(src.getDirection() == SequenceRenderContext.HORIZONTAL) {
      me.translatePoint(0, padding);
    } else {
      me.translatePoint(padding, 0);
    }
    return sve;
  }

    public String toString() {
	return "PaddingRenderer(" + getRenderer().toString() + ")";
    }
}
