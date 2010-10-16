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

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Paint;
import java.awt.Shape;
import java.awt.event.MouseEvent;
import java.awt.geom.GeneralPath;

import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.FeatureHolder;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.symbol.Location;
import org.biojava.utils.AbstractChangeable;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;

/**
 * @author Thomas Down
 * @author Matthew Pocock
 * @author David Huen
 */
public class TickFeatureRenderer
extends AbstractChangeable
implements FeatureRenderer {
  public static final ChangeType FILL = new ChangeType(
    "The fill paint has changed",
    "org.biojava.bio.gui.sequence.TickFeatureRenderer",
    "FILL"
  );

  public static final ChangeType OUTLINE = new ChangeType(
    "The outline paint has changed",
    "org.biojava.bio.gui.sequence.TickFeatureRenderer",
    "OUTLINE"
  );

  public static final ChangeType DEPTH = new ChangeType(
    "The size of the arrow has changed",
    "org.biojava.bio.gui.sequence.TickFeatureRenderer",
    "DEPTH"
  );

  private Paint fill;
  private Paint outline;
  private double depth = 25.0;

  public TickFeatureRenderer() {
    fill = Color.blue;
    outline = Color.black;
  }

  public void setFill(Paint p)
  throws ChangeVetoException {
    if(hasListeners()) {
      ChangeSupport cs = getChangeSupport(SequenceRenderContext.REPAINT);
      synchronized(cs) {
        ChangeEvent ce = new ChangeEvent(
          this, SequenceRenderContext.REPAINT,
          null, null, new ChangeEvent(
            this, FILL, p, fill
          )
        );
        cs.firePreChangeEvent(ce);
        fill = p;
        cs.firePostChangeEvent(ce);
      }
    } else {
      fill = p;
    }
  }

  public Paint getFill() {
    return fill;
  }

  public void setOutline(Paint p)
  throws ChangeVetoException {
    if(hasListeners()) {
      ChangeSupport cs = getChangeSupport(SequenceRenderContext.REPAINT);
      synchronized(cs) {
        ChangeEvent ce = new ChangeEvent(
          this, SequenceRenderContext.REPAINT,
          null, null, new ChangeEvent(
            this, OUTLINE, p, outline
          )
        );
        cs.firePreChangeEvent(ce);
        outline = p;
        cs.firePostChangeEvent(ce);
      }
    } else {
      outline = p;
    }
  }

  public Paint getOutline() {
    return outline;
  }

  public void setDepth(double arrowSize)
  throws ChangeVetoException {
    if(hasListeners()) {
      ChangeSupport cs = getChangeSupport(SequenceRenderContext.LAYOUT);
      synchronized(cs) {
        ChangeEvent ce = new ChangeEvent(
          this, SequenceRenderContext.LAYOUT,
          null, null, new ChangeEvent(
            this, DEPTH, new Double(this.depth), new Double(arrowSize)
          )
        );
        cs.firePreChangeEvent(ce);
        this.depth = arrowSize;
        cs.firePostChangeEvent(ce);
      }
    } else {
      this.depth = arrowSize;
    }
  }

  public double getDepth() {
    return depth;
  }

  public void renderFeature(
    Graphics2D g,
    Feature f,
    SequenceRenderContext src
  ) {
    Shape s = null;
    Location loc = f.getLocation();
    float min = (float) src.sequenceToGraphics(loc.getMin());
    float max = (float) src.sequenceToGraphics(loc.getMax());
    float pos = (min + max) / 2;

    float fDepth = (float) depth;
    float fDepthByThree = fDepth / 3.0F;

    if (f instanceof StrandedFeature) {
      StrandedFeature.Strand strand = ((StrandedFeature) f).getStrand();
      if(src.getDirection() == SequenceRenderContext.HORIZONTAL) {
        if(strand == StrandedFeature.POSITIVE) {
          GeneralPath path = new GeneralPath();
          path.moveTo(pos, 0.0F);
          path.lineTo(pos, fDepth);
          path.lineTo(pos + fDepthByThree, fDepth);
          path.lineTo(pos, fDepth - fDepthByThree);
          path.closePath();
          s = path;
        } else if(strand == StrandedFeature.NEGATIVE) {
          GeneralPath path = new GeneralPath();
          path.moveTo(pos, 0.0F);
          path.lineTo(pos, fDepth);
          path.lineTo(pos - fDepthByThree, fDepth);
          path.lineTo(pos, fDepth - fDepthByThree);
          path.closePath();
          s = path;
        }
      } else { // vertical
        if(strand == StrandedFeature.POSITIVE) {
          GeneralPath path = new GeneralPath();
          path.moveTo(0.0F, pos);
          path.lineTo(fDepth, pos);
          path.lineTo(fDepth, pos + fDepthByThree);
          path.lineTo(fDepth - fDepthByThree, pos);
          path.closePath();
          s = path;
        } else if(strand == StrandedFeature.NEGATIVE) {
          GeneralPath path = new GeneralPath();
          path.moveTo(0.0F, pos);
          path.lineTo(fDepth, pos);
          path.lineTo(fDepth, pos - fDepthByThree);
          path.lineTo(fDepth - fDepthByThree, pos);
          path.closePath();
          s = path;
        }
      }
    }

    if(fill != null) {
      Paint prevPaint = g.getPaint();
      g.setPaint(fill);
      g.fill(s);
      g.setPaint(prevPaint);
    }
    if (outline != null) {
      Paint prevPaint = g.getPaint();
      g.setPaint(outline);
      g.draw(s);
      g.setPaint(prevPaint);
    }
  }

  public double getDepth(SequenceRenderContext src) {
      return depth;
  }

  public FeatureHolder processMouseEvent(
    FeatureHolder hits,
    SequenceRenderContext src,
    MouseEvent me
  ) {
    return hits;
  }
}
