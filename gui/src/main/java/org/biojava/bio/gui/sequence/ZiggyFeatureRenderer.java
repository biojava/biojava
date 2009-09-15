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
import java.awt.event.MouseEvent;
import java.awt.geom.Line2D;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.util.Iterator;

import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.FeatureHolder;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.symbol.Location;
import org.biojava.utils.AbstractChangeable;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeVetoException;

/**
 * A feature renderer that draws non-contiguous features as a set of boxes
 * joined by zig-zags.
 * <p>
 * This is applicable to rendering cds's or non-contiguous homologies for
 * example.
 *
 * @author Matthew Pocock
 */
public class ZiggyFeatureRenderer
extends AbstractChangeable
implements FeatureRenderer, java.io.Serializable {
  private Paint outline = Color.black;
  private Paint fill = Color.yellow;
  private double blockDepth = 10.0;

  public void setFill(Paint p)
  throws ChangeVetoException {
    if(hasListeners()) {
      ChangeSupport cs = getChangeSupport(SequenceRenderContext.REPAINT);
      synchronized(cs) {
        ChangeEvent ce = new ChangeEvent(
          this, SequenceRenderContext.REPAINT,
          p, this.fill
        );
        cs.firePreChangeEvent(ce);
        this.fill = p;
        cs.firePostChangeEvent(ce);
      }
    } else {
      this.fill = p;
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
        ChangeEvent ce = new ChangeEvent(this, SequenceRenderContext.REPAINT);
        cs.firePreChangeEvent(ce);
        this.outline = p;
        cs.firePostChangeEvent(ce);
      }
    } else {
      this.outline = p;
    }
  }

  public Paint getOutline() {
    return outline;
  }

  public void setBlockDepth(double depth)
  throws ChangeVetoException {
    if(hasListeners()) {
      ChangeSupport cs = getChangeSupport(SequenceRenderContext.LAYOUT);
      synchronized(cs) {
        ChangeEvent ce = new ChangeEvent(this, SequenceRenderContext.LAYOUT);
        cs.firePreChangeEvent(ce);
        this.blockDepth = depth;
        cs.firePostChangeEvent(ce);
      }
    } else {
      this.blockDepth = depth;
    }
  }

  public double getBlockDepth() {
    return blockDepth;
  }

  public double getDepth(SequenceRenderContext src) {
    return blockDepth + 1.0;
  }

  public void renderFeature(
    Graphics2D g, Feature f, SequenceRenderContext context
  ) {
    Location loc = f.getLocation();
    Iterator i = loc.blockIterator();
    Location last = null;
    if(i.hasNext()) {
      last = (Location) i.next();
      renderLocation(g, last, context);
    }
    while(i.hasNext()) {
      Location next = (Location) i.next();
      renderLink(g, f, last, next, context);
      renderLocation(g, next, context);
      last = next;
    }
  }

  private void renderLocation(
    Graphics2D g, Location loc, SequenceRenderContext context
  ) {
    Rectangle2D.Double block = new Rectangle2D.Double();
    double min = context.sequenceToGraphics(loc.getMin());
    double max = context.sequenceToGraphics(loc.getMax()+1);
    if(context.getDirection() == SequenceRenderContext.HORIZONTAL) {
      block.setFrame(
        min, 0.0,
        max - min, blockDepth
      );
    } else {
      block.setFrame(
        0.0, min,
        blockDepth, max - min
      );
    }
    if(fill != null) {
      g.setPaint(fill);
      g.fill(block);
    }
    if(outline != null) {
      g.setPaint(outline);
      g.draw(block);
    }
  }

    private StrandedFeature.Strand getStrand(Feature f) {
	if (f instanceof StrandedFeature) {
	    return ((StrandedFeature) f).getStrand();
	} else {
	    FeatureHolder fh = f.getParent();
	    if (fh instanceof Feature) {
		return getStrand((Feature) fh);
	    } else {
		return StrandedFeature.UNKNOWN;
	    }
	}
    }

  private void renderLink(
    Graphics2D g, Feature f, Location source, Location dest,
    SequenceRenderContext context
  ) {
    Line2D line = new Line2D.Double();
    Point2D startP;
    Point2D midP;
    Point2D endP;
    double half = blockDepth * 0.5;
    if(context.getDirection() == SequenceRenderContext.HORIZONTAL) {
      if(getStrand(f) == StrandedFeature.NEGATIVE) {
        double start = context.sequenceToGraphics(dest.getMin());
        double end = context.sequenceToGraphics(source.getMax()+1);
        double mid = (start + end) * 0.5;
        startP = new Point2D.Double(start, half);
        midP   = new Point2D.Double(mid,   blockDepth);
        endP   = new Point2D.Double(end,   half);
      } else {
        double start = context.sequenceToGraphics(source.getMax()+1);
        double end = context.sequenceToGraphics(dest.getMin());
        double mid = (start + end) * 0.5;
        startP = new Point2D.Double(start, half);
        midP   = new Point2D.Double(mid,   0.0);
        endP   = new Point2D.Double(end,   half);
      }
    } else {
      if (getStrand(f) == StrandedFeature.NEGATIVE) {
        double start = context.sequenceToGraphics(dest.getMin()+1);
        double end = context.sequenceToGraphics(source.getMax());
        double mid = (start + end) * 0.5;
        startP = new Point2D.Double(half,       start);
        midP   = new Point2D.Double(blockDepth, mid);
        endP   = new Point2D.Double(half,       end);
      } else {
        double start = context.sequenceToGraphics(source.getMax());
        double end = context.sequenceToGraphics(dest.getMin()+1);
        double mid = (start + end) * 0.5;
        startP = new Point2D.Double(half, start);
        midP   = new Point2D.Double(0.0,  mid);
        endP   = new Point2D.Double(half, end);
      }
    }
    g.setPaint(getOutline());
    line.setLine(startP, midP);
    g.draw(line);
    line.setLine(midP, endP);
    g.draw(line);
  }

  public FeatureHolder processMouseEvent(
    FeatureHolder hits,
    SequenceRenderContext src,
    MouseEvent me
  ) {
    return hits;
  }
}
