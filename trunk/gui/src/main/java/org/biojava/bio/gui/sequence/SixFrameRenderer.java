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
import java.util.List;
import java.util.Vector;

import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.RangeLocation;
import org.biojava.utils.AbstractChangeable;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeVetoException;
 
/**
 * Class that handles drawing in six frames for other
 * classes.
 *
 * Suitable for use with CDS features particularly.
 * <p>
 * Partly based on Matthew Pocock's ZiggyFeatureRenderer.
 *
 * @author David Huen
 */

public class SixFrameRenderer 
    extends AbstractChangeable 
    implements SequenceRenderer {
  // this really ought to derive from a different class since
  // it doesn't implement a functional paint method.
  private double bandWidth = 25.0;
  private double blockWidth= 20.0;
  private Paint outline = Color.blue;
  private Paint fill = Color.red;
  private Paint lines = Color.black;

  // the following ought to be exported elsewhere later
  // this is not threadsafe!
  StrandedFeature.Strand strand;
  boolean  zigPrecedesBlock, phaseKnown;
  Point2D.Double zigOriginP;
  
  int zigOriginS;
  int moduloFrame;
//  double midpointOffset;
  double offsetPrevBlock;
  public SixFrameRenderer() {
  }

  public double getDepth(SequenceRenderContext src) {
    return 6.0 * bandWidth + 1;
  }

  public double getMinimumLeader(SequenceRenderContext src) {
    return 0.0;
  }
 
  public double getMinimumTrailer(SequenceRenderContext src) {
    return 0.0;
  }

  private double frameOffset(
                   int moduloFrame, 
                   StrandedFeature.Strand strand) {
    // computes the offset for a given frame
    if (strand == StrandedFeature.NEGATIVE) {
      return bandWidth * (moduloFrame + 3);
    }
    else {
      // default is positive strand too
      return bandWidth * moduloFrame;
    }
  }

  public void setFill(Paint p)
  throws ChangeVetoException {
    if(hasListeners()) {
      ChangeSupport cs = getChangeSupport(SequenceRenderContext.REPAINT);
      synchronized(cs) {
        ChangeEvent ce = new ChangeEvent(this, SequenceRenderContext.REPAINT);
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
 
  public void setBlockWidth(double width)
  throws ChangeVetoException {
    if(hasListeners()) {
      ChangeSupport cs = getChangeSupport(SequenceRenderContext.LAYOUT);
      synchronized(cs) {
        ChangeEvent ce = new ChangeEvent(this, SequenceRenderContext.LAYOUT);
        cs.firePreChangeEvent(ce);
        this.blockWidth = width;
        cs.firePostChangeEvent(ce);
      }
    } else {
      this.blockWidth = width;
    }
  }  

  public double getBlockWidth() {
    return blockWidth;
  }

/**
 * Used to initialise the spliced transcript renderer for
 * a CDS feature where the ends of the feature define the
 * frame of the feature.
 */
  public void startZiggy(StrandedFeature.Strand strand) {
   // initialise state variables
   System.out.println("startZiggy strand: " + strand);
   zigPrecedesBlock = false;
   phaseKnown = false;
   this.strand = strand;
  }

/**
 * This method is called to initialise the renderer for a
 * spliced transcript.
 *
 * @param strand the strand that the transcipt is on.
 * @param phase the initial translation phase of transcript coding region.
 *
 */
  public void startZiggy(StrandedFeature.Strand strand, int phase) {
   // initialise state variables
//   System.out.println("startZiggy strand: " + strand);
   zigPrecedesBlock = false;
   phaseKnown = true;
   moduloFrame = phase;
   this.strand = strand;
  }

/**
 * Render another "exon" in the correct translation frame.
 *
 */
  public void renderLocation(
           Graphics2D g,
           SequenceRenderContext src,
           Location loc) {

    int minS = loc.getMin();
    int maxS = loc.getMax();
    double minP = src.sequenceToGraphics(minS);
    double maxP = src.sequenceToGraphics(maxS);

    // handle the zig
    double midpointSeqCoord;
    double terminalSeqCoord;
    double terminalDrawCoord;
    double midpointOffset;
    double offset;

    if (zigPrecedesBlock) {
      // This is a continuation of the transcript.
      // do the ziggy
      // the apex is midway between blocks.
//      System.out.println("renderLocation continued moduloFrame: " + moduloFrame  + " " + minS + " " + maxS); 
        terminalSeqCoord = minP;

      // current phase is previous phase corrected by intron-induced
      // phase change.
      moduloFrame = (moduloFrame + minS - 1 - zigOriginS)%3;
      offset = frameOffset(moduloFrame, strand);

      // zigs for forward frames go up, reverse go down.
      if (strand == StrandedFeature.POSITIVE) {
        midpointOffset = Math.min(offset, offsetPrevBlock);
      }
      else {
        midpointOffset = Math.max(offset, offsetPrevBlock) + blockWidth;
      }

      terminalDrawCoord = offset + 0.5 * blockWidth;

      Point2D.Double midpoint, terminal;
      if (src.getDirection() == SequenceRenderContext.HORIZONTAL) {
        // horizontal axis
        midpointSeqCoord = 0.5*(zigOriginP.getX() + terminalSeqCoord);
        midpoint = new Point2D.Double(midpointSeqCoord, midpointOffset);
        terminal = new Point2D.Double(terminalSeqCoord, terminalDrawCoord);
      }
      else {
        midpointSeqCoord = 0.5*(zigOriginP.getY() + terminalSeqCoord);
        midpoint = new Point2D.Double(midpointOffset, midpointSeqCoord);
        terminal = new Point2D.Double(terminalDrawCoord, terminalSeqCoord);
      }

      // draw ziggy
      Paint prevPaint = g.getPaint();
      g.setPaint(outline);
      Line2D line = new Line2D.Double(zigOriginP, midpoint);
      g.draw(line);
      line = new Line2D.Double(midpoint, terminal);
      g.draw(line);
      g.setPaint(prevPaint);
    }
    else {
      // this is first block, there is no zig yet.
      // compute the frame.

      if (!phaseKnown)
          moduloFrame = minS%3;

      // compute offset for frame to be drawn
      offset = frameOffset(moduloFrame, strand);
      zigPrecedesBlock = true;
//      System.out.println("renderLocation 1st moduloFrame: " + moduloFrame + " " + minS + " " + maxS); 
    }

    // draw the block
    Rectangle2D.Double block = new Rectangle2D.Double();
    if (src.getDirection() == SequenceRenderContext.HORIZONTAL) {
      block.setFrame(
        minP, offset,
        maxP-minP, blockWidth);
    }
    else {
      block.setFrame(
        offset, minP,
        blockWidth, maxP-minP);
    }
    
    if (fill != null) {
      Paint prevPaint = g.getPaint();
      g.setPaint(fill);
      g.fill(block);
      g.setPaint(prevPaint);
    }

    if (outline != null) {
      Paint prevPaint = g.getPaint();
      g.setPaint(outline);
      g.draw(block);
      g.setPaint(prevPaint);
    }

    // update origin for next ziggy
    // origin of next ziggy on current block
    double seqCoord;
    double drawCoord = offset + 0.5 * blockWidth;

    // all features arrive in sequence order irrespective of strand
    // only the direction of zig changes.
    seqCoord = maxP;  

    // cache zig sequence origin and phase.
    zigOriginS = maxS;
    offsetPrevBlock = offset;
  
    if (src.getDirection() == SequenceRenderContext.HORIZONTAL) {
      zigOriginP = new Point2D.Double(seqCoord, drawCoord);
    }
    else {
      zigOriginP = new Point2D.Double(drawCoord, seqCoord);
    }
  }

  public java.util.List sequenceExtentOfPixels(
                          SequenceRenderContext src) {
    // returns a list giving the extents represented by
    // pixels in drawing area.
    RangeLocation rangeS = src.getRange();
    int minLocS = rangeS.getMin();
    int maxLocS = rangeS.getMax();
    int minP = (int) src.sequenceToGraphics(minLocS);
    int maxP = (int) src.sequenceToGraphics(maxLocS);
    // Vector to hold ranges for each pixel.  There are maxP-minP+1 pixels.
    java.util.List extentList = new Vector(maxP - minP +1, 20);
    for (int currP = minP; currP < maxP; currP++) {
      // compute the extent of represented by each pixel
      int minRange = Math.max((int) src.graphicsToSequence(currP), minLocS);
      int maxRange = Math.min(((int) src.graphicsToSequence(currP + 1) - 1), maxLocS);
      extentList.add(new RangeLocation(minRange, maxRange));
    }
    return extentList;
  }

/** 
 * draws required bar in correct translation frame.
 *
 * @param base     the sequence coordinate of the bar.
 * @param strand   strand that the bar annotates.
 */
  public void drawLine(
                Graphics2D g,
                SequenceRenderContext src,
                int base, 
                StrandedFeature.Strand strand) {

    Paint prevPaint = g.getPaint();
    g.setPaint(lines);

    // compute the frame to use.
    int moduloFrame = base%3;

//    System.out.println("drawLine: base,strand,modulo" + base + " " + strand + " " + moduloFrame);
    // get required offset for frame
    double offset = frameOffset(moduloFrame, strand);

    // compute position of line to be drawn
    int lineP = (int) src.sequenceToGraphics(base);

    // draw the line
    if (src.getDirection() == SequenceRenderContext.HORIZONTAL) {
      g.drawLine(lineP, (int) offset,
                 lineP, (int) (offset + blockWidth));
    }
    else {
      g.drawLine((int) offset, lineP,
                 (int)(offset + blockWidth), lineP);
    }
    g.setPaint(prevPaint);
  }

  public void paint(
    Graphics2D g,
    SequenceRenderContext src) {
    // this doesn't do anything
  }

  public SequenceViewerEvent processMouseEvent(
    SequenceRenderContext src,
    MouseEvent me,
    List path
  ) {
    path.add(this);
    int sPos = src.graphicsToSequence(me.getPoint());
    return new SequenceViewerEvent(this, null, sPos, me, path);
  }
}

