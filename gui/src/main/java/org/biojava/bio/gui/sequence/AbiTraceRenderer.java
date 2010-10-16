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
import java.awt.event.MouseEvent;
import java.awt.geom.Line2D;
import java.awt.geom.Rectangle2D;
import java.util.List;

import org.biojava.bio.BioError;
import org.biojava.bio.program.abi.ABITrace;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.utils.AbstractChangeable;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;

/**
 * Renders an ABI trace file as a chromatogram graph.
 *
 * @author Matthew Pocock
 * @author Mark Schreiber
 */
public class AbiTraceRenderer
extends AbstractChangeable
implements SequenceRenderer {
  public static final ChangeType TRACE = new ChangeType(
    "The trace has changed",
    AbiTraceRenderer.class,
    "TRACE",
    SequenceRenderContext.LAYOUT
  );

  public static final ChangeType DEPTH = new ChangeType(
    "The trace render depth has changed",
    AbiTraceRenderer.class,
    "DEPTH",
    SequenceRenderContext.LAYOUT
  );

  private ABITrace trace;
  private double depth;

  public AbiTraceRenderer() {
  }

  public void paint(Graphics2D g, SequenceRenderContext ctxt) {
    if(ctxt.getDirection() == SequenceRenderContext.VERTICAL || trace == null) {
      return;
    }

    try {
      Rectangle2D clip = g.getClip().getBounds2D();

      int min = Math.max(ctxt.getRange().getMin(), ctxt.graphicsToSequence(clip.getMinX()));
      int max = Math.min(ctxt.getRange().getMax(), ctxt.graphicsToSequence(clip.getMaxX()));;
      int[] baseCalls = trace.getBasecalls();
      int[] traceA = trace.getTrace(DNATools.a());
      int[] traceG = trace.getTrace(DNATools.g());
      int[] traceC = trace.getTrace(DNATools.c());
      int[] traceT = trace.getTrace(DNATools.t());

      g.setColor(Color.green);
      renderTrace(baseCalls, traceA, g, ctxt, min, max);
      g.setColor(Color.black);
      renderTrace(baseCalls, traceG, g, ctxt, min, max);
      g.setColor(Color.blue);
      renderTrace(baseCalls, traceC, g, ctxt, min, max);
      g.setColor(Color.red);
      renderTrace(baseCalls, traceT, g, ctxt, min, max);
    } catch (IllegalSymbolException ise) {
      throw new BioError("Can't process trace file", ise);
    }
  }

  private void renderTrace(int[] baseCalls, int[] trace, Graphics2D g, SequenceRenderContext ctxt, int min, int max) {
    double depthScale = depth / 2000.0; // assume X gredations
    Line2D line = new Line2D.Float();
    for(int pos = min; pos <= max; pos++) {
      int minT;
      int maxT;

      if(pos == 1) {
        minT = 0;
      } else {
        minT = (baseCalls[pos - 2] + baseCalls[pos - 1]) / 2;
      }

      if(pos == baseCalls.length) {
        maxT = trace.length - 1;
      } else {
        maxT = (baseCalls[pos - 1] + baseCalls[pos - 0]) / 2;
      }

      double scale = ctxt.getScale() / (double) (maxT - minT);

      double stg = ctxt.sequenceToGraphics(pos);
      for(int i = minT; i < maxT; i++) {
        int j = i - minT;
        line.setLine(
          stg + scale * (0.5 + j),
          depth - trace[i] * depthScale,
          stg + scale * (0.5 + j + 1),
          depth - trace[i + 1] * depthScale
        );

        g.draw(line);
      }
    }
  }

  public void setTrace(ABITrace trace)
  throws ChangeVetoException {
    ChangeSupport cs = getChangeSupport(TRACE);
    synchronized(cs) {
      ChangeEvent ce = new ChangeEvent(this, TRACE, trace, this.trace);
      cs.firePreChangeEvent(ce);
      this.trace = trace;
      cs.firePostChangeEvent(ce);
    }
  }

  public ABITrace getTrace() {
    return trace;
  }

  public void setDepth(double depth)
  throws ChangeVetoException {
    if(depth < 0.0) {
      throw new ChangeVetoException("Can't set depth to a negative number: " + depth);
    }

    ChangeSupport cs = getChangeSupport(DEPTH);
    synchronized(cs) {
      ChangeEvent ce = new ChangeEvent(this, DEPTH, new Double(depth), new Double(this.depth));
      cs.firePreChangeEvent(ce);
      this.depth = depth;
      cs.firePostChangeEvent(ce);
    }
  }

  public double getDepth(SequenceRenderContext src) {
    if(src.getDirection() == SequenceRenderContext.HORIZONTAL && trace != null) {
      return depth;
    } else {
      return 0.0;
    }
  }

  public double getMinimumLeader(SequenceRenderContext src) {
    return 0.0;
  }

  public double getMinimumTrailer(SequenceRenderContext src) {
    return 0.0;
  }

  public SequenceViewerEvent processMouseEvent(
    SequenceRenderContext src,
    MouseEvent me,
    List path
  ) {
    // don't do anything
    return null;
  }
}

