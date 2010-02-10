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
import java.util.Iterator;
import java.util.List;

import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.symbol.RangeLocation;
import org.biojava.bio.symbol.SymbolList;

// The graphics model in Java
// drawing space -> applet space -> device space
// All operations are cumulative, including translates
// translates will move drawing rightward/downward for any supplied value

/**
 * Compute sites of stop codons.  It uses a child renderer to
 * do actual drawing.
 *
 * @author David Huen
 */ 


public class StopRenderer implements SequenceRenderer {
    private double scaleThreshold = 0.005;
    private SixFrameRenderer pane;
    private int moduloFrame;
    private StrandedFeature.Strand strand;

    public StopRenderer(
             SixFrameRenderer pane,
             int moduloFrame, 
             StrandedFeature.Strand strand) {
    this.pane = pane;
    this.moduloFrame = moduloFrame;
    this.strand = strand;
    }

    public double getDepth(SequenceRenderContext src) {
      // an arbitrary limit is set here to prevent excessive sequence
      // download.
      if (src.getScale() < scaleThreshold) return 0;
        else return pane.getDepth(src);
    }

    public double getMinimumLeader(SequenceRenderContext src) {
      return 0.0;
    }

    public double getMinimumTrailer(SequenceRenderContext src) {
      return 0.0;
    }

    private boolean isStop(SymbolList seq,
      int base,
      StrandedFeature.Strand strand) {
      // tests whether there is a stop at given location.
      // the triplet is either base, +1, +2 or -1, -2
      // depending on the strand searched
      if (strand == StrandedFeature.POSITIVE) {
        // check that search does not exceed bounds
        if (base + 2 > seq.length()) return false;

        // search top strand
        // first base must be t
        if (seq.symbolAt(base) != DNATools.t()) return false;

        // second base cannot be c or t
        if (seq.symbolAt(base+1) == DNATools.c()) return false;
        if (seq.symbolAt(base+1) == DNATools.t()) return false;

        // if second base is g, the third must be a
        if (seq.symbolAt(base+1) == DNATools.g()) {
          if (seq.symbolAt(base+2) != DNATools.a()) return false;
        }
        else {
          // second base is a: third must be a or g.
          if (seq.symbolAt(base+2) == DNATools.c()) return false;
          if (seq.symbolAt(base+2) == DNATools.t()) return false;
        }

        // oh well, must be a stop, innit?
        return true;

      } else {
        // check bounds
        if (base - 2 < 1) return false;

        // search bottom strand
        // first base must be t
        if (seq.symbolAt(base) != DNATools.a()) return false;

        // second base cannot be c or t on reverse strand
        if (seq.symbolAt(base-1) == DNATools.a()) return false;
        if (seq.symbolAt(base-1) == DNATools.g()) return false;

        // if second base is g, the third must be a
        if (seq.symbolAt(base-1) == DNATools.c()) {
          if (seq.symbolAt(base-2) != DNATools.t()) return false;
        }
        else {
          // second base is a: third must be a or g.
          if (seq.symbolAt(base-2) == DNATools.a()) return false;
          if (seq.symbolAt(base-2) == DNATools.g()) return false;
        }

        // ach! a stop!
        return true;
      }
    }

    private void renderOneFrame(
      Graphics2D g,
      SequenceRenderContext src,
      RangeLocation range,
      boolean onceOnly) {
      // method to draw by checking succeeding triplets for
      // stop codons.
      // write it for horizontal rendering first.
      SymbolList seq = src.getSymbols();

      // get extent of sequence to render
      // hope it agrees with clip region!
      int minS = range.getMin();
      int maxS = range.getMax();

      // we start at the first triplet whose first base is within
      // the range.
      if (minS%3 > moduloFrame) {
        // first triplet of my frame is in next mod-zero triplet
        minS = (minS/3 + 1) * 3 + moduloFrame;
      }
      else if (minS%3 != moduloFrame) {
        // first triplet is in current mod-zero triplet
        minS = (minS/3) * 3 + moduloFrame;
      }

      // now we search every triplet from minS upward seeking stops.
      for (int base = minS; base <= maxS; base += 3) {
        // check for stop
        if (!isStop(seq, base, strand)) continue;

        // we have a stop, render a line
        pane.drawLine(g, src, base, strand);

        // do I call it quits now?
        if (onceOnly) return;
      }
    }

    public void paint(
      Graphics2D g,
      SequenceRenderContext src
    ) {
      double scale = src.getScale();

      // this is a completely arbitrary limit to stop my viewer
      // from attempting to trigger the download of HUGE amounts 
      // of sequence.
      if (scale < scaleThreshold) return;

      // could we get more than one stop per pixel at the current
      // scale?      
      if (scale < 0.05) {
        // yes, we can. Iterate thru' graphics space.
        Iterator extentsI = pane.sequenceExtentOfPixels(src).iterator();

        // check each extent for stops
        while (extentsI.hasNext()) {
          RangeLocation range = (RangeLocation) extentsI.next();
          renderOneFrame(g, src, range, true);
        }
      }
      else {
        // no we can't. Iterate thru' sequence.
        renderOneFrame(g, src, src.getRange(), false);
      }
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
