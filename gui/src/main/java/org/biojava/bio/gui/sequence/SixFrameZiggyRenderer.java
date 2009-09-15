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
import java.awt.Paint;
import java.awt.event.MouseEvent;
import java.util.Arrays;
import java.util.Iterator;

import org.biojava.bio.Annotation;
import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.bio.seq.ComponentFeature;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.FeatureFilter;
import org.biojava.bio.seq.FeatureHolder;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SimpleAssembly;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.RangeLocation;
import org.biojava.utils.AbstractChangeable;
import org.biojava.utils.ChangeVetoException;

/**
 * A feature renderer that computes the data necessary to render
 * multi-exon transcripts without CDS data.
 * <p>
 * The actual drawing is done by a child renderer.  In this case,
 * SixFrameRenderer is used, which can use data from this renderer
 * to display transcripts in the correct translation frames.
 *
 * @author David Huen
 */

public class SixFrameZiggyRenderer
  extends AbstractChangeable
  implements FeatureRenderer, java.io.Serializable {
  private SixFrameRenderer pane;
  public SixFrameZiggyRenderer(SixFrameRenderer pane) {
    this.pane = pane;
  }

  public void setFill(Paint p)
    throws ChangeVetoException {
    pane.setFill(p);
  }

  public Paint getFill() {
    return pane.getFill();
  }

  public void setOutline(Paint p)
    throws ChangeVetoException {
    pane.setOutline(p);
  }

  public Paint getOutline() {
    return pane.getOutline();
  }

  public void setBlockDepth(double depth)
    throws ChangeVetoException {
    pane.setBlockWidth(depth);
  }

  public double getBlockDepth() {
    return pane.getBlockWidth();
  }

  public double getDepth(SequenceRenderContext src) {
    return pane.getDepth(src);
  }

  private boolean isStop(
                    Sequence seq,
                    int base,
                    StrandedFeature.Strand strand) {
    // tests whether there is a stop at given location.
    // the triplet is either base, +1, +2 or -1, -2
    // depending on the strand searched
    if (strand == StrandedFeature.POSITIVE) {
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

  private int findORF(
                Sequence seq,
                StrandedFeature.Strand strand) {
    // finds in a SymbolList the specified phase with
    // longest ORF and returns the phase.
    int[] lastStop = {0, 0, 0};
//    int[] longestORF = {0, 0, 0};
    int bestPhase = 0;
    int highestORFSize = 0;

    // scan thru' the sequence looking for stops
    int length = seq.length();
    if (length < 4) return 0;

    // set limits of search
    int startSearch, endSearch;
    if (strand == StrandedFeature.POSITIVE) {
      startSearch = 1;
      endSearch = length - 2;
    }
    else {
      startSearch = 3;
      endSearch = length;
    }

    for (int i=startSearch; i <= endSearch; i++) {
      if (isStop(seq, i, strand)) {
        // stop found
        int phase = i%3;
        int currORFSize = i - lastStop[phase];

        // is this a candidate for best phase?
        if (currORFSize > highestORFSize) {
          bestPhase = phase;
          highestORFSize= currORFSize;
//          longestORF[phase] = currORFSize;
        }
        lastStop[phase] = i;
//        System.out.println("findORF i phase, largest: " + i + " "
//           + phase + " " + currORFSize);
      }
    }

    // there is always the possibility that there are a few stops
    // near the beginning then no more.
    // The best phase will then be misdetected.
    // Assume closure at end of frame.
    for (int i=0; i < 3; i++) {
      int currORFSize = endSearch - lastStop[i];
      if (currORFSize > highestORFSize) {
        bestPhase = i;
        highestORFSize= currORFSize;
//        longestORF[phase] = currORFSize;
      }
    }

    return bestPhase;
  }

  private Sequence assembleFusedSequence(Feature [] block, Sequence seq) {
    // assembles a fused sequence from component features
    // only assembles in the forward direction but will
    // sort exons as necessary.

    SimpleAssembly sa = new SimpleAssembly("temp", "temp");
    ComponentFeature.Template cft = new ComponentFeature.Template();
    cft.annotation = Annotation.EMPTY_ANNOTATION;
    cft.strand = StrandedFeature.POSITIVE;
    cft.componentSequence = seq;

    int last = 0;
    for (int j= 0; j < block.length; j++) {
      // fuse all "exons" irrespective of orientation.
      Feature thisExon = block[j];

      cft.componentLocation = thisExon.getLocation();
      int length = cft.componentLocation.getMax() -
                     cft.componentLocation.getMin() + 1;
      cft.location = new RangeLocation(last+1, last+length);
      last += length;
//      System.out.println("assemble: " + cft.componentLocation.getMin() + " " + cft.componentLocation.getMax());

      try {
        sa.createFeature(cft);
      }
      catch (BioException be) {
        throw new BioError(
          "Couldn't merge exons.", be
          );
      }
      catch (ChangeVetoException cve) {
        throw new BioError(

          "Couldn't merge exons.",cve
          );
      }
    }
    return sa;
  }

  public void renderFeature(
                Graphics2D g,
                Feature f,
                SequenceRenderContext context) {
//    System.out.println("SixFrameZiggyRenderer called");

    if (!(f instanceof StrandedFeature)) return;
    // create a fused version of the transcript
    // this solution is ugly as hell, a botched abortion of a fix
    //
    // the algorithm is hideously simple.  Irrespective of the
    // strandedness of the transcript, a fused sequence will be
    // generated in the forward direction.
    //
    // this "transcript" will then be searched for the longest
    // ORF in the correct strand and the phase of the largest ORF
    // returned.  It really doesn't matter whether the min sequence
    // end is the 5' or 3' of the transcript as phase is consistent
    // thru' an ORF.
    //
    // By just passing the best phase over to SixFrameRenderer, the
    // the phase of successive exons can be computed from just the
    // previous exon phase and the preceding intron size.

    //filter for only the exons
    FeatureFilter filt = new FeatureFilter.ByType("exon");
    FeatureHolder exons = f.filter(filt, false);

    // sort the returned exons in ascending order
    // disappointment...
    int featureCount = exons.countFeatures();
    Feature[] orderedExons = new Feature[featureCount];
    int i=0;
    for (Iterator fi=exons.features(); fi.hasNext();) {
      orderedExons[i++] = (Feature) fi.next();
    }
    Arrays.sort(orderedExons, new Feature.ByLocationComparator());

    Sequence fused = assembleFusedSequence(orderedExons, f.getSequence());

    StrandedFeature.Strand strand = ((StrandedFeature) f).getStrand();

    // findORF will find the best phase within the "ORF" but that
    //  needs to be corrected for the phase in which the ORF is
    // embedded into the sequence
    int phase = findORF(fused, strand);

//    System.out.println("fused length, phase, strand: " + fused.length() + " "
//                        + phase + " " + strand);
//    System.out.println("sequence is :- " + fused.seqString());

    // Iterate over exon child features: these are already ordered.
    Location loc = null;
    for (i = 0; i < orderedExons.length; i++) {
      loc = ((Feature) orderedExons[i]).getLocation();
      if (i == 0) {
        // first exon
        pane.startZiggy(strand, (2 + loc.getMin() + phase)%3);
        pane.renderLocation(g, context, loc);
//        System.out.println("block value is " + loc);
      }
      else {
        pane.renderLocation(g, context, loc);
//        System.out.println("block value is " + loc);
      }
    }

  }

  public FeatureHolder processMouseEvent(
    FeatureHolder hits,
    SequenceRenderContext src,
    MouseEvent me
  ) {
    return hits;
  }
}
