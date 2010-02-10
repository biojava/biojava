package org.biojava.bio.gui.sequence;

import java.awt.Graphics2D;
import java.util.Iterator;

import org.biojava.bio.seq.FeatureHolder;
import org.biojava.bio.symbol.GappedSymbolList;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.RangeLocation;
import org.biojava.bio.symbol.SymbolList;

/**
 * A renderer that will display a gapped sequence as a discontinuous series of
 * regions.
 *
 * <p>
 * Each ungapped block in the gapped symbol list will be displayed as a
 * contiguous region by this renderer. Where there are gaps, this renderer
 * will display nothing. Then, when the gaps are over, it will continue to
 * render the ungapped sequence. This has the effect of snapping the image
 * of the ungapped sequence where there are gaps, so as to allow it to be
 * viewed in the gapped co-ordinate system.
 * </p>
 *
 * @author Matthew Pocock
 */
public class GappedRenderer
extends SequenceRendererWrapper {
  public GappedRenderer() {
    super();
  }

  public GappedRenderer(SequenceRenderer renderer) {
    super(renderer);
  }

  public double getDepth(SequenceRenderContext src) {
    if(src.getSymbols() instanceof GappedSymbolList) {
      GappedSymbolList gsym = (GappedSymbolList) src.getSymbols();
      Location ungapped = gsym.getUngappedLocation();
      double depth = 0.0;

      Iterator bi = ungapped.blockIterator();
      while(bi.hasNext()) {
        RangeLocation loc = (RangeLocation) bi.next();
        depth = Math.max(depth, super.getDepth(makeContext(src, loc)));
      }

      return depth;
    } else {
      return super.getDepth(src);
    }
  }

  public double getMinimumLeader(SequenceRenderContext src) {
    if(src.getSymbols() instanceof GappedSymbolList) {
      GappedSymbolList gsym = (GappedSymbolList) src.getSymbols();
      Location ungapped = gsym.getUngappedLocation();
      Iterator bi = ungapped.blockIterator();
      if(bi.hasNext()) {
        return super.getMinimumLeader(
              makeContext(src, (RangeLocation) bi.next()));
      } else {
        return 0.0;
      }
    } else {
      return super.getMinimumLeader(src);
    }
  }

  public double getMinimumTrailer(SequenceRenderContext src) {
    if(src.getSymbols() instanceof GappedSymbolList) {
      GappedSymbolList gsym = (GappedSymbolList) src.getSymbols();
      Location ungapped = gsym.getUngappedLocation();
      Iterator bi = ungapped.blockIterator();
      RangeLocation loc = null;
      while(bi.hasNext()) {
        loc = (RangeLocation) bi.next();
      }
      if(loc == null) {
        return 0.0;
      } else {
        return super.getMinimumTrailer(makeContext(src, loc));
      }
    } else {
      return super.getMinimumTrailer(src);
    }
  }

  public void paint(
          Graphics2D g,
          SequenceRenderContext src
          ) {
    if(src.getSymbols() instanceof GappedSymbolList) {
      GappedSymbolList gsym = (GappedSymbolList) src.getSymbols();
      Location ungapped = gsym.getUngappedLocation();
      Iterator bi = ungapped.blockIterator();
      while(bi.hasNext()) {
        RangeLocation loc = (RangeLocation) bi.next();
        super.paint(g, makeContext(src, loc));
      }
    } else {
      super.paint(g, src);
    }
  }

  protected SequenceRenderContext makeContext(SequenceRenderContext src,
          RangeLocation loc) {
    GappedSymbolList gsl = (GappedSymbolList) src.getSymbols();
    RangeLocation sourceLoc = new RangeLocation(
            gsl.viewToSource(loc.getMin()),
            gsl.viewToSource(loc.getMax())
    );
    int trans = loc.getMin() - sourceLoc.getMin();
    SymbolList ugsl = gsl.getSourceSymbolList();

    return new SubSequenceRenderContext(
            src,
            ugsl,
            (ugsl instanceof FeatureHolder) ? (FeatureHolder) ugsl : null,
            sourceLoc,
            trans);
  }
}
