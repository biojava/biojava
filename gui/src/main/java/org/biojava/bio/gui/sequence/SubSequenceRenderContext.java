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

import java.awt.Font;
import java.awt.geom.Point2D;

import org.biojava.bio.seq.FeatureHolder;
import org.biojava.bio.symbol.RangeLocation;
import org.biojava.bio.symbol.SymbolList;

/**
 * Allows a new renderer to "wrap" another one, replacing one or more values.
 *
 * <p>
 * Use this when implementing SequenceRenderer classes that modify the data
 * that is passed on to delegate renderers e.g. filtering the features,
 * transforming the sequence or moving the rendering co-ordinates.
 * </p>
 *
 * @author Matthew Pocock
 */
public class SubSequenceRenderContext
implements SequenceRenderContext {
  private final SequenceRenderContext src;
  private final SymbolList symbols;
  private final FeatureHolder features;
  private final RangeLocation range;
  private final int symOffset;

  public SubSequenceRenderContext(
    SequenceRenderContext src,
    SymbolList symbols,
    FeatureHolder features,
    RangeLocation range
  ) {
    this(src, symbols, features, range, 0);
  }

  public SubSequenceRenderContext(
          SequenceRenderContext src,
          SymbolList symbols,
          FeatureHolder features,
          RangeLocation range,
          int symOffset
  ) {
    this.src = src;
    this.symbols = symbols;
    this.features = features;
    this.range = range;
    this.symOffset = symOffset;
  }

  public int getDirection() {
    return src.getDirection();
  }

  public double getScale() {
    return src.getScale();
  }

  public double sequenceToGraphics(int i) {
    return src.sequenceToGraphics(i + symOffset);
  }

  public int graphicsToSequence(double d) {
    return src.graphicsToSequence(d) - symOffset;
  }

  public int graphicsToSequence(Point2D point) {
    return src.graphicsToSequence(point) - symOffset;
  }

  public SymbolList getSymbols() {
    if(symbols == null) {
      return src.getSymbols();
    } else {
      return symbols;
    }
  }

  public FeatureHolder getFeatures() {
    if(features == null) {
      return src.getFeatures();
    } else {
      return features;
    }
  }

  public RangeLocation getRange() {
    if(range == null) {
      return src.getRange();
    } else {
      return range;
    }
  }

  public SequenceRenderContext.Border getLeadingBorder() {
    return src.getLeadingBorder();
  }

  public SequenceRenderContext.Border getTrailingBorder() {
    return src.getTrailingBorder();
  }

  public Font getFont() {
    return src.getFont();
  }
}
