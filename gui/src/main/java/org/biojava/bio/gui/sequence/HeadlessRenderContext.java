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
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.symbol.RangeLocation;
import org.biojava.bio.symbol.SymbolList;

/**
 * <p>
 * A stand-alone SequenceRenderContext to make it easy to render to an image.
 * </p>
 *
 * <p>
 * This class makes it very easy to render sequence information into an
 * arbitrary graphics object without the need to fuss about with AWT or
 * Swing components. You chose the width of the image and the region of the
 * sequence to render. It will calculate the scale factor to ensure that the
 * whole region of the sequence fits into that width. You can then use the
 * context to render any number of SequenceRenderer instances to any Graphics2D
 * instance you want, for example, to an image that's to be written out by a
 * servlet.
 * </p>
 *
 * <h2>Example</h2>
 *
 * <pre>
 * HeadlessRenderContext ctxt = new HeadlessRenderContext(
 *   seq,   // the sequence to render
 *   range, // a RangeLocation giving the block you want to render
 *   width  // an int specifying the image width in pixles
 * );
 *
 * BufferedImage img = new BufferedImage(
 *   width,                                   // image width
 *   (int) Math.ceil(seqRend.getDepth(ctxt),  // calculated height
 *   BufferedImage.TYPE_INT_RGB               // let's use RGB
 * );
 *
 * // set stuff up
 * Graphics2D graph = img.createGraphics();
 * graph.setPaint(Color.WHITE);
 * graph.fillRect(0, 0, img.getWidth(), img.getHeight());
 *
 * // and now render the sequences
 * sequenceRenderer.paint(graph, ctxt);
 *
 * // let's dump this out as a png
 * ImageIO.write(image, "png", myFile);
 * </pre>
 *
 * @since 1.3
 * @author Matthew Pocock
 */

public class HeadlessRenderContext
implements SequenceRenderContext {
  private static final Font FONT = new Font(null, Font.PLAIN, 10);

  private final RangeLocation range;
  private final double scale;
  private final Sequence seq;
  private final double offset;

  public HeadlessRenderContext(Sequence seq, RangeLocation range, int width) {
    this.seq = seq;
    this.range = range;
    this.scale = (double) width /
    (double) (range.getMax() - range.getMin() + 1);
    offset = -( (double) range.getMin() * scale);
  }

  public int getDirection() {
    return HORIZONTAL;
  }

  public FeatureHolder getFeatures() {
    return seq;
  }

  public Font getFont() {
    return FONT;
  }

  public SequenceRenderContext.Border getLeadingBorder() {
    return new SequenceRenderContext.Border();
  }

  public RangeLocation getRange() {
    return range;
  }

  public double getScale() {
    return scale;
  }

  public SymbolList getSymbols() {
    return seq;
  }

  public SequenceRenderContext.Border getTrailingBorder() {
    return new SequenceRenderContext.Border();
  }

  public double sequenceToGraphics(int i) {
    return ((double) (i - 1)) * scale + offset;
  }

  public int graphicsToSequence(Point2D point) {
    return graphicsToSequence(point.getX());
  }

  public int graphicsToSequence(double d) {
    return ((int) ((d - offset) / scale)) + 1;
  }
}
