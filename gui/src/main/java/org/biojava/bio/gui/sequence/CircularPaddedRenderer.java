package org.biojava.bio.gui.sequence;

import java.awt.Graphics2D;

/**
 *
 *
 * @author Matthew Pocock
 */
public class CircularPaddedRenderer
implements CircularRenderer {
  private CircularRenderer delegate;
  private double prePadding;
  private double postPadding;

  public CircularPaddedRenderer() {
    this(null, 0.0, 0.0);
  }

  public CircularPaddedRenderer(CircularRenderer delegate) {
    this(delegate, 0.0, 0.0);
  }

  public CircularPaddedRenderer(double prePadding, double postPadding) {
    this(null, prePadding, postPadding);
  }

  public CircularPaddedRenderer(CircularRenderer delegate,
                                double prePadding,
                                double postPadding)
  {
    this.delegate = delegate;
    this.prePadding = prePadding;
    this.postPadding = postPadding;
  }

  public double getDepth(CircularRendererContext crc) {
    CircularRendererContext subCtxt = new SubCircularRendererContext(
            crc,
            null,
            null,
            crc.getRadius() + prePadding);
    return prePadding + delegate.getDepth(subCtxt) + postPadding;
  }

  public void paint(Graphics2D g2, CircularRendererContext crc) {
    CircularRendererContext subCtxt = new SubCircularRendererContext(
            crc,
            null,
            null,
            crc.getRadius() + prePadding);
    delegate.paint(g2, subCtxt);
  }
}
