package org.biojava.bio.gui.sequence;

import java.awt.Graphics2D;

import org.biojava.bio.seq.FeatureFilter;

/**
 *
 *
 * @author Matthew Pocock
 */
public class CircularFeatureFilteringRenderer
implements CircularRenderer {
  private boolean recurse;
  private FeatureFilter filter;
  private CircularRenderer renderer;

  public CircularFeatureFilteringRenderer(CircularRenderer renderer,
                                          FeatureFilter filter,
                                          boolean recurse)
  {
    this.renderer = renderer;
    this.filter = filter;
    this.recurse = recurse;
  }

  public double getDepth(CircularRendererContext crc) {
    CircularRendererContext subCtxt = new SubCircularRendererContext(
            crc,
            null,
            crc.getFeatures().filter(filter, recurse),
            Double.NaN );
    return renderer.getDepth(subCtxt);
  }

  public void paint(Graphics2D g2, CircularRendererContext crc) {
    CircularRendererContext subCtxt = new SubCircularRendererContext(
            crc,
            null,
            crc.getFeatures().filter(filter, recurse),
            Double.NaN );
    renderer.paint(g2, subCtxt);
  }
}
