package org.biojava.bio.gui.sequence;

import java.awt.Graphics2D;

import org.biojava.bio.seq.Feature;

/**
 *
 *
 * @author Matthew Pocock
 */
public interface CircularFeatureRenderer {
  public void renderFeature(
    Graphics2D g,
    Feature f,
    CircularRendererContext context );

  public double getDepth(CircularRendererContext crc);
}
