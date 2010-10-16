package org.biojava.bio.gui.sequence;

import java.awt.Graphics2D;
import java.awt.Shape;
import java.awt.geom.AffineTransform;
import java.util.Iterator;

import org.biojava.bio.seq.Feature;

/**
 *
 *
 * @author Matthew Pocock
 */
public class CircularFeaturesRenderer
implements CircularRenderer {
  private CircularFeatureRenderer renderer;

  public CircularFeaturesRenderer() {}

  public CircularFeaturesRenderer(CircularFeatureRenderer renderer) {
    this.renderer = renderer;
  }

  public CircularFeatureRenderer getRenderer() {
    return renderer;
  }

  public void setRenderer(CircularFeatureRenderer renderer) {
    this.renderer = renderer;
  }

  public double getDepth(CircularRendererContext crc) {
    return renderer.getDepth(crc);
  }

  public void paint(Graphics2D g2, CircularRendererContext crc) {
    for (Iterator i = crc.getFeatures().features(); i.hasNext();) {
      Shape clip = g2.getClip();
      AffineTransform at = g2.getTransform();

      Feature f = (Feature) i.next();
      renderer.renderFeature(g2, f, crc);

      g2.setTransform(at);
      g2.setClip(clip);
    }
  }
}
