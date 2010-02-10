package org.biojava.bio.gui.sequence;

import java.awt.Graphics2D;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;


/**
 * Renders multiple renderers, each in their own concentric rings.
 *
 * @author Matthew Pocock
 * @since 1.4
 */
public class CircularMLR
implements CircularRenderer {
  private List renderers = new ArrayList();

  public void addRenderer(CircularRenderer renderer) {
    renderers.add(renderer);
  }

  public void removeRenderer(CircularRenderer renderer) {
    renderers.remove(renderer);
  }

  public double getDepth(CircularRendererContext crc) {
    double depth = 0;

    for(Iterator i = renderers.iterator(); i.hasNext(); ) {
      CircularRenderer rend = (CircularRenderer) i.next();
      CircularRendererContext subCtxt = new SubCircularRendererContext(
              crc,
              null,
              null,
              crc.getRadius() + depth);
      depth += rend.getDepth(subCtxt);
    }

    return depth;
  }

  public void paint(Graphics2D g2, CircularRendererContext crc) {
    double depth = 0.0;

    for(Iterator i = renderers.iterator(); i.hasNext(); ) {
      CircularRenderer rend = (CircularRenderer) i.next();
      CircularRendererContext subCtxt = new SubCircularRendererContext(
              crc,
              null,
              null,
              crc.getRadius() + depth);
      rend.paint(g2, subCtxt);
      depth += rend.getDepth(crc);
    }
  }
}
