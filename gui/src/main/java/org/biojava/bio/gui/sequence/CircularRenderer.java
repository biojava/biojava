package org.biojava.bio.gui.sequence;

import java.awt.Graphics2D;

/**
 * Render information from a CircularRendererContext onto a graphics context.
 *
 * <p>
 * Every CircularRenderer paints information about a sequence, it's symbols or features
 * or some other property, into a ring. The depth of the ring is given by getDepth().
 * </p>
 *
 * @author Matthew Pocock
 * @since 1.4
 */
public interface CircularRenderer {
  /**
   * Get the depth needed for this renderer.
   *
   * @param crc   the CircularRendererContext to render information from
   * @return      the depth required to render the context
   */
  public double getDepth(CircularRendererContext crc);

  /**
   * Paint this renderer.
   *
   * @param g2    the graphics to paint to
   * @param crc   the context giving the data to paint
   */
  public void paint(Graphics2D g2, CircularRendererContext crc);
}
