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

/**
 * The interface for things that can render labels for a line of information
 * about a sequence.
 * <p>
 * Renderers are always activated within the context of a particular
 * SequenceRenderContext.
 * A single LabelRenderer can be shared among many sequence panels, or added
 * multiple times to the same panel. The renderer is required to request how
 * much leading and trailing space it requires, as well as the depth (space
 * orthogonal to the direction that the sequence is rendered).
 *
 * @author Matthew Pocock
 */
public interface LabelRenderer {
  
  /**
   * Render a label for the information for sp to g.
   *
   * @param g the Graphics2D to render to
   * @param sp the SequencePanel that encapsulates the information to render
   * @param min the minimum symbol to render (inclusive)
   * @param max the maximum symbol to render (inclusive)
   */
  void paint(
    Graphics2D g, SequenceRenderContext sp,
    int min, int max, SequenceRenderContext.Border border
  );

  /**
   * Retrieve the minimum space required to render the label.
   *
   * @param sp the SequencePanel to return info for
   * @return the leading distance of the renderer for that sequence panel
   */
  double getMinimumWidth(SequenceRenderContext sp);
    
  public static LabelRenderer RENDER_NOTHING = new RenderNothing();
  
  static class RenderNothing implements LabelRenderer {
    public void paint(
         Graphics2D g, SequenceRenderContext sp,
         int min, int max, SequenceRenderContext.Border border
    ) {}
    
    public double getMinimumWidth(SequenceRenderContext src) {
      return 0.0;
    }
  };
}
