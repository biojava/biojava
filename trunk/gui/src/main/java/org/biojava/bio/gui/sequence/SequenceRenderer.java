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
import java.awt.event.MouseEvent;
import java.util.List;

import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeForwarder;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeType;

/**
 * The interface for things that can render a line of information about a
 * sequence.
 * <p>
 * Renderers are always activated within the context of a particular sequence
 * panel. A single Renderer can be shared among many sequence panels, or added
 * multiple times to the same panel. The renderer is required to request how
 * much leading and trailing space it requires, as well as the depth (space
 * orthogonal to the direction in which the sequence is rendered).
 * <p>
 * The leading and trailing distances are the number of pixels overhang needed
 * to cleanly render any line of sequence information. For example, a ruler will
 * need trailing space to render the total sequence length at the end.
 *
 * @author Thomas Down
 * @author Matthew Pocock
 * @author Keith James
 */
public interface SequenceRenderer {
  
  /**
   * Render a portion (possibly all) of the information for src to g, displaying
   * all of the data that would fall within seqBox.
   *
   * @param g the Graphics2D to render to
   * @param src the SequenceRenderContext that encapsulates the information to render
   */
  void paint(Graphics2D g, SequenceRenderContext src);
  
  /**
   * Retrieve the depth of this renderer when rendering src.
   * <p>
   * The depth may vary between sequence panels - for example based upon
   * sequence length. Each line of information in the SequenceRendererContext
   * only renders a region of the sequence. The depth for one complete line may
   * be different from that for another due to the sequence having more or less
   * information in that region to show. For example, a feature renderer
   * implementation may chose to collapse down to a depth of zero pixels if
   * there are no features to render within a region.
   *
   * @param src the SequenceRenderContext to return info for
   * @return the depth of the renderer for that sequence panel
   */
  double getDepth(SequenceRenderContext src);

  /**
   * Retrieve the minimum leading distance for this renderer when rendering src.
   * <p>
   * The leading distance may vary between sequence panels - for example based
   * upon sequence length.
   *
   * @param src the SequenceRenderContext to return info for
   * @return the leading distance of the renderer for that sequence panel
   */
  double getMinimumLeader(SequenceRenderContext src);
    
  /**
   * Retrieve the minimum trailing distance for this renderer when rendering src.
   * <p>
   * The trailing distance may vary between sequence panels - for example based
   * upon sequence length.
   *
   * @param src the SequenceRenderContext to return info for
   * @return the trailing distance of the renderer for that sequence panel
   */
  double getMinimumTrailer(SequenceRenderContext src);
  
  /**
   * Produce a SequenceViewerEvent in response to a mouse gesture.
   * <p>
   * A SequenceRenderer that performs any form of coordinate remapping should
   * ensure that it appropriately transforms the mouse event. However, in the
   * SequenceViewerEvent returned, the MouseEvent should be in untransformed
   * coordinates.
   * <p>
   * The SequenceRenderer implementation should append itself to the path list
   * before constructing the SequenceViewerEvent.
   *
   * @param src the SequenceRenderContext currently in scope
   * @param me  a MouseEvent that caused this request
   * @param path the List of SequenceRenderer instances passed through so far
   * @return a SequenceViewerEvent encapsulating the mouse gesture
   *
   * @since 1.2
   */
  SequenceViewerEvent processMouseEvent(
    SequenceRenderContext src,
    MouseEvent me,
    List path
  );
  
  public static class RendererForwarder extends ChangeForwarder {
    public RendererForwarder(SequenceRenderer source, ChangeSupport cs) {
      super(source, cs);
    }
    
    public ChangeEvent generateEvent(ChangeEvent ce) {
      ChangeType cType = ce.getType();
      ChangeType newType;
      if(cType.isMatchingType(SequenceRenderContext.LAYOUT)) {
        newType = SequenceRenderContext.LAYOUT;
      } else if(cType.isMatchingType(SequenceRenderContext.REPAINT)) {
        newType = SequenceRenderContext.REPAINT;
      } else {
        return null;
      }
      return new ChangeEvent(getSource(), newType, null, null, ce);
    }
  }
}
