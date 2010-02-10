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
import java.awt.Shape;
import java.awt.event.MouseEvent;
import java.awt.geom.AffineTransform;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.FeatureHolder;
import org.biojava.utils.AbstractChangeable;
import org.biojava.utils.ChangeType;

/**
 * Allows you to stack multiple feature renderers up (for example a label renderer and
 * a beaded renderer) and have them treated as a single renderer for layout.
 *
 * @author Matthew Pocock
 */
public class StackedFeatureRenderer
extends AbstractChangeable
implements FeatureRenderer {
  public static final ChangeType RENDERERS = new ChangeType(
    "The renderers have changed",
    StackedFeatureRenderer.class,
    "RENDERERS",
    SequenceRenderContext.LAYOUT
  );
  
  private List renderers = new ArrayList();
  
  public StackedFeatureRenderer() {}
  
  // fixme: events?
  public void addRenderer(FeatureRenderer renderer) {
    renderers.add(renderer);
  }
  
  public void removeRenderer(FeatureRenderer renderer) {
    renderers.remove(renderer);
  }
  
  public double getDepth(SequenceRenderContext src) {
    double depth = 0.0;
    
    for(Iterator i = renderers.iterator(); i.hasNext(); ) {
      FeatureRenderer fr = (FeatureRenderer) i.next();
      depth += fr.getDepth(src);
    }
    
    return depth;
  }
  
  public void renderFeature(
    Graphics2D g,
    Feature f, 
    SequenceRenderContext src
  ) {
    AffineTransform at = g.getTransform();
    Shape oldClip = g.getClip();
    
    for(Iterator i = renderers.iterator(); i.hasNext(); ) {
      FeatureRenderer fr = (FeatureRenderer) i.next();
      fr.renderFeature(g, f, src);

      double depth = fr.getDepth(src);
      g.translate(0.0, depth);
    }
    
    g.setClip(oldClip);
    g.setTransform(at);
  }
  
  public FeatureHolder processMouseEvent(
    FeatureHolder hits,
    SequenceRenderContext src,
    MouseEvent me
  ) {
    return hits;
  }
}
