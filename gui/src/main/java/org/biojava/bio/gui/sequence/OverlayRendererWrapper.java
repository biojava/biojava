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
import java.io.Serializable;
import java.util.List;

import org.biojava.utils.AbstractChangeable;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeForwarder;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.Changeable;

/**
 * This class wraps SequenceRenderer classes to present
 * an Overlay Marker interface that will indicate to 
 * LayeredRenderer that subsequent renderers should overdraw
 * the same space as the wrapped renderer.
 * <p>
 * Extensively plagiarised from SequenceRendererWrapper.java.
 *
 * @author     David Huen
 */
public class OverlayRendererWrapper
extends AbstractChangeable
implements SequenceRenderer, Serializable, OverlayMarker {
  public static ChangeType RENDERER = new ChangeType(
    "The renderer used to render the filtered features has changed",
    "org.biojava.bio.gui.sequence.OverlayRendererWrapper",
    "RENDERER",
    SequenceRenderContext.LAYOUT
  );
  
  private SequenceRenderer renderer;
  private transient ChangeForwarder rendForwarder;
  
  public OverlayRendererWrapper() {}
  public OverlayRendererWrapper(SequenceRenderer renderer) {
    this.renderer = renderer;
  }
  
  protected ChangeSupport getChangeSupport(ChangeType ct) {
    ChangeSupport cs = super.getChangeSupport(ct);
    
    if(rendForwarder == null) {
      rendForwarder = new SequenceRenderer.RendererForwarder(this, cs);
      if((renderer != null) && (renderer instanceof Changeable)) {
        Changeable c = (Changeable) this.renderer;
        c.addChangeListener(
          rendForwarder,
          SequenceRenderContext.REPAINT
        );
      }
    }
    
    return cs;
  }
  
  public void setRenderer(SequenceRenderer renderer)
  throws ChangeVetoException {
    if(hasListeners()) {
      ChangeEvent ce = new ChangeEvent(
        this, RENDERER,
        renderer, this.renderer
      );
      ChangeSupport cs = getChangeSupport(RENDERER);
      synchronized(cs) {
        cs.firePreChangeEvent(ce);
        if((this.renderer != null) && (this.renderer instanceof Changeable)) {
          Changeable c = (Changeable) this.renderer;
          c.removeChangeListener(rendForwarder);
        }
        this.renderer = renderer;
        if(renderer instanceof Changeable) {
          Changeable c = (Changeable) renderer;
          c.removeChangeListener(rendForwarder);
        }
        cs.firePostChangeEvent(ce);
      }
    } else {
      this.renderer = renderer;
    }
  }
  
  public SequenceRenderer getRenderer() {
    return this.renderer;
  }
  
  public double getDepth(SequenceRenderContext src) {
    SequenceRenderer sr = getRenderer();
    if(sr == null) {
      return 0.0;
    } else {
      return sr.getDepth(src);
    }
  }
  
  public double getMinimumLeader(SequenceRenderContext src) {
    SequenceRenderer sr = getRenderer();
    if(sr == null) {
      return 0.0;
    } else {
      return sr.getMinimumLeader(src);
    }
  }
  
  public double getMinimumTrailer(SequenceRenderContext src) {
    SequenceRenderer sr = getRenderer();
    if(sr == null) {
      return 0.0;
    } else {
      return sr.getMinimumTrailer(src);
    }
  }
  
  public void paint(
    Graphics2D g,
    SequenceRenderContext src
  ) {
    SequenceRenderer sr = getRenderer();
    if(sr != null) {
      sr.paint(g, src);
    }
  }
  
  public SequenceViewerEvent processMouseEvent(
    SequenceRenderContext src,
    MouseEvent me,
    List path
  ) {
    path.add(this);
    return getRenderer().processMouseEvent(src, me, path);
  }
}

