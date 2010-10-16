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
 * An implementation of SequenceRenderer that delegates rendering to another
 * renderer.
 * 
 * <p>
 * This takes care of all event notification and method invocation for you.
 * Subclass this and over-ride methods, and then possibly call the method on
 * super to forward the call on to the wrapped renderer.
 * </p>
 *
 * @author Matthew Pocock
 */
public class SequenceRendererWrapper
extends AbstractChangeable
implements SequenceRenderer, Serializable {
  public static ChangeType RENDERER = new ChangeType(
    "The renderer used to render the filtered features has changed",
    "org.biojava.bio.gui.sequence.FilteringRenderer",
    "RENDERER",
    SequenceRenderContext.LAYOUT
  );
  
  private SequenceRenderer renderer;
  private transient ChangeForwarder rendForwarder;

  /**
   *  Create a new renderer with no wrapped renderer.
   *
   * It is important that you call setRenderer() on this instance before
   * trying to use it to display anything.
   */
  public SequenceRendererWrapper() {}

  /**
   * Create a new wrapper with a wrapped renderer
   *
   * @param renderer  the SequenceRenderer to wrap up
   */
  public SequenceRendererWrapper(SequenceRenderer renderer) {
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
          c.addChangeListener(rendForwarder);
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

