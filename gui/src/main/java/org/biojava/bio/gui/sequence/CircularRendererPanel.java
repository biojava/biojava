package org.biojava.bio.gui.sequence;

import java.awt.Graphics;
import java.awt.Graphics2D;

import javax.swing.JComponent;

import org.biojava.bio.seq.FeatureHolder;
import org.biojava.bio.symbol.SymbolList;

/**
 * Renders a sequence as a circle using a CircularRenderer.
 *
 * <p>
 * This component will first transform the graphic coordinates so that 0,0 is at
 * the centre of the circle. The size of the circle is estimated from the radius
 * property and the depth of the renderer.
 * </p>
 *
 * <p>
 * All angles are measured in radians. Some java gui classes use radians and
 * some use degrees. Be carefull to use the right one. Math has a couple
 * of methods for conversions.
 * </p>
 *
 * @author Matthew Pocock
 * @since 1.4
 */
public class CircularRendererPanel
extends JComponent {
  private final CircularRendererContext ctxt;

  {
    ctxt = new CTXT();
  }

  private SymbolList symList;
  private double radius;
  private CircularRenderer renderer;
  private double offset;

  public double getRadius() {
    return radius;
  }

  public void setRadius(double radius) {
    this.radius = radius;
  }

  public SymbolList getSequence() {
    return symList;
  }

  public void setSequence(SymbolList symList) {
    this.symList = symList;
  }

  public double getOffset() {
    return offset;
  }

  public void setOffset(double offset) {
    this.offset = offset;
  }

  public CircularRenderer getRenderer() {
    return renderer;
  }

  public void setRenderer(CircularRenderer renderer) {
    this.renderer = renderer;
  }

  public synchronized void paintComponent(Graphics g) {
    super.paintComponent(g);
    if(!isActive()) return;

    double depth = renderer.getDepth(ctxt);

    Graphics2D g2 = (Graphics2D) g;
    g2.translate((depth + radius), (depth + radius));

    renderer.paint(g2, ctxt);
  }

  private boolean isActive() {
    return renderer != null;
  }

  private final class CTXT
  implements CircularRendererContext {
    public double getOffset() {
      return offset;
    }

    public double getAngle(int indx) {
      return ((double) indx) * 2.0 * Math.PI / ((double) symList.length()) + offset;
    }

    public int getIndex(double angle) {
      return (int) ( (angle - offset) * ((double) symList.length()) / (2.0 * Math.PI));
    }

    public double getRadius() {
      return radius;
    }

    public SymbolList getSymbols() {
      return symList;
    }

    public FeatureHolder getFeatures() {
      if(symList instanceof FeatureHolder) {
        return (FeatureHolder) symList;
      } else {
        return FeatureHolder.EMPTY_FEATURE_HOLDER;
      }
    }
  }
}
