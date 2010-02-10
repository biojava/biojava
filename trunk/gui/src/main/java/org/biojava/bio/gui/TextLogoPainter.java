/*
 *                    BioJava development code
 *
 * This code may be freely disttributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be disttributed with the code.  If you do not have a copy,
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


package org.biojava.bio.gui;

import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.Rectangle;
import java.awt.Shape;
import java.awt.font.FontRenderContext;
import java.awt.font.GlyphVector;
import java.awt.geom.AffineTransform;
import java.awt.geom.Rectangle2D;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.beans.PropertyChangeSupport;
import java.util.Comparator;
import java.util.Iterator;
import java.util.SortedSet;
import java.util.TreeSet;

import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.bio.BioRuntimeException;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;

/**
 * A logo painter that paints in stacked letters.
 * The total height of the letters is
 * proportional to the total informaton in the state. The height of each letter
 * is proportional to its emission probability. The most likely letter is drawn
 * highest.
 *
 * @author Matthew Pocock
 */
public class TextLogoPainter implements LogoPainter {
  /**
   * A comparator to set up our letters & information scores nicely.
   */
  private static final Comparator COMP = new ResValComparator();

  /**
   * Supports the bean property logoFont.
   */
  private PropertyChangeSupport pcs;

  /**
   * The property for the logoFont.
   */
  private Font logoFont;

  /**
   * Retrieve the current font.
   *
   * @return the current logo font
   */
  public Font getLogoFont() {
    return logoFont;
  }

  /**
   * Set the current logo font.
   *
   * @param logoFont the new Font to render the logo letters in
   */
  public void setLogoFont(Font logoFont) {
    firePropertyChange("logoFont", this.logoFont, logoFont);
    this.logoFont = logoFont;
  }

  public void addPropertyChangeListener(PropertyChangeListener listener) {
    pcs.addPropertyChangeListener(listener);
  }

  public void removePropertyChangeListener(PropertyChangeListener listener) {
    pcs.removePropertyChangeListener(listener);
  }

  public void addPropertyChangeListener(String propertyName, PropertyChangeListener listener) {
    pcs.addPropertyChangeListener(propertyName, listener);
  }

  public void removePropertyChangeListener(String propertyName, PropertyChangeListener listener) {
    pcs.removePropertyChangeListener(propertyName, listener);
  }

  public void firePropertyChange(String propertyName, Object oldValue, Object newValue) {
    pcs.firePropertyChange(propertyName, oldValue, newValue);
  }

  public void firePropertyChange(String propertyName, int oldValue, int newValue) {
    pcs.firePropertyChange(propertyName, oldValue, newValue);
  }

  public void firePropertyChange(String propertyName, boolean oldValue, boolean newValue) {
    pcs.firePropertyChange(propertyName, oldValue, newValue);
  }

  public void firePropertyChange(PropertyChangeEvent evt) {
    pcs.firePropertyChange(evt);
  }

  public boolean hasListeners(String propertyName) {
    return pcs.hasListeners(propertyName);
  }

  public void paintLogo(LogoContext ctxt) {
    Graphics2D g2 = ctxt.getGraphics();
    Distribution dist = ctxt.getDistribution();
    SymbolStyle style = ctxt.getStyle();
    SymbolTokenization toke = null;
    try {
        toke = dist.getAlphabet().getTokenization("token");
    } catch (BioException ex) {
        throw new BioRuntimeException(ex);
    }

    Rectangle bounds = ctxt.getBounds();
    double width = bounds.getWidth();
    double height = bounds.getHeight();
    double base = bounds.getY() + bounds.getHeight();
    
    /* This used to have some built-in scaling support, but I've disabled this because
     * DistributionLogo does scaling too!
     */
    
    // double scale = height * (
    //  DistributionLogo.totalInformation(dist) /
    //  DistributionLogo.totalBits(dist)
    // );

    SortedSet info = new TreeSet(COMP);

    try {
      for(
        Iterator i = ((FiniteAlphabet) dist.getAlphabet()).iterator();
        i.hasNext();
      ) {
        Symbol s = (Symbol) i.next();
        info.add(new ResVal(s, dist.getWeight(s) * height));
      }
    } catch (IllegalSymbolException ire) {
      throw new BioError("Symbol distsapeared from dist alphabet", ire);
    }

    FontRenderContext frc = g2.getFontRenderContext();
    for(Iterator i = info.iterator(); i.hasNext();) {
      ResVal rv = (ResVal) i.next();

      String s = null;
      try {
          s = toke.tokenizeSymbol(rv.getToken());
      } catch (IllegalSymbolException ex) {
          throw new BioRuntimeException(ex);
      }
      GlyphVector gv = logoFont.createGlyphVector(frc, s);
      Shape outline = gv.getOutline();
      Rectangle2D oBounds = outline.getBounds2D();

      AffineTransform at = new AffineTransform();
      at.setToTranslation(0.0, base-rv.getValue());
      at.scale(width / oBounds.getWidth(), rv.getValue() / oBounds.getHeight());
      at.translate(-oBounds.getMinX(), -oBounds.getMinY());
      outline = at.createTransformedShape(outline);

      try {
        g2.setPaint(style.fillPaint(rv.getToken()));
      } catch (IllegalSymbolException ire) {
        g2.setPaint(Color.black);
      }
      g2.fill(outline);

      try {
        g2.setPaint(style.outlinePaint(rv.getToken()));
      } catch (IllegalSymbolException ire) {
        g2.setPaint(Color.gray);
      }
      g2.draw(outline);

      base -= rv.getValue();
    }
  }

  public TextLogoPainter() {
    pcs = new PropertyChangeSupport(this);
    logoFont = new Font("Tahoma", Font.PLAIN, 12);
  }

  /**
   * A symbol/information tuple.
   */
  private static class ResVal {
    private Symbol symbol;
    private double value;

    public final Symbol getToken() {
      return symbol;
    }

    public final double getValue() {
      return value;
    }

    public ResVal(Symbol sym, double val) {
      symbol = sym;
      value = val;
    }
  }

  /**
   * The comparator for comparing symbol/information tuples.
   */
  private static class ResValComparator implements Comparator {
    public final int compare(Object o1, Object o2) {
      ResVal rv1 = (ResVal) o1;
      ResVal rv2 = (ResVal) o2;

      double diff = rv1.getValue() - rv2.getValue();
      if(diff < 0) return -1;
      if(diff > 0) return +1;
      return rv1.getToken().getName().compareTo(rv2.getToken().getName());
    }
  }
}
