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


package org.biojava.bio.gui;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Rectangle;
import java.awt.RenderingHints;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.util.Iterator;

import javax.swing.JComponent;

import org.biojava.bio.BioError;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;

/**
 * The GUI component for rendering a DistributionLogo. By default, this uses the text
 * logo style - with letters stacked on top of one another, scaled by the total
 * information in the dist, and uses a PlainStyle colorer that outlines in
 * black, and fills in grey.
 *
 * @author Matthew Pocock
 */
public class DistributionLogo extends JComponent {
  /**
   * The default logo painter.
   */
  private static final LogoPainter DEFAULT_LOGO_PAINTER = new TextLogoPainter();
  private static final BlockPainter DEFAULT_BLOCK_PAINTER = new PlainBlock();

  /**
   * A usefull constant to keep arround.
   */
  private static double bits = Math.log(2.0);

  /**
   * The dist to render.
   */
  private Distribution dist;

  /**
   * The logoPainter property.
   */
  private LogoPainter logoPainter = DEFAULT_LOGO_PAINTER;

  /**
   * The style property.
   */
  private SymbolStyle style = new PlainStyle(Color.black, Color.gray);

  /**
   * The information/plain property
   */
  boolean scaleByInformation = true;

  private BlockPainter blockPainter = DEFAULT_BLOCK_PAINTER;
  private RenderingHints renderingHints = null;

  /**
   * Retrieve the currently rendered dist.
   *
   * @return  a Distribution
   */
  public Distribution getDistribution() {
    return dist;
  }

  /**
   * <p>
   * Set the dist to render.
   * </p>
   *
   * <p>
   * The dist must be over a FiniteAlphabet so that we can draw the numbers
   * for each Symbol.
   * </p>
   *
   * @param dist the new Distribution to render
   */
  public void setDistribution(Distribution dist)
  throws IllegalAlphabetException {
    firePropertyChange("dist", this.dist, dist);
    this.dist = dist;
  }

  /**
   * Retrieve the current logo painter.
   *
   * @return  the LogoPainter used to render the dist
   */
  public LogoPainter getLogoPainter() {
    return logoPainter;
  }

  /**
   * <p>
   * Set the logo painter.
   * </p>
   *
   * <p>
   * This will alter the way that the dist is rendered to screen.
   * </p>
   *
   * @param logoPainter the new logoPainter
   */
  public void setLogoPainter(LogoPainter logoPainter) {
    firePropertyChange("logoPainter", this.logoPainter, logoPainter);
    this.logoPainter = logoPainter;
  }

  /**
   * Retrieve the current style.
   *
   * @return the current SymbolStyle
   */
  public SymbolStyle getStyle() {
    return style;
  }

  /**
   * <p>
   * Set the symbol style.
   * </p>
   *
   * <p>
   * This will change the outline and fill paints for the logos
   * </p>
   *
   * @param style the new SymbolStyle to use
   */
  public void setStyle(SymbolStyle style) {
    firePropertyChange("style", this.style, style);
    this.style = style;
  }

  public boolean isScaleByInformation() {
    return scaleByInformation;
  }

  public void setScaleByInformation(boolean scale) {
    this.scaleByInformation = scale;
  }

  public BlockPainter getBlockPainter() {
    return blockPainter;
  }

  public void setBlockPainter(BlockPainter blockPainter) {
    this.blockPainter = blockPainter;
  }

  /**
   * Create a new DistributionLogo object. It will set up all the properties except the
   * dist to render.
   */
  public DistributionLogo() {
    this.addPropertyChangeListener(new PropertyChangeListener() {
      public void propertyChange(PropertyChangeEvent pce) {
        String name = pce.getPropertyName();
        if(name.equals("dist") ||
           name.equals("logoPainter") ||
           name.equals("style") )
        {
          repaint();
        }
      }
    });

    Dimension d = new Dimension(20, 20);
    setMinimumSize(d);
    setPreferredSize(d);
  }

  /**
   * Calculate the information content of a symbol in bits.
   *
   * @param s the symbol to calculate for
   * @param dist the <code>Distribution</code> that the symbol comes from and
   * in whose context the information content will be calculated.
   * @throws IllegalSymbolException if r is not within the dist.
   */
  public static double entropy(Distribution dist, Symbol s) throws IllegalSymbolException {
    double p = dist.getWeight(s);
    if (p == 0.0) {
        return 0;
    }
    double lp = Math.log(p);

    return -p * lp / bits;
  }

  /**
   * Retrieve the maximal number of bits possible for this type of dist.
   *
   * @return maximum bits as a double
   */
  public static double totalBits(Distribution dist) {
    return Math.log(((FiniteAlphabet) dist.getAlphabet()).size()) / bits;
  }

  /**
   * <p>
   * Calculates the total information of the dist in bits.
   * </p>
   *
   * <p>
   * This calculates <code>totalBits - sum_r(entropy(r))</code>
   * </p>
   *
   * @return  the total information in the dist
   */
  public static double totalInformation(Distribution dist) {
    double inf = totalBits(dist);

    for(
      Iterator i = ((FiniteAlphabet) dist.getAlphabet()).iterator();
      i.hasNext();
    ) {
      Symbol s = (Symbol) i.next();
      try {
        inf -= entropy(dist, s);
      } catch (IllegalSymbolException ire) {
        throw new BioError(
        "Symbol evaporated while calculating information", ire);
      }
    }

    return inf;
  }

  /**
   * Transforms the graphics context so that it is in bits space,
   * and then requests the logo painter to fill the area.
   */
  public void paintComponent(Graphics g) {
    final Graphics2D g2 = (Graphics2D) g;

    if(renderingHints != null){
      g2.setRenderingHints(renderingHints);
    }

    Rectangle clip = g2.getClipBounds();
    if(isOpaque()) {
      g2.clearRect(clip.x, clip.y, clip.width, clip.height);
    }
    if(getDistribution() == null) {
      return;
    }

    final Rectangle bounds = getBounds();
    if(isScaleByInformation()) {
      int height = bounds.height;
      double scale = height * (totalInformation(getDistribution()) / totalBits(getDistribution()));
      bounds.height = (int) scale;
      bounds.y += (int) (height - scale);
    } 

    LogoContext ctxt = new LogoContext() {
      public Graphics2D getGraphics() {
        return g2;
      }
      public Distribution getDistribution() {
        return DistributionLogo.this.getDistribution();
      }
      public Rectangle getBounds() {
        return bounds;
      }
      public SymbolStyle getStyle() {
        return DistributionLogo.this.getStyle();
      }
      public BlockPainter getBlockPainter() {
        return DistributionLogo.this.getBlockPainter();
      }
    };
    getLogoPainter().paintLogo(ctxt);
  }

  public RenderingHints getRenderingHints() {
    return renderingHints;
  }

  public void setRenderingHints(RenderingHints renderingHints) {
    this.renderingHints = renderingHints;
  }
}
