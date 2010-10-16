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
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.Shape;
import java.awt.font.FontRenderContext;
import java.awt.font.GlyphVector;
import java.awt.geom.AffineTransform;
import java.awt.geom.Rectangle2D;

import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.bio.BioRuntimeException;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.bio.symbol.AtomicSymbol;
import org.biojava.bio.symbol.IllegalSymbolException;

/**
 * A BlockPainter that renders letters in proportion to the size of the signal.
 *
 * @author Matthew Pocock
 * @author Thomas Down
 */
public class TextBlock implements BlockPainter {
  private Font logoFont = new Font("Tahoma", Font.PLAIN, 12);

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
    this.logoFont = logoFont;
  }

  public void paintBlock(LogoContext ctxt, Rectangle2D block, AtomicSymbol sym) {
    Graphics2D g2 = ctxt.getGraphics();
    SymbolStyle style = ctxt.getStyle();
    Distribution dist = ctxt.getDistribution();
    SymbolTokenization toke = null;
    try {
        toke = dist.getAlphabet().getTokenization("token");
    } catch (BioException ex) {
        throw new BioRuntimeException(ex);
    }

    FontRenderContext frc = g2.getFontRenderContext();
    try {
        GlyphVector gv = logoFont.createGlyphVector(frc, toke.tokenizeSymbol(sym));

        Shape outline = gv.getOutline();
        Rectangle2D oBounds = outline.getBounds2D();

        AffineTransform at = new AffineTransform();
        at.setToTranslation(block.getX(), block.getY());
        at.scale(
                 block.getWidth() / oBounds.getWidth(),
                 block.getHeight() / oBounds.getHeight()
                 );
        at.translate(-oBounds.getMinX(), -oBounds.getMinY());
        outline = at.createTransformedShape(outline);

        try {
            g2.setPaint(style.fillPaint(sym));
        } catch (IllegalSymbolException ire) {
            g2.setPaint(Color.black);
        }
        g2.fill(outline);

        try {
            g2.setPaint(style.outlinePaint(sym));
        } catch (IllegalSymbolException ire) {
            g2.setPaint(Color.gray);
        }
        g2.draw(outline);
    } catch (IllegalSymbolException ex) {
        throw new BioError("Couldn't tokenize symbol", ex);
    }
  }
}
