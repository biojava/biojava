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

import java.awt.Rectangle;
import java.awt.geom.Rectangle2D;
import java.util.Iterator;

import org.biojava.bio.BioError;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.symbol.AtomicSymbol;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalSymbolException;

/**
 * A logo painter that paints in bars. The total height of the bars is
 * proportional to the total informaton in the state.
 *
 * @author Matthew Pocock
 */
public class BarLogoPainter implements LogoPainter {
  public void paintLogo(LogoContext lCtxt) {
    Distribution dis = lCtxt.getDistribution();
    BlockPainter blockPainter = lCtxt.getBlockPainter();

    Rectangle bounds = lCtxt.getBounds();
    double width = bounds.getWidth();
    double stepWidth = width / (double) ((FiniteAlphabet) dis.getAlphabet()).size();
    double height = bounds.getHeight();

    double w = 0.0;
    for(
      Iterator i = ((FiniteAlphabet) dis.getAlphabet()).iterator();
      i.hasNext();
    ) {
      AtomicSymbol s = (AtomicSymbol) i.next();
      double rh = 0.0;

      try {
        rh = dis.getWeight(s) * height;
      } catch (IllegalSymbolException ire) {
        throw new BioError("State alphabet has changed while painting", ire);
      }

      Rectangle2D outline = new Rectangle2D.Double(
        bounds.getX() + w,
        bounds.getY() + height - rh,
        stepWidth,
        rh
      );

      blockPainter.paintBlock(lCtxt, outline, s);

      w += stepWidth;
    }
  }

  public BarLogoPainter() {}
}
