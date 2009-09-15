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
import java.awt.Graphics2D;
import java.awt.geom.Rectangle2D;

import org.biojava.bio.symbol.AtomicSymbol;
import org.biojava.bio.symbol.IllegalSymbolException;

/**
 * A painter that just draws a block (or bar).
 *
 * @author Matthew Pocock
 */
public class PlainBlock implements BlockPainter {
  public void paintBlock(LogoContext ctxt, Rectangle2D block, AtomicSymbol sym) {
    Graphics2D g2 = ctxt.getGraphics();
    SymbolStyle style = ctxt.getStyle();
    
    try {
      g2.setPaint(style.fillPaint(sym));
    } catch (IllegalSymbolException ire) {
      g2.setPaint(Color.black);
    }
    g2.fill(block);
    
    try {
      g2.setPaint(style.outlinePaint(sym));
    } catch (IllegalSymbolException ire) {
      g2.setPaint(Color.gray);
    }
    g2.draw(block);
  }
}
