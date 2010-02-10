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

import java.awt.Paint;

import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;

/**
 * The interface for things that say how to paint a symbol.
 * <p>
 * Given a symbol, this allows you to get the color to outline or fill the
 * glyphs for rendering the symbol. This may be something as simple as colouring
 * dots on a scatter-plot, or labeling a key, or it may be as complicated as
 * sequence logos.
 *
 * @author Matthew Pocock
 */
public interface SymbolStyle {
  /**
   * Return the outline paint for a symbol.
   *
   * @param s the symbol to outline
   * @return the Paint to use
   * @throws IllegalSymbolException if this SymbolStyle can not handle the
   *         symbol
   */
  Paint outlinePaint(Symbol s) throws IllegalSymbolException;

  /**
   * Return the fill paint for a symbol.
   *
   * @param s the symbol to fill
   * @return the Paint to use
   * @throws IllegalSymbolException if this SymbolStyle can not handle the
   *         symbol
   */
  Paint fillPaint(Symbol s) throws IllegalSymbolException;
}
