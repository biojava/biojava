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

package org.biojava.bio.symbol;

import java.io.Serializable;
import java.util.Collections;

import org.biojava.bio.BioError;

/**
 * A view of windows onto another SymbolList.
 * <p>
 * In practice, what this means is that you can view a DNA sequence as codons which
 * do not overlap.
 *
 * @author Matthew Pocock
 */
class WindowedSymbolList
extends AbstractSymbolList implements Serializable {
  /**
   * The source sequence that we will transliterate.
   */
  private final SymbolList source;

  /**
   * The alphabet of window Symbols.
   */
  private final Alphabet alpha;

  /**
   * The width of the window.
   */
  private final int width;

  /**
   * Retrieve the underlying SymbolList being viewed.
   *
   * @return the source SymbolList
   */
  public SymbolList getSource() {
    return source;
  }

  /**
   * Create a WindowedSymbolList with the given window width.
   */
  public WindowedSymbolList(SymbolList source, int width)
  throws IllegalArgumentException {
    if( (source.length() % width) != 0 ) {
      throw new IllegalArgumentException(
        "The source length must be divisible by the window width: " +
        source.length() + " % " + width + " = " + (source.length() % width)
      );
    }
    this.source = source;
    Alphabet a = source.getAlphabet();
    this.alpha = AlphabetManager.getCrossProductAlphabet(
      Collections.nCopies(width, a)
    );
    this.width = width;
  }

  public Alphabet getAlphabet() {
    return alpha;
  }

  public int length() {
    return source.length() / width;
  }

  public Symbol symbolAt(int index)
  throws IndexOutOfBoundsException {
    if(index < 1 || index > length()) {
      throw new IndexOutOfBoundsException(
        "index must be within (1 .. " +
        length() + "), not " + index
      );
    }

    index = (index-1)*width + 1;

    try {
      return alpha.getSymbol(source.subList(index, index+width-1).toList());
    } catch (IllegalSymbolException iae) {
      throw new BioError("Alphabet changed underneath me", iae);
    }
  }
}
