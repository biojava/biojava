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

/**
 * An reverse view of another SymbolList.
 *
 * @author Matthew Pocock
 */
class ReverseSymbolList
extends AbstractSymbolList implements Serializable {
  /**
   * The source sequence that we will transliterate.
   */
  private final SymbolList source;
  /**
   * Retrieve the underlying SymbolList being viewed.
   *
   * @return the source SymbolList
   */
  public SymbolList getSource() {
    return source;
  }
  
  /**
   * Create a reverse view of source.
   */
  public ReverseSymbolList(SymbolList source) {
    this.source = source;
  }

  public Alphabet getAlphabet() {
    return source.getAlphabet();
  }

  public int length() {
    return source.length();
  }
  
  public Symbol symbolAt(int index)
  throws IndexOutOfBoundsException {
    return source.symbolAt(length() - index + 1);
  }
}
