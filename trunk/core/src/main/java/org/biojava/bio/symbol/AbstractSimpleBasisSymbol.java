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

import java.util.Collections;
import java.util.List;

import org.biojava.bio.Annotation;
import org.biojava.bio.BioError;
import org.biojava.utils.ListTools;

/**
 * A basic implementation of BasisSymbol.
 *
 * If you wish to construct new Symbols, you should normally do so via utility methods
 * on <code>AlphabetManager</code>.
 *
 * This may be a useful base class for custom implementations.
 *
 * @author Matthew Pocock
 */

abstract class AbstractSimpleBasisSymbol
        extends SimpleSymbol
        implements BasisSymbol
{
  protected List symbols;

  protected AbstractSimpleBasisSymbol(
    Annotation annotation, List symbols
  ) throws IllegalSymbolException {
    this(annotation);
    if(symbols == null) {
      throw new NullPointerException("symbols can't be null");
    }
    if(symbols.size() == 0) {
      throw new IllegalSymbolException(
        "Can't create BasisSymbol for an empty list. Use the Gap symbol."
      );
    }
    this.symbols = ListTools.createList(symbols);
  }

  protected AbstractSimpleBasisSymbol(
    Annotation annotation
  ) {
    super(annotation);
  }

  public AbstractSimpleBasisSymbol(
    Annotation annotation,
    Alphabet matches
  ) {
    this(annotation);
    this.matches = matches;
    this.symbols = Collections.nCopies(1, this);
  }

  public AbstractSimpleBasisSymbol(
    Annotation annotation,
    List symbols,
    Alphabet matches
  ) throws IllegalSymbolException {
    this(annotation, symbols);
    this.matches = matches;
  }

  public final List getSymbols() {
    if(symbols == null) {
      symbols = createSymbols();
    }
    if(symbols.size() == 0) {
      throw new BioError(
        "Assertion Failure: symbols array is of length 0 in " + this +
        "\n\tname: " + getName() +
        "\n\tsymbols: " + this.symbols +
        "\n\tmatches: " + this.matches
      );
    }
    return symbols;
  }

  protected List createSymbols() {
    throw new BioError("Assertion Failure: Symbols list is null");
  }
}
