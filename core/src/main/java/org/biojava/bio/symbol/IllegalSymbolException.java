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

import org.biojava.bio.BioException;

/**
 * The exception to indicate that a symbol is not valid within a context.
 * <p>
 * The usual reason for throwing an IllegalSymbolException is that you are
 * trying to add a symbol to a sequence with an alpabet that does not contain
 * the symbol. This is the sequence/alphabet equivalent of a ClassCastException
 * for objects.
 * <p>
 * Frequently, these excepions are actualy generated from Alphabet.validate.
 *
 * @author Matthew Pocock
 */
public class IllegalSymbolException extends BioException {
  private final Symbol sym;

  /**
   * Retrieve the symbol that caused this exception, or null.
   */
  public Symbol getSymbol() {
    return sym;
  }

  /**
   * Make the exception with a message.
   */
  public IllegalSymbolException(String message) {
    this(null, null, message);
  }

  /**
   * Make the exception with a message and a symbol.
   */
  public IllegalSymbolException(Symbol sym, String message) {
    this(null, sym, message);
  }

  public IllegalSymbolException(Throwable cause, String message) {
    this(cause, null, message);
  }

  public IllegalSymbolException(Throwable cause, Symbol sym, String message) {
    super(message, cause);
    this.sym = sym;
  }
}
