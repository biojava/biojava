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


package org.biojava.bio.seq.io;

import org.biojava.bio.Annotatable;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;

/**
 * Encapsulate a mapping between BioJava Symbol objects and
 * some string representation.
 *
 * @author Thomas Down
 * @since 1.2
 */

public interface SymbolTokenization extends Annotatable {
  public final static class TokenType {
    private String type;

    private TokenType(String type) {
      this.type = type;
    }

    public String toString()
    {
      return "TokenType:" + type;
    }
  }

  public final static TokenType CHARACTER = new TokenType("CHARACTER");
  public final static TokenType FIXEDWIDTH = new TokenType("FIXEDWIDTH");
  public final static TokenType SEPARATED = new TokenType("SEPARATED");
  public final static TokenType UNKNOWN = new TokenType("UNKNOWN");

  /**
   * The alphabet to which this tokenization applies.
   */

  public Alphabet getAlphabet();

  /**
   * Determine the style of tokenization represented by this object.
   */

  public TokenType getTokenType();


  /**
   * Returns the symbol for a single token.
   * <p>
   * The Symbol will be a member of the alphabet. If the token is not recognized
   * as mapping to a symbol, an exception will be thrown.
   *
   * @param token the token to retrieve a Symbol for
   * @return the Symbol for that token
   * @throws IllegalSymbolException if there is no Symbol for the token
   */

  public Symbol parseToken(String token)
          throws IllegalSymbolException;

  /**
   * Return an object which can parse an arbitrary character stream into
   * symbols.
   *
   * @param listener The listener which gets notified of parsed symbols.
   */

  public StreamParser parseStream(SeqIOListener listener);

  /**
   * Return a token representing a single symbol.
   *
   * @param sym The symbol
   * @throws IllegalSymbolException if the symbol isn't recognized.
   */

  public String tokenizeSymbol(Symbol sym) throws IllegalSymbolException;

  /**
   * Return a string representation of a list of symbols.
   *
   * @param symList A SymbolList
   * @throws IllegalAlphabetException if alphabets don't match
   */

  public String tokenizeSymbolList(SymbolList symList) throws IllegalAlphabetException, IllegalSymbolException;
}
