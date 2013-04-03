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

import java.util.List;
import java.util.NoSuchElementException;
import java.util.Set;

import org.biojava.bio.Annotatable;
import org.biojava.bio.BioException;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.utils.ChangeType;

/**
 * <p>
 * The set of AtomicSymbols which can be concatenated together to make a
 * SymbolList.
 * </p>
 *
 * <p>
 * A non-atomic symbol is considered to be contained within this alphabet if
 * all of the atomic symbols that it could match are members of this alphabet.
 * </p>
 *
 * @author Matthew Pocock
 * @author Thomas Down
 */
 
public interface Alphabet extends Annotatable {
  /** 
   * <p>
   * This ChangeType indicates that some symbols have been added or removed from
   * the alphabet. The current and previous fields should indicate what symbols
   * were there originally, and what they have been replaced with.
   * <p>
   *
   * <p>
   * If the alphabet wishes to propagate that the symbol has changed state, then
   * previous and current should be null, but the chainedEvent property should
   * rever to the ChangeEvent on the unerlying Symbol.
   * </p>
   */
  public static ChangeType SYMBOLS = new ChangeType(
    "The set of symbols in this alphabet has changed.",
    "org.biojava.bio.symbol.Alphabet",
    "SYMBOLS"
  );
  
  /**
   * This signals that the available parsers have changed. If a parser is added,
   * it will appear in getChanged(). If it is removed, it will appear in
   * getPrevious().
   */
  public static ChangeType PARSERS = new ChangeType(
    "The set of parsers has changed.",
    "org.biojava.bio.symbol.Alphabet",
    "PARSERS"
  );
  
  /**
   * Get the name of the alphabet.
   *
   * @return  the name as a string.
   */
  String getName();

  /**
   * Return an ordered List of the alphabets which make up a
   * compound alphabet.  For simple alphabets, this will return
   * a singleton list of itself. The returned list should be immutable.
   *
   * @return a List of alphabets
   */
  List<Alphabet> getAlphabets();

  /**
   * <p>
   * Get a symbol from the Alphabet which corresponds
   * to the specified ordered list of symbols.
   * </p>
   *
   * <p>
   * The symbol at i in the list must be a member of the i'th alphabet in
   * getAlphabets. If all of the symbols in rl are atomic, then the resulting
   * symbol will also be atomic. If any one of them is an ambiguity symbol then
   * the resulting symbol will be the appropriate ambiguity symbol.
   * </p>
   *
   * @param rl A list of Symbol instances
   * @throws IllegalSymbolException if the members of rl are
   *            not Symbols over the alphabets returned from
   *            <code>getAlphabets</code>
   */
  Symbol getSymbol(List<Symbol> rl) 
    throws IllegalSymbolException;

  /**
   * <p>
   * Get a symbol that represents the set of symbols in syms.
   * </p>
   *
   * <p>
   * Syms must be a set of Symbol instances each of which is contained within
   * this alphabet. This method is used to retrieve ambiguity symbols.
   * </p>
   *
   * @param syms  the Set of Symbols that will be found in getMatches of the
   *            returned symbol
   * @return a Symbol (possibly fly-weighted) for the Set of symbols in syms
   */
  Symbol getAmbiguity(Set<Symbol> syms)
  throws IllegalSymbolException;
  
  /**
   * <p>
   * Get the 'gap' ambiguity symbol that is most appropriate for this alphabet.
   * </p>
   *
   * <p>
   * In general, this will be a BasisSymbol that represents a list of
   * AlphabetManager.getGapSymbol() the same length as the getAlphabets list.
   * </p>
   *
   * @return the appropriate gap Symbol instance
   */
  Symbol getGapSymbol();
  
  /**
   * <p>
   * Returns whether or not this Alphabet contains the symbol.
   * </p>
   *
   * <p>
   * An alphabet contains an ambiguity symbol iff the ambiguity symbol's
   * getMatches() returns an alphabet that is a proper sub-set of this
   * alphabet. That means that every one of the symbols that could match the
   * ambiguity symbol is also a member of this alphabet.
   * </p>
   *
   * @param s the Symbol to check
   * @return  boolean true if the Alphabet contains the symbol and false otherwise
   */
  boolean contains(Symbol s);

  /**
   * <p>
   * Throws a precanned IllegalSymbolException if the symbol is not contained
   * within this Alphabet.
   * </p>
   *
   * <p>
   * This function is used all over the code to validate symbols as they enter
   * a method. Also, the code is littered with catches for
   * IllegalSymbolException. There is a preferred style of handling this,
   * which should be covererd in the package documentation.
   * </p>
   *
   * @param s the Symbol to validate
   * @throws  IllegalSymbolException if r is not contained in this alphabet
   */
  void validate(Symbol s) throws IllegalSymbolException;
  
  /**
   * <p>
   * Get a SymbolTokenization by name.
   * </p>
   *
   * <p>
   * The parser returned is guaranteed to return Symbols and SymbolLists that
   * conform to this alphabet.
   * </p>
   *
   * <p>
   * Every alphabet should have a SymbolTokenzation under the name 'token' that
   * uses the symbol token characters to translate a string into a
   * SymbolList. Likewise, there should be a SymbolTokenization under the name
   * 'name' that uses symbol names to identify symbols. Any other names may
   * also be defined, but the behavior of the returned SymbolTokenization is
   * not defined here.
   * </p>
   * <p>
   * A SymbolTokenization under the name 'default' should be defined for all
   * sequences, that determines the behavior when printing out a
   * sequence. Standard behavior is to define the 'token' SymbolTokenization
   * as default if it exists, else to define the 'name' SymbolTokenization as
   * the default, but others are possible.
   * </p>
   *
   * @param name  the name of the parser
   * @return  a parser for that name
   * @throws NoSuchElementException if the name is unknown
   * @throws BioException if for any reason the tokenization could not be built
   * @since 1.2
   */
    
    public SymbolTokenization getTokenization(String name) throws BioException;
  
  /**
   * A really useful static alphabet that is always empty.
   */
  static final FiniteAlphabet EMPTY_ALPHABET = new EmptyAlphabet();
}
