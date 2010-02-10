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

import java.util.Iterator;

import org.biojava.utils.ChangeVetoException;

/**
 * An alphabet over a finite set of Symbols.
 * <p>
 * This interface makes the distinction between an alphabet over a finite (and
 * possibly small) number of symbols and an Alphabet over an infinite
 * (or extremely large) set of symbols. Within a FiniteAlphabet, the == operator
 * should be sufficient to decide upon equality for all AtomicSymbol instances.
 * <p>
 * The alphabet functions as the repository of objects in the fly-weight design
 * pattern. Only symbols within an alphabet should appear in object that claim
 * to use the alphabet - otherwise something is in error.
 *
 * @author Matthew Pocock
 */
public interface FiniteAlphabet extends Alphabet {
  /**
   * The number of symbols in the alphabet.
   *
   * @return the size of the alphabet
   */
  int size();

  /**
   * Retrieve an Iterator over the AtomicSymbols in this FiniteAlphabet.
   * <p>
   * Each AtomicSymbol as for which this.contains(as) is true will be returned
   * exactly once by this iterator in no specified order.
   *
   * @return an Iterator over the contained AtomicSymbol objects
   */
  Iterator<Symbol> iterator();


  /**
   * Adds a symbol to this alphabet.
   * <p>
   * If the symbol matches more than one AtomicSymbol, then each symbol matching
   * it will be added.
   *
   * @param s the Symbol to add
   * @throws IllegalSymbolException if the symbol is null, or if for any reason
   *         it can't be added
   * @throws ChangeVetoException  if either the alphabet doesn't allow symbols
   *         to be added, or the change was vetoed
   */
  public void addSymbol(Symbol s)
  throws IllegalSymbolException, ChangeVetoException;


  /**
   * Remove a symbol from this alphabet.
   * <p>
   * If the symbol matches multiple AtomicSymbols, then each matching symbol it
   * will be removed.
   *
   * @param s the Symbol to removeintGot
   * @throws IllegalSymbolException if the symbol is null, or if for any reason
   *         it can't be removed
   * @throws ChangeVetoException  if either the alphabet doesn't allow symbols
   *         to be added, or the change was vetoed
   */
  public void removeSymbol(Symbol s)
  throws IllegalSymbolException, ChangeVetoException;
}
