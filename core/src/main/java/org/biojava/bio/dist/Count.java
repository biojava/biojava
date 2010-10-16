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


package org.biojava.bio.dist;

import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.AtomicSymbol;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.Changeable;

/**
 * <p>
 * An encapsulation of a count over the Symbols within an alphabet.
 * </p>
 *
 * <p>
 * A Count is effectively a vector of counts from an Alphabet.
 * </p>
 *
 * @author Matthew Pocock
 * @since 1.1
 */
public interface Count extends Changeable {
  /**
   * <p>
   * Whenever a component count changes the values that would be returned by
   * getCount, they should fire a ChangeEvent with this object as the type.
   * </p>
   *
   * <p>
   * If the whole count changes, then the change and previous fields of
   * the ChangeEvent should be left null. If only a single weight is modified,
   * then change should be of the form Object[] { symbol, new Double(newVal) }
   * and previous should be of the form Object[] { symbol, new Double(oldVal) }
   * </p>
   */
  public static final ChangeType COUNTS = new ChangeType(
    "distribution weights changed",
    "org.biojava.bio.dist.Count",
    "COUNTS"
  );
  
  /**
   * The alphabet from which this Count is over.
   *
   * @return  the Alphabet associated with this Count
   */
  Alphabet getAlphabet();
    
  /**
   * Return the counts for a given Symbol.
   *
   * @param s the Symbol
   * @return  the number of counts for this symbol
   * @throws IllegalSymbolException if s is not from this Count's alphabet
   */
  double getCount(AtomicSymbol s) throws IllegalSymbolException;
  
  /**
   * Set the count for the Symbol s.
   *
   * @param s the Symbol emitted
   * @param c  the new count for the Symbol
   * @throws IllegalSymbolException if s is not from this state's alphabet, or
   *         if it is an ambiguity symbol and the implementation can't handle
   *         this case
   * @throws ChangeVetoException if this distribution does not allow counts
   *         to be tampered with, or if one of the listeners vetoed this change
   */
  void setCount(AtomicSymbol s, double c)
  throws IllegalSymbolException, ChangeVetoException;
  
  
  /**
   * Set the probability or odds that Symbol s is emitted by this state.
   *
   * @param s the Symbol emitted
   * @param c  the delta to add to the count for the Symbol
   * @throws IllegalSymbolException if s is not from this state's alphabet, or
   *         if it is an ambiguity symbol and the implementation can't handle
   *         this case
   * @throws ChangeVetoException if this Count does not allow counts
   *         to be tampered with, or if one of the listeners vetoed this change
   */
  void increaseCount(AtomicSymbol s, double c)
  throws IllegalSymbolException, ChangeVetoException;
  
  /**
   * Set the counts in this Counts to be equal to the counts in c.
   *
   * @param c  the Count object to copy the counts from
   * @throws IllegalAlphabetException if c has a different Alphabet to this
   *         Count
   * @throws ChangeVetoException if this Count does not allow the counts to be
   *         tampered with, or if one of the listeners vetoed this change
   */
  void setCounts(Count c)
  throws IllegalAlphabetException, ChangeVetoException;
  
  /**
   * Reset all the counts to zero.
   *
   * @throws ChangeVetoException if this Count does not allow the counts to be
   *         tampered with, or if one of the listeners vetoed this change
   */
  void zeroCounts()
  throws ChangeVetoException;
}
