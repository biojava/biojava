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

/**
 * <p>
 * An encapsulation of the way symbols map to bit-patterns.
 * </p>
 *
 * <p>
 * A packing will encapsulate the process of converting between symbols and
 * bit-patterns. For example, in DNA you could use 00, 01, 10 and 11 to
 * represent the four bases (a, g, c, t). Many applications may require a
 * specific packing. You may need to store full ambiguity information, or
 * perhaps you can discard this capability to reduce stoorage space. You may
 * care about the bit-pattern produced because you need interoperability or
 * an algorithm needs to be fed correctly, or you may not care about the
 * packing at all. This interface is here to allow you to chose the most
 * appropreate packing for your task.
 * </p>
 *
 * @author Matthew Pocock
 */
public interface Packing {
  /**
   * The FiniteAlphabet this packing is for.
   *
   * @return  the FiniteAlphabet that we can pack
   */
  FiniteAlphabet getAlphabet();
  
  /**
   * <p>
   * Return a byte representing the packing of a symbol.
   * The bits will be from 1 >> 0 through to 1 >> (wordSize - 1).
   * </p>
   *
   * @param sym  the Symbol to pack
   * @return  a byte containing the packed symbol
   * @throws IllegalSymbolException if sym is not in getAlphabet().
   */
  byte pack(Symbol sym)
  throws IllegalSymbolException;
  
  /**
   * <p>
   * Return the symbol for a packing.
   * </p>
   *
   * @param packed  the byte pattern for a Symbol
   * @return the Symbol that was packed
   * @throws IllegalSymbolException if the packing doesn't represent a valid
   *         Symbol
   */
  Symbol unpack(byte packed)
  throws IllegalSymbolException;
  
  /**
   * <p>
   * The number of bits required to pack a symbol.
   * </p>
   *
   * @return the word size as a byte
   */
  byte wordSize();
  
  /**
   * <p>
   * Flag to state if ambiguities are stored.
   * </p>
   *
   * <p>
   * Packings are free to either store ambiguity information or to discard
   * it (presumably converting all ambiguities to a standard AtomicSymbol
   * and then packing that). You can check wether ambiguities are handled
   * by calling this method.
   * </p>
   *
   * @return true if ambiguities are stored, false otherwise
   */
  boolean handlesAmbiguity();
}
