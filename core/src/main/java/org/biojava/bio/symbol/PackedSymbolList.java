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

import org.biojava.utils.AssertionFailure;

/**
 * <p>
 * A SymbolList that stores symbols as bit-patterns in an array of longs.
 * </p>
 *
 * <p>
 * Bit-packed symbol lists are space efficient compared to the usual pointer
 * storage model employed by implementations like SimpleSymbolList. This
 * comes at the cost of encoding/decoding symbols from the storage. In
 * practice, the decrease in memory when storing large sequences makes
 * applications go quicker because of issues like page swapping.
 * </p>
 *
 * <p>
 * Symbols can be mapped to and from bit-patterns. The Pattern interface
 * encapsulates this. A SymbolList can then be stored by writing these
 * bit-patterns into memory. This implementation stores the bits
 * in the long elements of an array. The first symbol will be packed into
 * bits 0 through packing.wordLength()-1 of the long at index 0.
 * <p>
 *
 * <h2>Example Usage</h2>
 * <pre>
 * SymbolList symL = ...;
 * SymbolList packed = new PackedSymbolList(
 *   PackingFactory.getPacking(symL.getAlphabet(), true),
 *   symL
 * );
 * </pre>
 *
 * @author Matthew Pocock
 * @author David Huen (new constructor for Symbol arrays and some speedups)
 */
public class PackedSymbolList
  extends
    AbstractSymbolList
  implements
    java.io.Serializable
{
  private static final byte BITS_PER_ELEMENT = 64;

  private final Packing packing;
  private final long[] syms;
  private final int length;
  private final byte symsPerElement;
  private final byte mask;

  // scratch area for optimisations
  // WARNING: these variables constitute an opportunity
  // for things to go wrong when doing multithreaded access
  // via symbolAt().  Keep SymbolAt() synchronized so they
  // don't get changed during a lookup! Naaasssty.
  private int currentMin = Integer.MAX_VALUE;
  private int currentMax = Integer.MIN_VALUE;
  private long currentWord;
  private int wordsize;

  public Alphabet getAlphabet() {
    return packing.getAlphabet();
  }

  public int length() {
    return length;
  }

  /**
   * <p>
   * Create a new PackedSymbolList directly from a bit pattern.
   * </p>
   *
   * <p>
   * <em>Warning:</em> This is a risky developer method.
   * You must be sure that the syms array is packed in a
   * way that is consistent with the packing. Also, it is your
   * responsibility to ensure that the length is sensible.</em>
   * </p>
   *
   * @param packing the Packing used
   * @param syms a long array containing already packed symbols
   * @param length the length of the sequence packed in symbols
   */
  public PackedSymbolList(Packing packing, long[] syms, int length) {
    this.symsPerElement = (byte) (BITS_PER_ELEMENT / packing.wordSize());
    this.packing = packing;
    this.syms = syms;
    this.length = length;
    this.mask = calcMask(packing);
    wordsize = packing.wordSize();
  }

  /**
   * <p>
   * Create a new PackedSymbolList as a packed copy of another symbol list.
   * </p>
   *
   * <p>
   * This will create a new and independent symbol list that is a copy of
   * the symbols in symList. Both lists can be modified independently.
   * </p>
   *
   * @param packing the way to bit-pack symbols
   * @param symList the SymbolList to copy
   */
  public PackedSymbolList(Packing packing, SymbolList symList)
  throws IllegalAlphabetException {
    if(packing.getAlphabet() != symList.getAlphabet()) {
      throw new IllegalAlphabetException(
        "Can't pack with alphabet " + packing.getAlphabet() +
        " and symbol list " + symList.getAlphabet()
      );
    }

    try {
      this.symsPerElement = (byte) (BITS_PER_ELEMENT / packing.wordSize());
      this.packing = packing;
      this.length = symList.length();
      this.syms = new long[
        length / symsPerElement +
        ((length % symsPerElement == 0) ? 0 : 1)
      ];
      this.mask = calcMask(packing);
      wordsize = packing.wordSize();

      // pack the body of the sequence
      int ii = 0;
      for(int i = 0; i < (syms.length - 1); i++) {
//        int ii = i * symsPerElement;
        long l = 0;
        int jj = 0;
        for(int j = 0; j < symsPerElement; j++) {
//          int jj = j * packing.wordSize();
          long p = packing.pack(symList.symbolAt(ii + j + 1));
          l |= (long) ((long) p << (long) jj);
          jj += wordsize;
        }
        syms[i] = l;
        ii += symsPerElement;
      }

      // pack the final word
      if(syms.length > 0) {
        long l = 0;
        ii = (syms.length - 1) * symsPerElement;
        int jMax = symList.length() % symsPerElement;
        if(jMax == 0) {
          jMax = symsPerElement;
        }
        for(int j = 0; j < jMax; j++) {
          int jj = j * packing.wordSize();
          long p = packing.pack(symList.symbolAt(ii + j + 1));
          l |= (long) ((long) p << (long) jj);
        }
        syms[syms.length - 1] = l;
      }
    } catch (IllegalSymbolException ise) {
      throw new AssertionFailure("Assertion Failure: Symbol got lost somewhere", ise);
    }
  }

  /**
   * <p>
   * Create a new PackedSymbolList from an array of Symbols.
   * </p>
   *
   * <p>
   * This will create a new and independent SymbolList formed from the
   * the symbol array.
   * </p>
   *
   * @param packing the way to bit-pack symbols
   * @param symbols an array of Symbols
   * @param length the number of Symbols to process from symbols
   * @param alfa the alphabet from which the Symbols are drawn
   */
  public PackedSymbolList(Packing packing, Symbol [] symbols, int length, Alphabet alfa)
  throws IllegalAlphabetException,IllegalArgumentException {

    // verify that the alphabet is one I can deal with.
    if(packing.getAlphabet() != alfa) {
      throw new IllegalAlphabetException(
        "Can't pack with alphabet " + packing.getAlphabet() +
        " and symbol list " + alfa
      );
    }

    // check that array length makes sense
    if (symbols.length < length) {
      throw new IllegalArgumentException(
        "Symbol array size is too small to get " + length +
        "symbols from."
      );
    }

    try {
      this.symsPerElement = (byte) (BITS_PER_ELEMENT / packing.wordSize());
      this.packing = packing;
      this.length = length;
      this.syms = new long[
        length / symsPerElement +
        ((length % symsPerElement == 0) ? 0 : 1)
      ];
      this.mask = calcMask(packing);
      wordsize = packing.wordSize();

      // pack the body of the sequence
      int ii = 0;
      for(int i = 0; i < (syms.length - 1); i++) {
        long l = 0;
        int jj = 0;
        for(int j = 0; j < symsPerElement; j++) {
          long p = packing.pack(symbols[ii + j]);
          l |= (long) ((long) p << (long) jj);
          jj += wordsize;
        }
        syms[i] = l;
        ii += symsPerElement;
      }

      // pack the final word
      if(syms.length > 0) {
        long l = 0;
        ii = (syms.length - 1) * symsPerElement;
        int jMax = length % symsPerElement;
        if(jMax == 0) {
          jMax = symsPerElement;
        }
        for(int j = 0; j < jMax; j++) {
          int jj = j * packing.wordSize();
          long p = packing.pack(symbols[ii + j]);
          l |= (long) ((long) p << (long) jj);
        }
        syms[syms.length - 1] = l;
      }
    } catch (IllegalSymbolException ise) {
      throw new AssertionFailure("Assertion Failure: Symbol got lost somewhere",ise);
    }
  }

  public Symbol symbolAt(int indx) {
    indx--;

    int word;
    int offset;
    long l;

    synchronized(this) {
      if ((indx < currentMin) || (indx > currentMax)) {
        word = indx / symsPerElement;
        offset = indx % symsPerElement;

        currentMin = indx - offset;
        currentMax = currentMin + symsPerElement - 1;
        currentWord = syms[word];
      }
      else {
        offset = indx - currentMin;
      }

      l = currentWord;
    }

    int jj = offset * wordsize;
    try {
      return packing.unpack((byte) ((l >> (long) jj) & mask));
    } catch (IllegalSymbolException ise) {
      throw new AssertionFailure("Could not unpack " + indx + " at " + "word" + "," + offset, ise);
    }
  }


  private static byte calcMask(Packing packing) {
    byte mask = 0;
    for(int i = 0; i < packing.wordSize(); i++) {
      mask |= 1 << i;
    }
    return mask;
  }

  /**
   * <p>
   * Return the long array within which the symbols are bit-packed.
   * </p>
   *
   * <p>
   * <em>Warning:</em> This is a risky developer method.
   * This is the actual array that this object uses to store the bits
   * representing symbols. You should not modify this in any way. If you do,
   * you will modify the symbols returned by symbolAt(). This methd is
   * provided primarily as an easy way for developers to extract the
   * bit pattern for storage in such a way as it could be fetched later and
   * fed into the appropriate constructor.
   * </p>
   *
   * @return the actual long array used to store bit-packed symbols
   */
  public long[] getSyms() {
    return syms;
  }
}
