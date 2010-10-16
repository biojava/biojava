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

import java.util.HashMap;
import java.util.Map;

import org.biojava.bio.seq.DNATools;

/**
 * <p>
 * A factory that is used to maintain associations between alphabets and
 * preferred bit-packings for them.
 * </p>
 *
 * <p>
 * There are many ways to pack the symbols for an alphabet as binary.
 * Different applications will wish to have different representations for
 * reasons such as integration with external formats, whether to store
 * ambiguity or not and what algorithms may be used on the bit-packed
 * representation. Also, it has utility methods to arrange the bit-strings
 * for symbols within a Java int primitive.
 * </p>
 *
 * @author Matthew Pocock
 * @author Thomas Down
 */
public class PackingFactory {
  private final static Map packForAlpha;

  static {
    packForAlpha = new HashMap();

  }

  /**
   * Get the default packing for an alphabet.
   *
   * @param alpha  the FiniteAlphabet that will be bit-packed
   * @param ambiguity  true if the packing should store ambiguity and false
   *                   if it can discard ambiguity information
   * @throws IllegalAlphabetException if this combination of alphabet and
   *                   ambiguity can not be handled
   **/
  public static Packing getPacking(FiniteAlphabet alpha, boolean ambiguity)
  throws IllegalAlphabetException {
    Packing pack = (Packing) packForAlpha.get(alpha);
    if(pack == null) {
      if(alpha == DNATools.getDNA()) {
        if(ambiguity) {
          pack = new DNAAmbPack();
        } else {
          pack = new DNANoAmbPack(DNATools.a());
        }
      } else {
        if(ambiguity) {
          // can't handle this in the general case
          throw new IllegalAlphabetException();
        } else {
          pack = new IndexedNoAmbPack(AlphabetManager.getAlphabetIndex(alpha));
        }
      }
    }
    return pack;
  }

  public static int primeWord(
    SymbolList symList,
    int wordLength,
    Packing packing
  ) throws IllegalSymbolException {
    int word = 0;
    for(int i = 0; i < wordLength; i++) {
      int p = packing.pack(symList.symbolAt(i+1));
      word |= (int) ((int) p << (int) (i * packing.wordSize()));
    }
    return word;
  }

  public static int nextWord(
    SymbolList symList,
    int word,
    int offset,
    int wordLength,
    Packing packing
  ) throws IllegalSymbolException {
    word = word >> (int) packing.wordSize();
    int p = packing.pack(symList.symbolAt(offset));
    word |= (int) p << ((int) (wordLength - 1) * packing.wordSize());
    return word;
  }

  public static void binary(long val) {
    for(int i = 63; i >= 0; i--) {
      System.out.print( ((((val >> i) & 1) == 1) ? 1 : 0) );
    }
    System.out.println();
  }
  public static void binary(int val) {
    for(int i = 31; i >= 0; i--) {
      System.out.print( ((((val >> i) & 1) == 1) ? 1 : 0) );
    }
    System.out.println();
  }
}

