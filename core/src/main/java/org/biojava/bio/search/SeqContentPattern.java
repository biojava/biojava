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

package org.biojava.bio.search;

import org.biojava.bio.symbol.AlphabetIndex;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.AtomicSymbol;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.SymbolList;

/**
 * A pattern that can be used to find regions with given sequence content.
 *
 * <p>
 * Regular expressions can be used to find sequence patterns. However, some
 * things can't be easily expressed as a regular expression. For example,
 * a region of length 10 that contains at least 8 Gs and up to two Ts and no
 * other symbols. A SeqContentPattern can be built that does represent this.
 * <p>
 *
 * <code><pre>
 * SeqContentPattern scp = new SeqContentPattern(DNATools.getDNA());
 * scp.setLength(10);
 * scp.setMinCounts(DNATools.g(), 8);
 * scp.setMaxCounts(DNATools.t(), 2);
 * scp.setMaxCounts(DNATools.c(), 0);
 * scp.setMaxCounts(DNATools.a(), 0);
 * </pre></code>
 *
 * <p>
 * The minimum counts default to 0, and the maximum counts default to the
 * length. If you have not manually set the maximum count for a symbol, it will
 * continue to adjust while you change the length. Once you have set it, it will
 * not vary, even if you do set the length. To re-set a maximum count to track
 * the length, set it to -1.
 * </p>
 *
 * <p>
 * All regions of the defined length for which all constraints are satisfied
 * will potentialy be found. At the moment we have not defined what will
 * happen for multiple regions that overlap all of which satisfy the
 * constraints.
 * </p>
 *
 * @author Matthew Pocock
 * @since 1.4
 */
public class SeqContentPattern implements BioPattern {
  private final AlphabetIndex index;
  private final int[] minCounts;
  private final int[] maxCounts;
  private int length;

  /**
   * Create a new SeqContentPattern over an alphabet.
   *
   * @param alpha  the FiniteAlphabet for this pattern
   */
  public SeqContentPattern(FiniteAlphabet alpha) {
    index = AlphabetManager.getAlphabetIndex(alpha);
    this.minCounts = new int[alpha.size()];
    this.maxCounts = new int[alpha.size()];

    for(int i = 0; i < minCounts.length; i++) {
      minCounts[i] = 0;
      maxCounts[i] = -1;
    }
  }

  /**
   * Get the current length.
   *
   * @return the length
   */
  public int getLength() {
    return length;
  }

  /**
   * Set the pattern length.
   *
   * @param length  the new length
   */
  public void setLength(int length) {
    this.length = length;
  }

  /**
   * Set the minimum counts required for a symbol.
   *
   * @param as  the AtomicSymbol to check
   * @param count  the minimum number of counts it must have
   * @throws IllegalSymbolException  if as is not known in this alphabet
   */
  public void setMinCounts(AtomicSymbol as, int count)
  throws IllegalSymbolException {
    minCounts[index.indexForSymbol(as)] = count;
  }

  /**
   * Get the minimum counts required for a symbol.
   *
   * @param as  the AtomicSymbol to check
   * @return the minimum number of counts it must have
   * @throws IllegalSymbolException  if as is not known in this alphabet
   */
  public int getMinCounts(AtomicSymbol as)
  throws IllegalSymbolException {
    return minCounts[index.indexForSymbol(as)];
  }

  /**
   * Set the maximum counts required for a symbol.
   * Use -1 to reset it to track the length.
   *
   * @param as  the AtomicSymbol to check
   * @param count  the maximum number of counts it must have
   * @throws IllegalSymbolException  if as is not known in this alphabet
   */
  public void setMaxCounts(AtomicSymbol as, int count)
  throws IllegalSymbolException {
    maxCounts[index.indexForSymbol(as)] = count;
  }

  /**
   * Get the maximum counts required for a symbol.
   *
   * @param as  the AtomicSymbol to check
   * @return the maximum number of counts it must have
   * @throws IllegalSymbolException  if as is not known in this alphabet
   */
  public int getMaxCounts(AtomicSymbol as)
  throws IllegalSymbolException {
    int c = maxCounts[index.indexForSymbol(as)];
    if(c == -1) {
      return length;
    } else {
      return c;
    }
  }

  public BioMatcher matcher(SymbolList symList)
  throws IllegalAlphabetException {
    if(symList.getAlphabet() != index.getAlphabet()) {
      throw new IllegalAlphabetException(
        "Attempted to search symbol list over " + symList.getAlphabet() +
        " but the search parameters only accept " + index.getAlphabet() );
    }

    int[] minCounts = new int[this.minCounts.length];
    int[] maxCounts = new int[this.maxCounts.length];
    for(int i = 0; i < minCounts.length; i++) {
      minCounts[i] = this.minCounts[i];

      int c = this.maxCounts[i];
      maxCounts[i] = (c == -1) ? length : c;
    }

    return new SeqContentMatcher(
      symList,
      index,
      minCounts,
      maxCounts,
      length );
  }
}

