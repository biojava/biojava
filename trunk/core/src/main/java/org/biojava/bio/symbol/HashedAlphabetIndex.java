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

import java.lang.ref.Reference;
import java.lang.ref.WeakReference;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Iterator;

import org.biojava.bio.BioError;
import org.biojava.utils.AbstractChangeable;
import org.biojava.utils.ChangeListener;

/**
 * Uses Arrays.binarySearch to retrieve indecies for symbols. To save on CPU,
 * an array of symbol hash codes is searched, avoiding the need to multipuly
 * calculate the hash codes of the alphabet symbols.
 *
 * @author Matthew Pocock
 * @since 1.1
 */
class HashedAlphabetIndex
extends AbstractChangeable implements AlphabetIndex, java.io.Serializable {
  private static final Comparator cmp = new HashComparator();

  private final Reference alphaRef;
  private final Symbol[] symbols;
  private final int[] hashes;

  public FiniteAlphabet getAlphabet() {
    return (FiniteAlphabet) alphaRef.get();
  }

  public int indexForSymbol(Symbol s)
  throws IllegalSymbolException {
    int indx = Arrays.binarySearch(hashes, s.hashCode());
    if(indx < 0) {
      getAlphabet().validate(s);
      throw new BioError(
        "Assertion Failure: " +
        "Symbol " + s.getName() + " was not an indexed member of the alphabet " +
        getAlphabet().getName() + " despite being in the alphabet."
      );
    }

    // we hit the correct symbol first time
    if(symbols[indx] == s) {
      return indx;
    }

    // it may have the same hash code and be after
    for(
      int i = indx;
      i < symbols.length && hashes[i] == hashes[indx];
      i++
    ) {
      if(symbols[i].equals(s)) {
        return i;
      }
    }

    // in some strange parallel universe, it may have the same hashcode and
    // be before
    for(
      int i = indx-1;
      i >= 0 && hashes[i] == hashes[indx];
      i--
    ) {
      if(symbols[i].equals(s)) {
        return i;
      }
    }

    // it has the same hash code, but isn't in the alphabet
    getAlphabet().validate(s);
    if(s instanceof AtomicSymbol) {
      throw new BioError(
        "Assertion Failure: " +
        "Symbol " + s.getName() + " was not an indexed member of the alphabet " +
        getAlphabet().getName() + " despite being in the alphabet."
      );
    } else {
      throw new IllegalSymbolException("Symbol must be atomic to be indexed.");
    }
  }

  public Symbol symbolForIndex(int i) throws IndexOutOfBoundsException {
    return symbols[i];
  }

  public HashedAlphabetIndex(FiniteAlphabet alpha) {
    alpha.addChangeListener(ChangeListener.ALWAYS_VETO, Alphabet.SYMBOLS);
    this.alphaRef = new WeakReference(alpha);
    symbols = new Symbol[alpha.size()];
    hashes = new int[alpha.size()];

    int i = 0;
    Iterator s = alpha.iterator();
    while(s.hasNext()) {
      symbols[i++] = (Symbol) s.next();
    }
    Arrays.sort(symbols, cmp);

    for(i = 0; i < symbols.length; i++) {
      hashes[i] = symbols[i].hashCode();
    }
  }

  private static class HashComparator implements Comparator {
    public boolean equals(Object o) {
      return o instanceof HashComparator;
    }

    public int compare(Object a, Object b) {
      return a.hashCode() - b.hashCode();
    }
  }
}
