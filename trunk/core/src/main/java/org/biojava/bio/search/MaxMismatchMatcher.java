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

import org.biojava.bio.symbol.BasisSymbol;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;

/**
 * A BioMatcher class returned by MaxMismatchPattern.matcher() that implements
 * searching of a SymbolList.
 * <p>
 * This class is public only to allow access to the mismatchCount() method.
 *
 * @author Matthew Pocock (wrote original MaxMissmatchMatcher class)
 * @author David Huen (debugging and extension of functionality)
 */
public class MaxMismatchMatcher
implements BioMatcher {
    // primary data
    private final Symbol [] patternSymbol; // an array containing the pattern symbols in reversed order
    private final SymbolList seq;

    // precomputed constants
    private final int patLength;
    private final int maxPatternSymbolIdx;
    private final int seqLength;
    private final int minMatches;
    private final boolean [] ambiguousPosition; // indicates if the corresponding position in patternSymbol represents
                                                // an ambiguity
    private boolean hasAmbiguity =false;

    // working numbers
    private final int[] matches; // this is a modulo'd rotating register that keeps track of the match count
    private int pos; // this is the symbol in sequence that will be added for consideration

    MaxMismatchMatcher(SymbolList seq,
                      SymbolList pattern,
                      int mismatches)
    {
        this.seq = seq;

        patLength = pattern.length();
        maxPatternSymbolIdx = patLength - 1;
        seqLength = seq.length();
        minMatches = patLength - mismatches;

        // construct reversed pattern array
        patternSymbol = new Symbol[patLength];
        ambiguousPosition = new boolean[patLength];
        for (int i=0; i < patLength; i++) {
            patternSymbol[i] = pattern.symbolAt(patLength - i);
            if (patternSymbol[i] instanceof BasisSymbol) {
                ambiguousPosition[i] = true;
                hasAmbiguity = true;
            }
        }

        // initialize matches
        matches = new int[patLength];
        for(int i = 0; i < matches.length; i++) matches[i] = 0;

        pos = 1;
    }

    public boolean find() 
    {

        if (pos >= seq.length()) {
            return false;
        }

        for (; pos <= seqLength; pos++) {
            if (addSymbol()) {
                pos++;
                return true;
            }
        }

        return false;
    }

    // computes consequences of symbol pointed to by pos
    private boolean addSymbol()
    {
        Symbol sym = seq.symbolAt(pos);
        // compute matches for continuation of match
        if (hasAmbiguity) {
            for (int i=0; i < maxPatternSymbolIdx; i++) {
                int idx = (pos + i) % patLength;
                if (ambiguousPosition[i]) {
                    if (patternSymbol[i].getMatches().contains(sym))
                        matches[idx]++;
                }
                else {
                    if (sym == patternSymbol[i])
                        matches[idx]++;
                }
            }
        }
        else {
            for (int i=0; i < maxPatternSymbolIdx; i++) {
                int idx = (pos + i) % patLength;
                if (sym == patternSymbol[i])
                    matches[idx]++;
            }
        }

        // initialise and compute initial match
        matches[(pos + maxPatternSymbolIdx) % patLength] = (sym == patternSymbol[maxPatternSymbolIdx])?1:0;

        return matches[pos % patLength] >= minMatches;
    }

    public int start() {
        // remember that pos is already incremented beyond the match
        return pos - matches.length;
    }

    public int end() {
        // remember that pos is already incremented beyond the match
        return pos - 1;
    }

    public SymbolList group() {
        return seq.subList(start(), end());
    }

    /**
     * Returns number of mismatches
     */
    public int mismatchCount() { return patLength - matches[(pos - 1) % patLength]; }
}


