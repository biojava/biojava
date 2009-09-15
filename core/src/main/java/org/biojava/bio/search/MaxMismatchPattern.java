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

import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.SymbolList;

/**
 * This class permits searching a SymbolList with another SymbolList while
 * permitting a specified number of mismatches.  The search pattern can
 * include ambiguity Symbols.
 *
 * @author Matthew Pocock (wrote original MaxMissmatchPattern class)
 * @author David Huen (debugging and extension to permit ambiguity symbols)
 */
public class MaxMismatchPattern
implements BioPattern {
  private int mismatches;
  private SymbolList pattern;

  public MaxMismatchPattern() {}

  public MaxMismatchPattern(SymbolList pattern, int mismatches) {
    this.pattern = pattern;
    this.mismatches = mismatches;
  }

  public int getMismatches() {
    return mismatches;
  }

  public void setMismatches(int mismatches) {
    this.mismatches = mismatches;
  }

  public SymbolList getPattern() {
    return pattern;
  }

  public void setPattern(SymbolList pattern) {
    this.pattern = pattern;
  }

  public BioMatcher matcher(SymbolList symList)
          throws IllegalAlphabetException {
    return new MaxMismatchMatcher(symList, pattern, mismatches);
  }
}


