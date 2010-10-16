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
 *
 *
 * @author Matthew Pocock
 */
public interface BioPattern {
  /**
   * Get a matcher that will use these parameters to search a SymbolList.
   *
   * <p>
   * The resulting BioMatcher is independant of this BioPattern.
   * In particular, calling any mutator methods on this pattern will not affect
   * the matcher.
   * </p>
   *
   * @param symList  the SymbolList to match against
   * @return a BioMatcher that will perform the search
   * @throws IllegalAlphabetException if symList is not over the right alphabet
   */
  BioMatcher matcher(SymbolList symList)
  throws IllegalAlphabetException;
}
