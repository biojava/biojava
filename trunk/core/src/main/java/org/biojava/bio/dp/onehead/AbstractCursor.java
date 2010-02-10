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


package org.biojava.bio.dp.onehead;

import java.util.Iterator;

import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.Symbol;

/**
 * An abstract instance of a single-head DP cursor.
 *
 * @author Matthew Pocock
 */
abstract class AbstractCursor implements DPCursor {
  private static final Symbol gap = AlphabetManager.getGapSymbol();
  private Iterator symIterator;
  
  private Symbol currentRes;
  private Symbol lastRes;
  
  public Symbol currentRes() {
    return currentRes;
  }
  
  public Symbol lastRes() {
    return lastRes;
  }
  
  public boolean canAdvance() {
    return symIterator.hasNext() || currentRes != gap;
  }
  
  public void advance() {
    lastRes = currentRes;
    currentRes = (symIterator.hasNext()) ? (Symbol) symIterator.next()
                                         : gap;
  }
  
  public AbstractCursor(Iterator symIterator) {
    this.symIterator = symIterator;
    this.currentRes = gap;
    this.lastRes = gap;
  }
  
  protected AbstractCursor() {
  }
}
