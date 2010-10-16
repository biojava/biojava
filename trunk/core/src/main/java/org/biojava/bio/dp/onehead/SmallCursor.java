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

import org.biojava.bio.dp.State;
import org.biojava.bio.symbol.SymbolList;

/**
 * Constant-memory implementation of single-head DP cursor.
 *
 * @author Matthew Pocock
 */
public class SmallCursor extends AbstractCursor {
  private final SymbolList symList;
  private double [] currentC;
  private double [] lastC;
  
  public SymbolList symList() {
    return symList;
  }
  
  public int length() {
    return symList.length();
  }
  
  public double [] currentCol() {
    return currentC;
  }
  
  public double [] lastCol() {
    return lastC;
  }

  public void advance() {
    super.advance();
    
    double [] v = lastC;
    lastC = currentC;
    currentC = v;
  }
  
  public SmallCursor(State [] states, SymbolList symList, Iterator symIterator) {
    super(symIterator);
    this.symList = symList;
    
    this.currentC = new double[states.length];
    this.lastC = new double[states.length];
  }
}
