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

import java.io.Serializable;

import org.biojava.bio.dp.DP;
import org.biojava.bio.dp.DPMatrix;
import org.biojava.bio.dp.MarkovModel;
import org.biojava.bio.dp.State;
import org.biojava.bio.symbol.SymbolList;

/**
 * The dynamic programming matrix for a single sequence.
 *
 * @author Matthew Pocock
 * @author Lukas Kall
 */
public class SingleDPMatrix implements DPMatrix, Serializable {
  protected final State [] states;
  protected final MarkovModel model;
  protected final SymbolList [] symList;
  public final double [][] scores; // [symbol][state]
  protected double score;
 
  public State [] states() {
    return states;
  }
  
  public MarkovModel model() {
    return model;
  }
  
  public SymbolList [] symList() {
    return symList;
  }
  
  public double getScore() {
    return score;
  }
  
  public void setScore(double score) {
    this.score = score;
  }
  
  public double getCell(int [] index)
  throws IndexOutOfBoundsException {
    if(index.length != 2) {
      throw new IndexOutOfBoundsException("index must be two-dimensional");
    }
    return scores[index[1]][index[0]];
  }
  
  public SingleDPMatrix(DP dp, SymbolList symList) {
    this.model = dp.getModel();
    this.states = dp.getStates();
    this.symList = new SymbolList [] { symList };
    this.score = Double.NaN;
    this.scores = new double[symList.length() + 2][states.length];
  }
}
