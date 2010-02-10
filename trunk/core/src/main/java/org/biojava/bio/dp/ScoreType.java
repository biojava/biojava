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


package org.biojava.bio.dp;
import java.io.Serializable;

import org.biojava.bio.dist.Distribution;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;

/**
 * This class computes the score that is used to be used
 * in a DP optimisation.
 *
 * @author Matthew Pocock
 */
public interface ScoreType {
  /**
   * Calculates the score associated with a distribution and a symbol.
   */
  double calculateScore(Distribution dist, Symbol sym)
  throws IllegalSymbolException;
  
  public final static ScoreType PROBABILITY = new Probability();
  
  /**
   * In this class, calculateScore returns the probability
   * of a Symbol being emitted.
   *
   * @author Matthew Pocock
   */
  public static class Probability implements ScoreType, Serializable {
    public double calculateScore(Distribution dist, Symbol sym)
    throws IllegalSymbolException {
      return dist.getWeight(sym);
    }
  };
  
  public final static ScoreType ODDS = new Odds();
  
  /**
   * In this class, calculateScore returns the odds ratio
   * of a symbol being emitted.  That is, the ratio of the
   * probability of a Symbol being emitted to it being
   * emitted by the null model.
   *
   * @author Matthew Pocock
   */
  public static class Odds implements ScoreType, Serializable {
    public double calculateScore(Distribution dist, Symbol sym)
    throws IllegalSymbolException {
      double d = dist.getWeight(sym);
      double n = dist.getNullModel().getWeight(sym);
      //System.out.println("Odds for " + sym.getName() + "\t= " + d + " / " + n);
      return d / n;
    }
  };
  
  public final static ScoreType NULL_MODEL = new NullModel();
  
  /**
   * In this class, calculateScore returns the probability of
   * a Symbol being emitted by the null model.
   *
   * @author Matthew Pocock
   */
  public static class NullModel implements ScoreType, Serializable {
    public double calculateScore(Distribution dist, Symbol sym)
    throws IllegalSymbolException {
      return dist.getNullModel().getWeight(sym);
    }
  };
}

