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


package org.biojava.bio.dist;

import org.biojava.bio.symbol.AtomicSymbol;
import org.biojava.bio.symbol.IllegalSymbolException;

/**
 * A distribution trainer that just ignores all counts.
 *
 * @author Matthew Pocock
 * @since 1.0
 */
public class IgnoreCountsTrainer implements DistributionTrainer {
  public void addCount(DistributionTrainerContext dtc, AtomicSymbol sym, double times) 
    throws IllegalSymbolException {}
  public double getCount(DistributionTrainerContext dtc, AtomicSymbol sym)
  throws IllegalSymbolException {
    return 0.0;
  }
  
  public void train(DistributionTrainerContext dtc, double weight) {}
  public void clearCounts(DistributionTrainerContext dtc) {}

  /**
   * Constructor intended for sub-classes.
   */
  protected IgnoreCountsTrainer() {}
  
  private static IgnoreCountsTrainer instance = new IgnoreCountsTrainer();

  /**
   * Returns the global singleton instance of the IgnoreCountsTrainer.
   *
   * @return the singleton instance of this trainer
   */
  public static IgnoreCountsTrainer getInstance() {
    return instance;
  }
}
