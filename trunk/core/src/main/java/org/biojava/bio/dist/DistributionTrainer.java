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
import org.biojava.utils.ChangeVetoException;

/**
 * <p>
 * An object that can be used to train a distribution up.
 * </p>
 *
 * <p>
 * This lets the distribution implementation handle counts or distributions
 * in the best way possible.
 * </p>
 *
 * @author Matthew Pocock
 * @since 1.0
 */
public interface DistributionTrainer {
  /**
   * <p>
   * Registers that sym was counted in this state.
   * </p>
   *
   * <p>
   * This method may be called multiple times with the same symbol. In this
   * case, the times should be summed.
   * </p>
   *
   * @param dtc the DistributionTrainerContext within which the count was added
   * @param sym the Symbol seen
   * @param times the number of times to add
   * @throws IllegalSymbolException if sym is not recognised
   */
  void addCount(DistributionTrainerContext dtc, AtomicSymbol sym, double times)
  throws IllegalSymbolException;
  
  /**
   * <p>
   * Get the current count for this state.
   * </p>
   *
   * <p>
   * This method may be called multiple times with the same symbol. Each time
   * it should return the agregate of the counts added with addCount since the
   * last invocation of clearCounts.
   * </p>
   *
   * @param dtc the DistributionTrainerContext within which the count was added
   * @param sym the Symbol seen
   * @return the agregate of the counts
   * @throws IllegalSymbolException  if sym is not recognised
   */
  double getCount(DistributionTrainerContext dtc, AtomicSymbol sym)
  throws IllegalSymbolException;
  
  /**
   * <p>
   * Trains the Distribution, given a null model.
   * </p>
   *
   * <p>
   * This will use the information collected with multiple addCount calls, and
   * the null model to generate the new weights.
   * </p>
   *
   * <p>
   * This method should not modify the underlying counts.
   * </p>
   *
   * @param dtc     the context to use
   * @param weight  how many lots of the null model to add
   * @throws ChangeVetoException if the distribution could not have its weights
   *         modified
   */
  void train(DistributionTrainerContext dtc, double weight)
  throws ChangeVetoException;
  
  /**
   * Clears all of the counts to zero.
   */
  void clearCounts(DistributionTrainerContext dtc);
}
