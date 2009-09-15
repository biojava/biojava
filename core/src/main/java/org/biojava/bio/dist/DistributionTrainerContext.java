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

import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;
import org.biojava.utils.ChangeVetoException;

/**
 * A context within a group of DistributionTrainers can be trained together.
 *
 * @author Matthew Pocock
 * @since 1.0
 */
public interface DistributionTrainerContext {
    /**
     * Return the number of pseudocounts added to the distribution when training.
     *
     * @return the null model weight
     */

  public double getNullModelWeight();
  
    /**
     * Set the number of pseudocounts to add when training the distribution.
     * These counts are added in proportion to the null model of the distribution
     * being trained.
     *
     * @param weight the new null model weight
     */

  public void setNullModelWeight(double weight);
  
  /**
   * <p>
   * Register a distribution object with this context.
   * </p>
   *
   * <p>
   * This method is a request to the context to register dist. If dist is already
   * registered then this method should do nothing. If it is not registered, then
   * it should invoke dist.registerWithTrainer
   * </p>
   *
   * @param dist the Distribution to register
   */
  void registerDistribution(Distribution dist);
  
  /**
   * <p>
   * Register a Distribution and an associated DistributionTrainer object.
   * </p>
   *
   * <p>
   * In the registerWithTrainer method of a Distribution, it should associate
   * itself with a trainer using this method.
   * </p>
   *
   * @param dist the distribution to be registered.
   * @param trainer the distribution's trainer object to be registered.
   */
  void registerTrainer(Distribution dist, DistributionTrainer trainer);
  
  /**
  * Return the Distribution trainer object from the current context.
  *
  * @param dist the Distribution whose trainer is required
  * @return the DistributionTrainer for the distribution
  */
   DistributionTrainer getTrainer(Distribution dist);
  
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
   * @param dist the Distribution that the symbol was associated with
   * @param sym the Symbol seen
   * @param times the number of times to add
   * @throws IllegalSymbolException if sym is not recognised by dist
   */
  void addCount(Distribution dist, Symbol sym, double times)
  throws IllegalSymbolException;
  
  /**
   * Return the number of counts of a particular symbol which will be used
   * to train the specified distribution.
   *
   * @param dist  the Distribution to return counts for
   * @param sym   the symbol to get the count for
   * @return the number of counts
   * @throws IllegalSymbolException if the symbol is not accepted by the
   *     distribution
   */
  
  double getCount(Distribution dist, Symbol sym)
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
   * @throws ChangeVetoException  if any of the distributions can't be trained
   */
  void train() throws ChangeVetoException;
  
  /**
   * Clears all of the counts to zero.
   */
  void clearCounts();
}
