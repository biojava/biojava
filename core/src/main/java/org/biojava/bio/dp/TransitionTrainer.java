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

import org.biojava.bio.symbol.IllegalSymbolException;

/**
 * An object that can be used to train the transitions within a MarkovModel.
 *
 * @author Matthew Pocock
 */
public interface TransitionTrainer {
  /**
   * Add 'count' to the transition from->to.
   * <p>
   * This method may be called multiple times with the same from,to pair in
   * which case, the times should be summed.
   *
   * @param from  the source state
   * @param to  the destination state
   * @param count   the number of counts to add
   */
  void addCount(State from, State to, double count)
  throws IllegalSymbolException, IllegalTransitionException;
  
  /**
   * Trains the transition, given an expected probability, and a weight for
   * that probability.
   * <p>
   * This is equivalent to adding a count of nullModel * weight to each
   * transition and then training with a weight of 0.
   *
   * @param nullModel the nullModel to use
   * @param weight  how many lots of the null model to add
   */
  void train(double nullModel, double weight) throws IllegalSymbolException;
  
  /**
   * Clears all of the counts to zero.
   */
  void clearCounts();
}
