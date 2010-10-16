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

import org.biojava.bio.dist.DistributionTrainerContext;

/**
 * Encapsulates the training of an entire model.
 *
 * @author Matthew Pocock
 */
public interface ModelTrainer extends DistributionTrainerContext {
  /**
   * Registers an HMM with this trainer.
   * <p>
   * If the model has been already registered, then this method should do
   * nothing. If it has not been registered, then this method should loop over
   * every state in the model and register the Distribution returned by
   * getWeight.
   *
   * @param model the MarkovModel to train
   */
  void registerModel(MarkovModel model);
}
