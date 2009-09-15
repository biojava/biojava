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

import org.biojava.bio.BioException;

/**
 * Flags an object as being able to register itself with a model trainer.
 *
 * @author Matthew Pocock
 */
public interface Trainable {
  /**
   * Perform any registration that is necessary with mt.
   * <p>
   * This may include registering handlers for transition or emission counts,
   * or registering other Trainable objects with the ModelTrainer.
   *
   * @param mt  the ModelTrainer that encapsulates the training environment
   */
  public void registerWithTrainer(ModelTrainer mt)
  throws BioException;
}
