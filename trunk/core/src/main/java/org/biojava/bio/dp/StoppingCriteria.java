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

/**
 * A callback that is invoked during the training of an HMM.
 *
 * @author Matthew Pocock
 */
public interface StoppingCriteria {
  /**
   * Decide if the training has completed.
   *
   * This can be on the basis of the scores published by ta. Or, it could be
   * set to expire at a given cycle, or due to user intervention. It is
   * perfectly acceptable for this method to have side effects such as logging.
   *
   * @param ta  the TrainingAlgorithm instance to check
   * @return true if the training is complete, false if it should be continued
   */
  boolean isTrainingComplete(TrainingAlgorithm ta);
}
