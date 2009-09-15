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
import org.biojava.bio.seq.db.SequenceDB;
import org.biojava.bio.symbol.IllegalSymbolException;

/**
 * @author Matthew Pocock
 * @author Thomas Down
 */
public interface TrainingAlgorithm {
  DP getDP();
  double getLastScore();
  double getCurrentScore();
  int getCycle();

  /**
   * Trains the sequences in db untill stopper says to finnish.
   */
  void train(SequenceDB db,
             double nullWeight, StoppingCriteria stopper)
  throws IllegalSymbolException, BioException;
}
