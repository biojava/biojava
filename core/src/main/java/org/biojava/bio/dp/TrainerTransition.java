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

/**
 * This is a small and ugly class for storing a trainer and a transition.
 * <p>
 * It is hash-code-able, and has a sensible equality operator.
 *
 * @author Matthew Pocock
 */
public class TrainerTransition implements Serializable {
  public TransitionTrainer trainer;
  public State from;
  public State to;
  
  /**
   * Two transitions are equal if they have the same trainer, from and to states.
   */
  public boolean equals(Object o)
  throws ClassCastException {
    TrainerTransition t = (TrainerTransition) o;
    return trainer == t.trainer && from == t.from && to == t.to;
  }
  
  /**
   * The hash code is model.hashCode() ^ from.hashCode() ^ to.hashCode()
   */
  public int hashCode() {
    return trainer.hashCode() ^ from.hashCode() ^ to.hashCode();
  }
  
  public TrainerTransition(TransitionTrainer trainer, State from, State to) {
    this.trainer = trainer;
    this.from = from;
    this.to = to;
  }
}
