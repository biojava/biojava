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
 * This exception indicates that there is no transition between two states.
 *
 * @author Matthew Pocock
 */
public class IllegalTransitionException extends BioException {
  private State from;
  private State to;
  
  public State getFrom() {
    return from;
  }
    
  public State getTo() {
    return to;
  }
  
  public IllegalTransitionException(State from, State to, String message) {
    super(message + "[" + from.getName() + " -> " + to.getName() + "]");
    this.from = from;
    this.to = to;
  }

  public IllegalTransitionException(State from, State to) {
    this(from, to, "");
  }
  
  public IllegalTransitionException() {
    super();
  }
}
