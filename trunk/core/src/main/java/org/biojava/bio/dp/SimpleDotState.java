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

import org.biojava.bio.Annotation;
import org.biojava.bio.symbol.FundamentalAtomicSymbol;

/**
 * A Dot state that you can make and use.
 * <p>
 * Dot states emit no sequence. They are there purely to make the wireing
 * of the model look neater, and to cut down the number of combinatorial
 * transitions that can so easily swamp models.
 *
 * @author Matthew Pocock
 * @author Thomas Down
 */
public class SimpleDotState
extends FundamentalAtomicSymbol implements DotState, Serializable {
  /**
   * Construct a new state with the specified name and annotation.
   * The token parameter is ignored but included for compatibility.
   *
   * @deprecated token is ignored since 1.2.  Use the 2-arg constructor instead.
   */
    
  public SimpleDotState(char token, String name, Annotation annotation) {
    super(name, annotation);
  }
  
  /**
   * Construct a new state with the specified name and annotation
   */
    
  public SimpleDotState(String name, Annotation annotation) {
    super(name, annotation);
  }
  
  public SimpleDotState(String name) {
    super(name, Annotation.EMPTY_ANNOTATION);
  }
}
