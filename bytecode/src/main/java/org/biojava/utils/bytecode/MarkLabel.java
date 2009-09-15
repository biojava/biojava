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
package org.biojava.utils.bytecode;

import java.util.*;

/**
 * A CodeGenerator that just marks a label that can be used for jumps.
 *
 * @author Matthew Pocock
 */

public class MarkLabel implements CodeGenerator {
  private final Label label;
  
  public MarkLabel(Label label) {
    this.label = label;
  }
  
  public void writeCode(CodeContext ctx) throws CodeException {
    ctx.markLabel(label);
  }
  
  public int stackDepth() {
    return 0;
  }
  
  public int stackDelta() {
    return 0;
  }
}

