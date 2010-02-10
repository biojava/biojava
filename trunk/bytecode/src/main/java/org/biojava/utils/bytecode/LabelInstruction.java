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

/**
 * Instructions which jump to a label.
 *
 * @author Thomas Down
 * @author Matthew Pocock
 */

class LabelInstruction implements Instruction {
  private final Label l;
  private final byte opcode;
  private final int delta;
  
  LabelInstruction(byte op, Label l, int delta) {
    if(l == null) {
      throw new NullPointerException("Label can not be null");
    }
    
    this.opcode = op;
    this.l = l;
    this.delta = delta;
  }
  
  public void writeCode(CodeContext ctx) throws CodeException {
    ctx.writeByte(opcode);
    ctx.writeLabel(l);
  }
  
  public int stackDepth() {
    return 0;
  }
  
  public int stackDelta() {
    return delta;
  }
}
