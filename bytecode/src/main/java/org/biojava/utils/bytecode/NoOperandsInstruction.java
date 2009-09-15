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
 * Instructions (like nop, dadd, etc.) which take no operands.
 *
 * @author Thomas Down
 * @author Matthew Pocock
 */

class NoOperandsInstruction implements Instruction {
  private final byte opcode;
  private final int delta;
  
  public NoOperandsInstruction(byte b, int delta) {
    this.opcode = b;
    this.delta = delta;
  }
  
  public void writeCode(CodeContext ctx) throws CodeException {
    ctx.writeByte(opcode);
  }
  
  public int stackDepth() {
    return Math.max(delta, 0);
  }
  
  public int stackDelta() {
    return delta;
  }
}
