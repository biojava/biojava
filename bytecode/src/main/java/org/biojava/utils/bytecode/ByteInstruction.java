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
 * Instructions which take a one-byte operand.
 *
 * @author Thomas Down
 * @author Matthew Pocock
 */

class ByteInstruction implements Instruction {
  private final byte opcode;
  private final byte val;
  private final int delta;
  
  ByteInstruction(byte opcode, byte val, int delta) {
    this.opcode = opcode;
    this.val = val;
    this.delta = delta;
  }
  
  public void writeCode(CodeContext ctx) throws CodeException {
    ctx.writeByte(opcode);
    ctx.writeByte(val);
  }
  
  public int stackDepth() {
    return Math.max(delta, 0);
  }
  
  public int stackDelta() {
    return delta;
  }
}
