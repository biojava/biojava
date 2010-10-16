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
 * Instructions which act on fields.
 *
 * @author Thomas Down
 * @author Matthew Pocock
 */

class FieldInstruction implements Instruction {
  private final CodeField field;
  private final byte opcode;
  
  FieldInstruction(byte op, CodeField f) {
    if(f == null) {
      throw new NullPointerException("CodeField can not be null");
    }
    this.opcode = op;
    this.field = f;
  }
  
  public void writeCode(CodeContext ctx) throws CodeException {
    ctx.writeByte(opcode);
    ctx.writeShort(ctx.getConstants().resolveField(field));
  }
  
  public int stackDepth() {
    return stackDelta();
  }
  
  public int stackDelta() {
    return 1;
  }
}
