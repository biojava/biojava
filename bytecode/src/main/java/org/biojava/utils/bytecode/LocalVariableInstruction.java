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
 * Instructions which load or store to a local variable.
 *
 * @author Thomas Down
 * @author Matthew Pocock
 */

class LocalVariableInstruction implements Instruction {
  private final LocalVariable var;
  private final byte opcode;
  private final byte specialOpCodeBase;
  
  LocalVariableInstruction(byte op, byte special, LocalVariable var) {
    if(var == null) {
      throw new NullPointerException("LocalVariable can not be null");
    }
    
    this.opcode = op;
    specialOpCodeBase = special;
    this.var = var;
  }
  
  public void writeCode(CodeContext ctx) throws CodeException {
    int slot = ctx.resolveLocal(var);
    if (slot <= 3) {
      ctx.writeByte((byte) (specialOpCodeBase + slot));
    } else if (slot < 256) {
      ctx.writeByte(opcode);
      ctx.writeByte((byte) slot);
    } else {
      ctx.writeByte(ByteCode.op_wide);
      ctx.writeByte(opcode);
      ctx.writeShort(slot);
    }
  }
  
  public int stackDepth() {
    return stackDelta();
  }
  
  public int stackDelta() {
    return 1;
  }
}
