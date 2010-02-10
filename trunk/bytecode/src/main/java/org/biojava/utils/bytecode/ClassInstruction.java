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
 * Instructions which use a Class.
 *
 * @author Matthew Pocock
 */

class ClassInstruction implements Instruction {
  private final CodeClass clazz;
  private final byte opcode;
  private final int delta;
  
  ClassInstruction(byte op, CodeClass c, int delta) {
    if(c == null) {
      throw new NullPointerException("CodeClass can not be null");
    }
    this.opcode = op;
    this.clazz = c;
    this.delta = delta;
  }
  
  public void writeCode(CodeContext ctx) throws CodeException {
    ctx.writeByte(opcode);
    ctx.writeShort(ctx.getConstants().resolveClass(clazz));
  }
  
  public int stackDepth() {
    return Math.max(delta, 0);
  }
  
  public int stackDelta() {
    return delta;
  }
}
