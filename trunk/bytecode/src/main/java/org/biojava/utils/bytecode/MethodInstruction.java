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
 * Instructions which call a method.
 *
 * @author Thomas Down
 * @author Matthew Pocock
 */

class MethodInstruction implements Instruction {
  private final CodeMethod meth;
  private final byte opcode;
  
  MethodInstruction(byte op, CodeMethod m) {
    if(m == null) {
      throw new NullPointerException("CodeMethod can not be null");
    }
    
    this.opcode = op;
    this.meth = m;
  }
  
  public void writeCode(CodeContext ctx) throws CodeException {
    ctx.writeByte(opcode);
    if (opcode == ByteCode.op_invokeinterface) {
      ctx.writeShort(ctx.getConstants().resolveInterfaceMethod(meth));
      int count = 1;
      for (int i = 0; i < meth.numParameters(); ++i) {
        CodeClass ptype = meth.getParameterType(i);
        count += CodeUtils.wordsForType(ptype);
      }
      ctx.writeByte((byte) count);
      ctx.writeByte((byte) 0);
    } else {
      ctx.writeShort(ctx.getConstants().resolveMethod(meth));
    }
  }
  
  public int stackDepth() {
    return 0;
  }
  
  public int stackDelta() {
    int popped = 0;
    
    if( (meth.getModifiers() & CodeUtils.ACC_STATIC) == 0) {
      popped++;
    }
    
    popped += meth.numParameters();
    
    int pushed = (meth.getReturnType() == CodeUtils.TYPE_VOID) ? 0 : 1;
    
    return pushed - popped;
  }
}
