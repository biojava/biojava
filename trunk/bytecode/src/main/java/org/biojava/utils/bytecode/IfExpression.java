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
 * A CodeGenerator that provides something semanticaly identical to if.
 * <P>
 * It generates code of the form:
 * <P>
 * conditional branch to trueLabel<br>
 * ifFalse instructions here<br>
 * branch to endLabel<br>
 * trueLabel<br>
 * ifTrue instructions here<br>
 * endLabel
 *
 *
 * @author Matthew Pocock
 */

public class IfExpression implements CodeGenerator {
  private Instruction ifInstruction;
  private CodeGenerator ifTrue;
  private CodeGenerator ifFalse;
  private Label trueLabel;
  private Label endLabel;
  private Instruction skipTrue;
  
  private InstructionVector instructions;
  
  public IfExpression(
    byte ifInstruction, 
    CodeGenerator ifTrue, 
    CodeGenerator ifFalse
  ) {
    this.trueLabel = new Label();
    this.endLabel = new Label();
    this.skipTrue = ByteCode.make_goto(endLabel);

    this.ifInstruction = ByteCode.make_if(ifInstruction, trueLabel);
    this.ifTrue = ifTrue;
    this.ifFalse = ifFalse;
    
    // lazyness - avoid work later
    instructions = new InstructionVector();
    instructions.add(this.ifInstruction);
    instructions.add(this.ifFalse);
    instructions.add(this.skipTrue);
    instructions.add(ByteCode.make_markLabel(trueLabel));
    instructions.add(this.ifTrue);
    instructions.add(ByteCode.make_markLabel(endLabel));
  }
  
  public Instruction getIfInstruction() {
    return ifInstruction;
  }
  
  public CodeGenerator getIfTrue() {
    return ifTrue;
  }
  
  public CodeGenerator getIfFalse() {
    return ifFalse;
  }
  
  public void writeCode(CodeContext ctx) throws CodeException {
    instructions.writeCode(ctx);
  }
  
  public int stackDepth() {
    // custom handlers needed because of jumps
    return
      ifInstruction.stackDepth() +
      Math.max(ifFalse.stackDepth(), ifTrue.stackDepth());
  }
  
  public int stackDelta() {
    // custom handler needed because of jumps
    return
      ifInstruction.stackDepth() +
      Math.max(ifFalse.stackDepth(), ifTrue.stackDepth()); // these should agree
  }
}
