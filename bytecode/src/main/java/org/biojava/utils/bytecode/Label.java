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
 * A Label used to mark a position in byte code.
 *
 * <p>
 * Labels are used as the targets for jumps, and for exception handlers. Labels
 * can be named. They implement CodeGenerator, which allows them to be added
 * to things like an InstructionVector. The writeCode method takes care of
 * marking the label with the context.
 * </p>
 *
 * @author Thomas Down
 * @author Matthew Pocock
 */
public class Label implements CodeGenerator {
  public String name;
  
  public Label() {
    name = null;
  }
  
  public Label(String name) {
    this.name = name;
  }
  
  public String toString() {
    if (name != null)
      return name;
    return super.toString();
  }
  
  public void writeCode(CodeContext ctx) throws CodeException {
    ctx.markLabel(this);
  }
  
  public int stackDepth() {
    return 0;
  }
  
  public int stackDelta() {
    return 0;
  }
}
