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
 * Interface for an object which can produce Java bytecode.
 *
 * @author Thomas Down
 * @author Matthew Pocock
 */

public interface CodeGenerator {
    /**
     * Write the byte or bytes for this CodeGenerator to a CodeContext.
     *
     * @param ctx  a CodeContext to write to
     * @throws CodeException if there was some failure in writing to the context
     */
    public void writeCode(CodeContext ctx) throws CodeException;
    
    /**
     * Return the total depth of the stack required by this CodeGenerator.
     *
     * <p>For single byte-code instructions, this will be the same as
     * stackDelta() if stackDelta() is positive, zero otherwise. For a
     * compound instruction, this will be the maximum stack depth required to
     * execute all sub-instructions.</p>
     *
     * @return the stack depth needed
     */
    public int stackDepth();
    
    /**
     * Return the change in the stack dept this generator will cause.
     *
     * <p>In the case of an instruction that adds items to the stack, stackDelta
     * will be positive. For instructions that removes items from the stack,
     * this will be negative.</p>
     *
     * @return the change between stack depth before and after execution of this
     *   code
     */
    public int stackDelta();
}

