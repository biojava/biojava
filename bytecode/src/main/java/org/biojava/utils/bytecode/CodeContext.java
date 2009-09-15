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
 * Interface which encapsulates the stream to which Java bytecode can
 * be written.
 *
 * <p>
 * The context takes care of all the book-keeping tasks associated with emitting
 * well-formed byte code. For example, the context manages jumps and local
 * variables.
 * </p>
 *
 * <p>
 * Most of the funcionality here is very low level. You will almost certainly
 * want to use CodeGenerator instances to manipulate a CodeContext, rather than
 * writing to it yourself.
 * </p>
 *
 * @author Thomas Down
 * @author Matthew Pocock
 */

public interface CodeContext {
  /**
   * Get the class for which a method is being generated.
   *
   * @return the current class
   */
  public CodeClass getCodeClass();

  /**
   * Get the method which is being generated.
   *
   * @return the current method
   */
  public CodeMethod getCodeMethod();

  /**
   * Get the constants pool for this context.
   *
   * @return the contant pool
   */
  public ConstantPool getConstants();

  // Write methods.

  /**
   * Write a single byte to the context.
   *
   * <p>
   * This can be used both to write opcodes and to write byte data to the
   * context.
   * </p>
   *
   * @param  b the byte to write
   */
  public void writeByte(byte b) throws CodeException;

  /**
   * Write a short (2 bytes) to the context.
   *
   * @param i  the short to write
   */
  public void writeShort(int i) throws CodeException;

  /**
   * Write the offset of a Label to the context.
   *
   * <p>This can be called before or after markLabel is invoked for the
   * corresponding label. The context will ensure that the offset is
   * correctly written before the method is fully emitted.</p>
   *
   * @param lab  the Label to write
   */
  public void writeLabel(Label lab) throws CodeException;

  /**
   * Resolve a local variable to the local variable slot assigned to it.
   *
   * <p>The context will ensure that local variables are stored in their
   * own bits of the local variable area. It may chose to re-use portions
   * of this area as local variables go out of scope.</p>
   *
   * @param lv  the LocalVariable to resolve
   * @return the index of the local variable slot
   */
  public int resolveLocal(LocalVariable lv) throws CodeException;

  /**
   * Mark a label at the current point in the stream.
   *
   * <p>This can be used as the target for branching instructions, such as
   * GOTO and IF.</p>
   *
   * @param lab the Label to mark
   * @throws CodeException if the label has previously been marked
   */
  public void markLabel(Label lab) throws CodeException;

  /**
   * Register a concrete type for a parametric type.
   *
   * <p>This is the mechanism where-by real CodeClass types are associated
   * with the virtual ParametricType types. If type pubishes that it
   * is a primative, an object or an array, then the concreteType must be
   * compattible. It's an error to bind the VOID type.</p>
   *
   * @for.developer
   *  You should probably call
   * <code>ParametricType.canAccept(concreteType)</code> to make sure of this.
   * This implementation will shield you from any modifications to the exact
   * semantics of ParametricType.
   *
   * @param type  ParametricType the parametric type to register
   * @param concreteType  the CodeClass that it resolves to
   * @throws CodeException if the type has already been registered or if the
   *   guarantees about type made in the parametric type are violated
   */
  public void registerParametricType(ParametricType type, CodeClass concreteType)
          throws CodeException;

  /**
   * Resolve a parametric type to a concrete class.
   *
   * <p>The type will be resolved by first searching through all those
   * registered with this context. If it found there, this value is returned.
   * If it is not, then the emediate parent context is searched. This parent
   * is responsible for searching its parent and so on. If a context has no
   * parent and the type is not registered with this context, it should raise
   * a CodeException.</p>
   *
   * @param type  the ParametricType to resolve
   * @return the  ColdeClass associated with that parametric type
   * @throws CodeException if the type has not been registered
   */
  public CodeClass resolveParametricType(ParametricType type)
          throws CodeException;

  /**
   * Open a sub context.
   *
   * <p>The sub context should inherit all the state of the parent context.
   * Modifications to the state of the child (for example, local variable
   * management or registered labels) should not be propogated to the parent.
   * Modifications to the parent should not be propogated to the child.
   * However, contexts should be used serialy, and only one context should be
   * accessible to a code generator at a time, so in practice, there should be
   * no way for a code generator using a child context to alter the state of
   * or discover the state of the parent context.</p>
   */
  public CodeContext subContext();

  /**
   * Open the context for writing.
   *
   * <p>This must be called before any code writing methods are called. It
   * can not be called more than once.</p>
   */
  public void open() throws CodeException;

  /**
   * Close the context for writing. It is at this point that any process
   * necisary for comitting the bytecode will be executed.
   *
   * <p>This must be called after all code writing methods have been called.
   * It can not be called more than once.</p>
   */
  public void close() throws CodeException;

  // Exception tables

  /**
   * Add an exception table entry.
   *
   * @param startHandled    the beginning of the try block
   * @param endHandled      the end of the try block
   * @param eClass          the exception class
   * @param handler         the beginning of the exception handler
   * @throws CodeException
   */ 
  public void addExceptionTableEntry(Label startHandled,
                                     Label endHandled,
                                     CodeClass eClass,
                                     Label handler)
          throws CodeException;

}
