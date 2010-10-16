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
 * A local variable.
 *
 * <p>
 * Local variables are used as identifiers for things that can be stored and
 * loaded from the local variable slots associated with each method. By using
 * LocalVariable intances, you are removed from the taudry task of book-keeping
 * these slots.
 * </p>
 *
 * <p>
 * To use a local variable, create an intance and then use it in a code
 * generator. The method will keep track of which local variables are in scope,
 * and will handle all the nastiness for you. You can re-use the same
 * local variable instance in different contexts, and it will be sanely
 * allocated different or the same slots.
 * </p>
 *
 * <p>
 * The JVM stores some things in single words, and others in pairs of words.
 * To hide this detail from you, local variables take a class that indicates
 * the type of thing they will store. Please populate this sensibly.
 * </p>
 *
 * @author Thomas Down
 * @author Matthew Pocock
 */
public final class LocalVariable {
  private final String name;
  private final CodeClass clazz;

  /**
   * Create a new local variable that will store values of a given type.
   *
   * @param clazz  the type of the values stored in this variable
   */
  public LocalVariable(CodeClass clazz) {
    this.clazz = clazz;
    this.name = null;
  }

  /**
   * Create a new local variable with a type and a name.
   *
   * <p>
   * The name may appear in debug output from the generator, and possibly in
   * stack-traces.
   * </p>
   *
   * @param clazz  the type of the values stored in this variable
   * @param name   the name of the variable
   */
  public LocalVariable(CodeClass clazz, String name) {
    this.clazz = clazz;
    this.name = name;
  }

  public String getName() {
    return name;
  }

  public int needSlots() {
    return (clazz == CodeUtils.TYPE_LONG || clazz == CodeUtils.TYPE_DOUBLE) ? 2 : 1;
  }

  public CodeClass getType() {
    return clazz;
  }

  public String toString() {
    return
            super.toString() +
            "[slots: " + needSlots() + ", name: " + name +
            ", class: " + clazz + "]";
  }
}
