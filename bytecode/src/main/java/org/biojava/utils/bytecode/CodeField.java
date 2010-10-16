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
 * Wrap up details about a field in a Java class file.
 *
 * <p>
 * Instances of this type will be instantiated by CodeClass instances, using
 * the getField() methods.
 * </p>
 *
 * @author Thomas Down
 * @author Matthew Pocock
 */
public final class CodeField {
  private final String name;
  private final CodeClass clazz;
  private final int modifiers;
  private final CodeClass container;

  CodeField(CodeClass container, String name, CodeClass clazz, int mods) {
    this.container = container;
    this.name = name;
    this.clazz = clazz;
    this.modifiers = mods;
  }

  /**
   * Get the name of the field.
   *
   * @return  the name of the field
   */
  public String getName() {
    return name;
  }

  /**
   * Get the fully qualified name of the field.
   *
   * @return the fully qualified name
   */
  public String getFullName() {
    return container.getName() + "." + getName();
  }

  /**
   * Get the class that contains this field.
   *
   * @return the containing class
   */
  public CodeClass getContainingClass() {
    return container;
  }

  /**
   * Get the type of the field.
   *
   * @return
   */
  public CodeClass getType() {
    return clazz;
  }

  /**
   * Get the moddifiers applied to this field.
   *
   * @return the modifiers
   */
  public int getModifiers() {
    return modifiers;
  }

  public String toString() {
    return super.toString() +
            " type: " + getType() +
            " class: " + clazz.getName() +
            " name: " + getName();
  }
}
