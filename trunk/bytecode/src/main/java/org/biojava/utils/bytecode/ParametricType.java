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
 * A template type.
 *
 * <p>Template types are resolved at code-generation type rather than at
 * Instruction generation type. They let you bind the concrete type for opcodes
 * at the last minute, so the same max conditional could be used for all
 * primative types, with the type only being bound at the last moment.</p>
 *
 * <p>Two ParametricType instances are the same if they are the same object,
 * regardless of their names.</p>
 *
 * @author Matthew Pocock
 */

public class ParametricType {
  private static CodeClass[] OBJECT_CC;
  
  static {
    OBJECT_CC = new CodeClass[] { IntrospectedCodeClass.forClass(Object.class) };
  }
  
  /**
   * Create a new ParametricType that claims nothing.
   *
   * @param name  the name given to this type
   * @return a new ParametricType instance with that name
   */
  public static ParametricType createType(String name) {
    return new ParametricType(name, false, false, false);
  }
  
  /**
   * Create a new ParametricType that claims to resolve to a primative type.
   *
   * @param name  the name given to this type
   * @return a new ParametricType instance with that name
   */
  public static ParametricType createPrimitiveType(String name) {
    return new ParametricType(name, true, false, false);
  }
  
  /**
   * Create a new ParametricType that claims to resolve to an object type.
   *
   * @param name  the name given to this type
   * @return a new ParametricType instance with that name
   */
  public static ParametricType createObjectType(String name) {
    return new ParametricType(name, false, true, false);
  }
  
  /**
   * Create a new ParametricType that claims to resolve to an array type. All
   * array types are object types.
   *
   * @param name  the name given to this type
   * @return a new ParametricType instance with that name
   */
  public static ParametricType createArrayType(String name) {
    return new ParametricType(name, false, true, true);
  }
  
  /**
   * Create a new ParametricType that claims to be castable to all the classes
   * in a list. Since neither Java nor bytecode support multiple inheritance,
   * the classes must either be interfaces, or classes that fall into an
   * inheritance path.
   *
   * @param name  the name given to this type
   * @param classes an array of Class objects that any bound type must be
   *   castable to
   * @return a new ParametricType that can bind to classes with these properties
   */
  public static ParametricType createType(
    String name,
    CodeClass[] classes
  ) {
    return new ParametricType(name, classes);
  }
  
  private final String name;
  private final boolean isPrimitive;
  private final boolean isObject;
  private final boolean isArray;
  private final CodeClass[] classes;
  
  private ParametricType(
    String name,
    boolean isPrimitive,
    boolean isObject,
    boolean isArray
  ) {
    this.name = name;
    this.isPrimitive = isPrimitive;
    this.isObject = isObject;
    this.isArray = isArray;
    if(isObject) {
      this.classes = OBJECT_CC;
    } else {
      this.classes = CodeUtils.EMPTY_LIST;
    }
  }
  
  private ParametricType(
    String name,
    CodeClass[] classes
  ) {
    this.name = name;
    this.classes = classes;
    this.isObject = true;
    this.isPrimitive = false;
    this.isArray = false;
  }
  
  /**
   * Get the name of this type.
   *
   * Names are not unique.
   *
   * @return the name given to this type
   */
  public String getName() {
    return name;
  }
  
  /**
   * Discover if this type must resolve to a primative.
   *
   * <p>It is an error for a parametric type to resolve to a non-primative if
   * this flag is set.</p>
   *
   * @return true if this is guaranteed to resolve to a primative
   */
  public boolean isPrimitive() {
    return isPrimitive;
  }
  
  public boolean isObject() {
    return isObject;
  }
  
  public boolean isArray() {
    return isArray;
  }
  
  public boolean canAccept(CodeClass cc) {
    if(cc.isArray() && this.isArray()) {
      return true;
    }
    
    if(!cc.isPrimitive() && this.isObject()) {
      return true;
    }
    
    if(cc.isPrimitive() && this.isPrimitive()) {
      return true;
    }
    
    return false;
  }
  
  public CodeClass[] getClasses() {
    return classes;
  }
  
  public String toString() {
    return "GenericType:" + name;
  }
}
