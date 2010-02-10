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
 * Interface for Java classes within the bytecode generation framework.
 * Any class (or interface) can be viewed as a CodeClass, whether it
 * is pre-existing or being generated.
 *
 * @author Thomas Down
 * @author Matthew Pocock
 */

public interface CodeClass {
  String getName();

  String getJName();

  String getDescriptor();

  CodeClass getSuperClass();

  List getInterfaces();

  /**
   * Get all methods declared by this class and its super classes, removing
   * all super class methods that are over ridden.
   *
   * <p>
   * This should return methods, regardless of their accessability.
   * </p>
   *
   * @return a Set containing all methods
   */
  Set getMethods();

  /**
   * Get the name of all methods that could be invoked through this class with
   * a given name.
   *
   * @param name  the name of the method
   * @return a Set of CodeMethod instances with that name
   */
  Set getMethodsByName(String name);

  /**
   * Get a method by name and argument list.
   *
   * @param name  the name of the method
   * @param args  the arguments it takes
   * @return      a matching method
   * @throws NoSuchMethodException  if there is no maching method
   */
  CodeMethod getMethod(String name, CodeClass[] args)
          throws NoSuchMethodException;

  /**
   * Get a constructor by argument list.
   *
   * @param args  the arguments it takes
   * @return      a matching constructor
   * @throws NoSuchMethodException  if there is no matching constructor
   */
  CodeMethod getConstructor(CodeClass[] args)
          throws NoSuchMethodException;

  /**
   * Get a field by its name.
   *
   * @param name  the field name
   * @return      a CodeField representing the field
   * @throws NoSuchFieldException if there is no field by that name accessible
   *   through this class
   */
  CodeField getFieldByName(String name)
          throws NoSuchFieldException;

  /**
   * Get all fields accessible through this class.
   *
   * @return  a Set of all accessible fields
   */
  Set getFields();

  /**
   * Get the modifiers associated with the class.
   *
   * @return  the modifier integer
   */
  int getModifiers();

  /**
   * Discover if the class represents a primitive type.
   *
   * @return  true if the class represents a primative type
   */
  public boolean isPrimitive();

  /**
   * Discover if the class is an array type.
   *
   * @return  true if the class is an array type
   */
  public boolean isArray();
}
