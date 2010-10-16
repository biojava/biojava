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
 * Wrap up details about a method in a Java class file
 *
 * @author Thomas Down
 * @author Matthew Pocock
 */

public interface CodeMethod {
  /**
   * The name of the method.
   *
   * @return the method name
   */
  public String getName();

  /**
   * The class that contains this method
   *
   * @return the containing class
   */
  public CodeClass getContainingClass();

  /**
   * The fully qualified name for this class
   *
   * @return the full name
   */
  public String getFullName();

  /**
   * A human-readable description of the class
   *
   * @return the class description
   */
  public String getDescriptor();

  /**
   * Get the modifiers, such as PUBLIC, ABSTRACT and so on
   *
   * @return the class modifiers
   */
  public int getModifiers();

  /**
   * Get the return type
   *
   * @return the return type
   */
  public CodeClass getReturnType();

  /**
   * Get the number of parameters taken by this method
   *
   * @return the number of parameters
   */
  public int numParameters();

  /**
   * Get the type of the parameter at a given position
   *
   * @param pos  the position to fetch the parameter type for
   * @return the type of the parameter at that position
   */
  public CodeClass getParameterType(int pos);
}
