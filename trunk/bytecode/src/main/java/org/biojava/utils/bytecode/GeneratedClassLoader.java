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

import java.io.*;
import java.util.*;

/**
 * A class loader that actually produces real Java classes from
 * GeneratedCodeClass instances.
 *
 * @author Thomas Down
 * @author Matthew Pocock
 */
public class GeneratedClassLoader extends ClassLoader {
  private Set generatedClasses = new HashSet();

  /**
   * Create a new loader with the default parent.
   */
  public GeneratedClassLoader() {
    super();
  }

  /**
   * Create a new loader with an explicitly set parent class loader.
   *
   * @param parent  the parent ClassLoader
   */
  public GeneratedClassLoader(ClassLoader parent) {
    super(parent);
  }

  /**
   * Define a class based upon a GeneratedCodeClass.
   *
   * @param cc  the GeneratedCodeClass to define
   * @return    the newly defined class
   * @throws CodeException  if there was a failure defining the class
   */
  public Class defineClass(GeneratedCodeClass cc) throws CodeException {
    try {
	    ByteArrayOutputStream bos = new ByteArrayOutputStream();
	    cc.createCode(bos);
	    byte[] code = bos.toByteArray();
      Class clazz = defineClass(cc.getName(), code, 0, code.length);
      generatedClasses.add(clazz.getName());
      return clazz;
    } catch (IOException ex) {
	    // Seems unlikely...
	    throw new CodeException();
    }
  }

  /**
   * Discover if a class for this name has already been defined by this class
   * loader.
   *
   * @param name  the name of the class
   * @return      true if the class has already been defined by this loader
   */
  public boolean hasGeneratedClass(String name) {
    return generatedClasses.contains(name);
  }
}
