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
package org.biojava.utils;

import java.io.BufferedInputStream;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Properties;
import java.util.StringTokenizer;

/**
 * a sub-class of java.util.Properties that provides the same constructors, adds two convenient load methods to load
 * the properties from files and, most importantly, adds getPropertyAsXXX() methods to get a property as an object of
 * type XXX.
 *
 * @author <A href="mailto:Gerald.Loeffler@vienna.at">Gerald Loeffler</A> for the 
 *         <A href="http://www.imp.univie.ac.at">IMP</A>
 */
public class TypedProperties extends Properties {
  /**
   * the default string of delimiter characters used by getAsStringList()
   */
  public static final String DEFAULT_DELIMITERS = ",;\t";
  
  /**
   * Creates an empty property list with no default values.
   */
  public TypedProperties() {
    super();
  }
  
  /**
   * Creates an empty property list with the specified defaults.
   * @param defaults the defaults.
   */
  public TypedProperties(Properties defaults) {
    super(defaults);
  }

  /**
   * Reads a property list (key and element pairs) from the file with the given file name.
   * @param fileName the file name. Not null.
   * @exception FileNotFoundException if the file does not exist, is a directory rather than a regular file, or for some
   *                                  other reason cannot be opened for reading.
   * @exception IOException if an error occurred when reading from the input stream created from the file with the given
   *                        name.
   */
  public void load(String fileName) throws FileNotFoundException, IOException {
    if (fileName == null) {
	throw new IllegalArgumentException("fileName not null");
    }

    load(new BufferedInputStream(new FileInputStream(fileName)));
  }

  /**
   * Reads a property list (key and element pairs) from the given
   * file which is interpreted as a resource of the given class. The difference
   * between a normal file and a resource file is the way in which the file is
   * located: with a normal file, the filename is taken literally to load the
   * file from the file system, whereas with a resource file, the given name is
   * used to ask the class loader of the given class to load the file
   * (see java.lang.Class.getResourceAsStream() and
   * java.lang.ClassLoader.getSystemResourceAsStream()).
   *
   * @see java.lang.Class
   * @see java.lang.ClassLoader
   *
   * @param clazz the class with which the resource identified by resourceName
   *               is taken to be associated with
   *               (java.lang.Class.getResourceAsStream() on this Class
   *               object is used to load the resource). If clazz is null, the
   *               resource is considered to be a system resource, and
   *               java.lang.ClassLoader.getSystemResourceAsStream() is used to
   *               load the resource.
   * @param resourceName the name of the resource from which to load the properties. It is a precondition that the
   *                     resource with this name exists (regardless whether it is interpreted as a system resource or a
   *                     class resource), otherwise an IllegalArgumentException is thrown.
   * @exception IOException if an error occurred when reading from the input stream created from the given resource.
   */
  public void load(Class clazz, String resourceName) throws IOException {
    InputStream is = null;
    if (clazz == null) {
      // load resource as system resource
      is = ClassLoader.getSystemResourceAsStream(resourceName);
      if (is == null) {
	  throw new IllegalArgumentException("system reource " + resourceName + " must exist");
      }
    } else {
      // load resource as class resource
      is = clazz.getResourceAsStream(resourceName);
      if (is == null) {
	  throw new IllegalArgumentException("resource " + resourceName + " associated with class " + clazz + " must exist");
      }
    }
    load(is);
  }

  /**
   * Searches for the property with the specified key in this property list. If the key is not found in this property
   * list, the default property list, and its defaults, recursively, are then checked. The method returns null if the
   * property is not found. If the property is found its value is parsed as an integer and returned. If parsing the 
   * value fails, a NumberFormatException is thrown.
   *
   * @param key the property key.
   * @return the integer value of the property with the given key or null if the given key is not associated with a
   *         property.
   * @exception NumberFormatException if the property associated with the given key does not have an integer value.
   */
  public Integer getPropertyAsInteger(String key) throws NumberFormatException {
    String v = getProperty(key);
    if (v == null) return null;
    return new Integer(v);
  }

  /**
   * Searches for the property with the specified key in this property list. If the key is not found in this property
   * list, the default property list, and its defaults, recursively, are then checked. The method returns null if the
   * property is not found. If the property is found its value is parsed as a long and returned. If parsing the 
   * value fails, a NumberFormatException is thrown.
   *
   * @param key the property key.
   * @return the long value of the property with the given key or null if the given key is not associated with a
   *         property.
   * @exception NumberFormatException if the property associated with the given key does not have an integer value.
   */
  public Long getPropertyAsLong(String key) throws NumberFormatException {
    String v = getProperty(key);
    if (v == null) return null;
    return new Long(v);
  }

  /**
   * Searches for the property with the specified key in this property list. If the key is not found in this property
   * list, the default property list, and its defaults, recursively, are then checked. The method returns null if the
   * property is not found. If the property is found its value is parsed as a double and returned. If parsing the 
   * value fails, a NumberFormatException is thrown.
   *
   * @param key the property key.
   * @return the double value of the property with the given key or null if the given key is not associated with a
   *         property.
   * @exception NumberFormatException if the property associated with the given key does not have an integer value.
   */
  public Double getPropertyAsDouble(String key) throws NumberFormatException {
    String v = getProperty(key);
    if (v == null) return null;
    return new Double(v);
  }

  /**
   * Searches for the property with the specified key in this property list. If the key is not found in this property
   * list, the default property list, and its defaults, recursively, are then checked. The method returns null if the
   * property is not found. If the property is found its value is parsed as an boolean and returned. If parsing the 
   * value fails, a RuntimeException is thrown.
   * <p>
   * If the property value is equal, ignoring case, to the string "true" or "yes" then the boolean value
   * returned from this method is true. If the property value is equal, ignoring case, to the string
   * "false" or "no"  then the boolean value returned from this method is false.
   *
   * @param key the property key.
   * @return the boolean value of the property with the given key or null if the given key is not associated with a
   *         property.
   * @exception RuntimeException if the property associated with the given key does not have an integer value.
   */
  public Boolean getPropertyAsBoolean(String key) throws RuntimeException {
    String v = getProperty(key);
    if (v == null) return null;
    if      (v.equalsIgnoreCase("true")  || v.equalsIgnoreCase("yes")) return Boolean.TRUE;
    else if (v.equalsIgnoreCase("false") || v.equalsIgnoreCase("no"))  return Boolean.FALSE;
    else throw new RuntimeException("property value " + v + " is not parseable as a boolean");
  }

  /**
   * Searches for the property with the specified key in this property list. If the key is not found in this property
   * list, the default property list, and its defaults, recursively, are then checked. The method returns null if the
   * property is not found. If the property is found its value is parsed as a list of strings and returned as a List
   * object that contains only String objects. Parsing the property value as a list of strings can not fail and so this
   * method does not throw an exception.
   * <p>
   * The property value is interpreted as String objects (tokens) separated by one or more (consecutive) separator
   * characters taken from the delims string. Any of these characters separates the tokens and can hence not be part of
   * any token! The tokens identified in this way are put into a List in the order in which they appear in the property 
   * value. White space at the beginning and end of each token are removed before storing the token as an element of the
   * list (this includes white space at the beginning and end of the complete property value)! Empty strings are also
   * never added to the list, i.e. if after removal of white space from a token a token is the empty string, it is not
   * stored in the list! All this results in a very natural conversion of the property value into a list of strings:
   * only "real" (non-white-space, non-white-space-bounded, non-delimiter-containing) sub-strings from the
   * property value are put as string elements into the list.
   *
   * @param key the property key.
   * @param delims the string of allowed delimiter characters (not null and not empty).
   * @return the List of strings for the property with the given key or null if the given key is not associated with a
   *         property. An empty list is returned if a property with the given key exists but its value is empty or
   *         consists only of white space.
   */
  public List getPropertyAsStringList(String key, String delims) {
    if (delims == null || delims.length() == 0) {
	throw new IllegalArgumentException("delims != null && delims.length() > 0");
    }

    String v = getProperty(key);
    if (v == null) return null;

    List l = new ArrayList();

    v = v.trim();
    if (v.length() > 0) {
      StringTokenizer st = new StringTokenizer(v, delims, false);
      while (st.hasMoreTokens()) {
        String token = (st.nextToken()).trim();
        if (token != null && token.length() > 0) l.add(token); // store only non-empty tokens
      }
    }

    return l;
  }
  
  /**
   * just like getPropertyAsStringList(String key, String delims) but uses ',' (comma), ';' (semicolon) and '\t' (tab)
   * as the possible delimiters.
   */
  public List getPropertyAsStringList(String key) {
    return getPropertyAsStringList(key, DEFAULT_DELIMITERS);
  }

  public String toString() {
    return "TypedProperties";
  }
  
  public boolean equals(Object o) {
    if (o == this) return true;
    
    // this class extends another class than Object:
    if (!super.equals(o)) return false;
    
    // this and that are identical if we made it 'til here
    return true;
  }
  
  public int hashCode() {
    // this class extends another class than Object:
    int hc = super.hashCode();

    return hc;
  }
  
  public Object clone() {
    TypedProperties o = (TypedProperties) super.clone();

    return o;
  }
}
