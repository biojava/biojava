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

package org.biojava.bio;

import java.util.LinkedHashMap;
import java.util.Map;

/**
 * <p>A no-frills implementation of Annotation that is just a wrapper
 * around a <code>LinkedHashMap</code>.</p>
 * 
 * <p>It will allow you to set any property, but will throw exceptions
 * if you try to retrieve a property that is not set. Because of the use of
 * a <code>LinkedHashMap</code> properties are iterated in the order they 
 * are entered.</p>
 * When you need to make a new Annotation instance and will want to populate
 * it with properties, use this class.
 * 
 * @author Matthew Pocock
 * @author Greg Cox
 * @author Mark Schreiber
 * @since 1.0
 * 
 * @see org.biojavax.SimpleRichAnnotation
 */
public class SimpleAnnotation extends AbstractAnnotation {
  /**
   * The properties map. This may be null if no property values have
   * yet been set.
   */
  private Map properties;

  /**
   * Get the property value pairs in this annotation
   * @return a <code>LinkedHashMap</code> of properties and values
   */
  protected final Map getProperties() {
    if(!propertiesAllocated()) {
      properties = new LinkedHashMap();
    }
    return properties;
  }

  /**
   * Test to see if any properties have been allocated
   * @return true if any have (even if they have since been removed),
   * false otherwise
   */
  protected final boolean propertiesAllocated() {
    return properties != null;
  }

  /**
   * Create a new, empty SimpleAnnotation instance.
   */
  public SimpleAnnotation() {
    super();
  }
  
  /**
   * Create a new SimpleAnnotation by copying the properties from another
   * one. The new Annotation instance will be independant of the original.
   * 
   * @param ann the Annotation to copy
   */
  public SimpleAnnotation(Annotation ann) {
    super(ann);
  }
  
  /**
   * Create a new SimpleAnnotation using the values in a Map.
   *
   * @param map  the Map to copy properties out of
   */
  public SimpleAnnotation(Map map) {
    super(map);
  }
}
