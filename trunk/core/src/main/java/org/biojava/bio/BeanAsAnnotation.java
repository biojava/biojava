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

import java.beans.IntrospectionException;
import java.util.Map;

import org.biojava.utils.BeanAsMap;

/**
 * Create an Annotation with properties matching those of a JavaBean instance.
 * <em>Note: this class is experimental and only partialy implemented.</em>
 *
 * @since 1.3
 * @author Matthew Pocock
 */

public class BeanAsAnnotation extends AbstractAnnotation {
  private Map properties;
  
  protected final Map getProperties() {
    return properties;
  }
  
  protected final boolean propertiesAllocated() {
    return true;
  }
  
  /**
   * Create a new BeanAsAnnotation for a bean.
   *
   * @param bean the JavaBean to view
   * @throws IntrospectionException  if the bean could not be introspected
   */
  public BeanAsAnnotation(Object bean)
  throws IntrospectionException {
    properties = new BeanAsMap(bean);
  }
}

