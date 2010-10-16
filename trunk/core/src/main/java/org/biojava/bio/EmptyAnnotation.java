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

import java.io.NotSerializableException;
import java.io.ObjectStreamException;
import java.io.Serializable;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Set;

import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.StaticMemberPlaceHolder;
import org.biojava.utils.Unchangeable;

/**
 * An always-empty Annotation.
 *
 * @author Matthew Pocock
 * @author Thomas Down
 * 
 * @since 1.0 as part of Annotation
 * @since 1.4 as top-level class
 * @see org.biojavax.EmptyRichAnnotation
 */
class EmptyAnnotation
extends Unchangeable
implements Annotation, Serializable {
  public Object getProperty(Object key) throws NoSuchElementException {
    throw new NoSuchElementException(
      "There are no keys in the Empty Annotation object: " +
      key
    );
  }
  
  public void setProperty(Object key, Object value)
  throws ChangeVetoException {
    throw new ChangeVetoException(
      "You can not add properties to the Empty Annotation object: " +
      key + " -> " + value
    );
  }
  
  public void removeProperty(Object key)
  throws ChangeVetoException 
  {
    throw new ChangeVetoException(
      "You cannot remove properties from the empty annotation (!)"
    );
  }
  
  public boolean containsProperty(Object key) {
    return false;
  }
  
  public Set keys() {
    return Collections.EMPTY_SET;
  }
  
  public Map asMap() {
    //return Collections.EMPTY_MAP; 1.3
    return new HashMap();
  }
  
  private Object writeReplace() throws ObjectStreamException {
    try {
      return new StaticMemberPlaceHolder(Annotation.class.getField("EMPTY_ANNOTATION"));
    } catch (NoSuchFieldException nsfe) {
      throw new NotSerializableException(nsfe.getMessage());
    }
  }
  
  public int hashCode() {
    return asMap().hashCode();
  }
  
  public boolean equals(Object o) {
    if (! (o instanceof Annotation)) {
      return false;
    }
    
    return ((Annotation) o).asMap().equals(asMap());
  }
}

