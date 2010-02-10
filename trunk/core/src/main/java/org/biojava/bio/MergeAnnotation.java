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

import java.io.Serializable;
import java.util.AbstractMap;
import java.util.AbstractSet;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Set;

import org.biojava.utils.AbstractChangeable;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeForwarder;
import org.biojava.utils.ChangeListener;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;

/**
 * Merged view onto a list of underlying Annotation objects.
 * Currently immutable (but reflects changes to underlying objects). Annotations
 * near the beginning of the list will have properties that take
 * precedence. It is possible to get the ordering of the annotations, or to
 * change it by removing and re-adding methods.
 * This Annotation implementation is immutable.
 *
 * @author Thomas Down
 * @author Matthew Pocock
 * @author Greg Cox
 * @author Francois Pepin
 * @since 1.2
 *
 * Use these when you have a list of Annotation instances that
 * need to be viewed as one. For example, if you have annotation for a feature
 * from a local database, in-memory objects and a web-page, you could build
 * three Annotation instances and merge them using a MergeAnnotation.
 */

public class MergeAnnotation
        extends
        AbstractChangeable
        implements
        Annotation,
        Serializable {
  private transient ChangeListener propertyForwarder = null;

  private List mergeSet;

  {
    mergeSet = new ArrayList();
  }

    /**
   * ChangeType of ChangeEvent fired before and after an annotation is added
   * to MergeAnnotation.
   *
   */
  public static final ChangeType ANNOTATION_CHANGED = new ChangeType(
    "annotation added",
    "org.biojava.bio.MergeAnnotation",
    "ANNOTATION_CHANGED"
  );
  
  /**
   * ChangeType of ChangeEvent fired before and after an annotation is added
   * to MergeAnnotation.
   *
   */
  public static final ChangeType ANNOTATION_ADD = new ChangeType(
    "annotation added from List",
    "org.biojava.bio.MergeAnnotation",
    "ANNOTATION_ADD",
    ANNOTATION_CHANGED
  );

    /**
   * ChangeType of ChangeEvent fired before and after an annotation is added
   * to MergeAnnotation.
   *
   */
  public static final ChangeType ANNOTATION_REMOVE = new ChangeType(
    "annotation deleted from List",
    "org.biojava.bio.MergeAnnotation",
    "ANNOTATION_REMOVE",
    ANNOTATION_CHANGED
  );


  
  /**
   * Add a new Annotation to to the end of the list to be merged.
   *
   * Use this to alter the Annotations being merged
   *
   * @param ann  the Annotation to add
   * @throws ChangeVetoException if the annotation could not be added
   */
  public void addAnnotation(Annotation ann)
          throws ChangeVetoException {
     if(!hasListeners())
       mergeSet.add(ann);
     else{
       ChangeEvent ce = new ChangeEvent(this,MergeAnnotation.ANNOTATION_ADD,ann);
       ChangeSupport changeSupport = super.getChangeSupport(MergeAnnotation.ANNOTATION_ADD);
       synchronized(changeSupport) {
        changeSupport.firePreChangeEvent(ce);
        mergeSet.add(ann);
        changeSupport.firePostChangeEvent(ce);
      }
     }
  }

  /**
   * Gets an unmodifiable view of the list of Annotations that are part of the
   * MergeAnnotation. Lower indices Annotation have precedence if 2
   * Annotations share the same property.
   * 
   * @return an unmodifiable <code>List</code> of the Annotations that form
   * this MergeAnnotation.
   */
  public List getAnnotations()
  {
    return Collections.unmodifiableList(mergeSet);
  }

  /**
   * Remove an Annotation from the list. This can be used to change the
   * ordering of the Annotations by re-adding it later.
   *
   * @param ann an <code>Annotation</code> to be removed.
   * @exception ChangeVetoException if an error occurs
   */
  public void removeAnnotation(Annotation ann)
    throws ChangeVetoException {
    if(!hasListeners())
       mergeSet.remove(ann);
     else{
       ChangeEvent ce = new ChangeEvent(this,MergeAnnotation.ANNOTATION_REMOVE,ann);
       ChangeSupport changeSupport = super.getChangeSupport(MergeAnnotation.ANNOTATION_REMOVE);
       synchronized(changeSupport) {
         changeSupport.firePreChangeEvent(ce);
         mergeSet.remove(ann);
         changeSupport.firePostChangeEvent(ce);
       }
     }
  }
  
  
  protected ChangeSupport getChangeSupport(ChangeType changeType) {
    ChangeSupport changeSupport = super.getChangeSupport(changeType);

    if (
            (Annotation.PROPERTY.isMatchingType(changeType) || changeType.isMatchingType(Annotation.PROPERTY))
            &&
            propertyForwarder == null
    ) {
      propertyForwarder = new PropertyForwarder(
              MergeAnnotation.this,
              changeSupport
      );
      for (Iterator i = mergeSet.iterator(); i.hasNext();) {
        Annotation a = (Annotation) i.next();

        a.addChangeListener(propertyForwarder, Annotation.PROPERTY);
      }
    }

    return changeSupport;
  }

  public void setProperty(Object key, Object value) throws ChangeVetoException {
    throw new ChangeVetoException("MergeAnnotations don't allow property setting at the moment");
  }

  public void removeProperty(Object key) throws ChangeVetoException {
    throw new ChangeVetoException("MergeAnnotations don't allow property removal at the moment");
  }

  public Object getProperty(Object key) {
    for (Iterator i = mergeSet.iterator(); i.hasNext();) {
      Annotation a = (Annotation) i.next();
      if (a.containsProperty(key)) {
        return a.getProperty(key);
      }
    }
    throw new NoSuchElementException("Can't find property " + key);
  }

  public boolean containsProperty(Object key) {
    for (Iterator i = mergeSet.iterator(); i.hasNext();) {
      Annotation a = (Annotation) i.next();
      if (a.containsProperty(key)) {
        return true;
      }
    }

    return false;
  }

  public Set keys() {
    Set s = new HashSet();
    for (Iterator i = mergeSet.iterator(); i.hasNext();) {
      Annotation a = (Annotation) i.next();
      s.addAll(a.keys());
    }
    return s;
  }

  public Map asMap() {
    return new MAMap();
  }

  private class MAEntrySet extends AbstractSet {
    private MAEntrySet() {
      super();
    }

    public Iterator iterator() {
      return new Iterator() {
        Iterator ksi = MergeAnnotation.this.keys().iterator();

        public boolean hasNext() {
          return ksi.hasNext();
        }

        public Object next() {
          Object k = ksi.next();
          Object v = getProperty(k);
          return new MAMapEntry(k, v);
        }

        public void remove() {
          throw new UnsupportedOperationException();
        }
      };
    }

    public int size() {
      return MergeAnnotation.this.keys().size();
    }
  }

  private class MAMapEntry implements Map.Entry {
    private Object key;
    private Object value;

    private MAMapEntry(Object key, Object value) {
      this.key = key;
      this.value = value;
    }

    public Object getKey() {
      return key;
    }

    public Object getValue() {
      return value;
    }

    public Object setValue(Object v) {
      throw new UnsupportedOperationException();
    }

    public boolean equals(Object o) {
      if (!(o instanceof Map.Entry)) {
        return false;
      }

      Map.Entry mo = (Map.Entry) o;
      return ((key == null ? mo.getKey() == null : key.equals(mo.getKey())) &&
              (value == null ? mo.getValue() == null : value.equals(mo.getValue())));
    }

    public int hashCode() {
      return (key == null ? 0 : key.hashCode()) ^ (value == null ? 0 : value.hashCode());
    }
  }

  private class MAMap extends AbstractMap {
    MAEntrySet es;

    private MAMap() {
      super();
      es = new MAEntrySet();
    }

    public Set entrySet() {
      return es;
    }

    public Set keySet() {
      return MergeAnnotation.this.keys();
    }

    public Object get(Object key) {
      try {
        return getProperty(key);
      } catch (NoSuchElementException ex) {
      }

      return null;
    }
  }

  /**
   * Listener used to forward changes for any of the underlying annotations to
   * listeners on this annotation.
   *
   * @author Thomas Down
   * @author Matthew Pocock
   * @since 1.2
   */
  protected class PropertyForwarder extends ChangeForwarder {
    /**
     * Create a new forwarder on behalf of a source using the change support.
     * @param source  the new source of events
     * @param cs      the ChangeSupport used to manage listeners
     */
    public PropertyForwarder(Object source, ChangeSupport cs) {
      super(source, cs);
    }

    public ChangeEvent generateEvent(ChangeEvent ce) {
      ChangeType ct = ce.getType();
      if (ct == Annotation.PROPERTY) {
        Object curVal = ce.getChange();
        if (curVal instanceof Object[]) {
          Object[] cur = (Object[]) curVal;
          if (cur.length == 2) {
            Object key = cur[0];
            Object value = cur[0];
            if (getProperty(key) != value) {
              return new ChangeEvent(
                      getSource(),
                      Annotation.PROPERTY,
                      curVal,
                      ce.getPrevious(),
                      ce
              );
            }
          }
        }
      }
      return null;
    }
  }
}
