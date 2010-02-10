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

import java.beans.BeanInfo;
import java.beans.IntrospectionException;
import java.beans.Introspector;
import java.beans.PropertyDescriptor;
import java.lang.reflect.InvocationTargetException;
import java.util.AbstractMap;
import java.util.AbstractSet;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

/**
 * @author Matthew Pocock
 */
public class BeanAsMap
  extends
    AbstractMap
{
  // fixme: does this class work at all?
  
  private static Map beanInfoCache;
  private static Object[] NO_PARAMS;
  
  static {
    beanInfoCache = new HashMap();
    NO_PARAMS = new Object[] {};
  }
  
  private static BeanInfo getBeanInfo(Class clazz)
  throws IntrospectionException {
    BeanInfo bi = (BeanInfo) beanInfoCache.get(clazz);
    
    if(bi == null) {
      beanInfoCache.put(clazz, bi = Introspector.getBeanInfo(clazz));
    }
    
    return bi;
  }
  
  private final BeanInfo beanInfo;
  private final Object bean;
  private final PropertyDescriptor[] descriptors;
  private final Set entrySet;
  
  public BeanAsMap(Object bean)
  throws IntrospectionException {
    this.beanInfo = getBeanInfo(bean.getClass());
    this.bean = bean;
    this.descriptors = this.beanInfo.getPropertyDescriptors();
    this.entrySet = new PropertySet();
  }
  
  public int size() {
    return entrySet.size();
  }
  
  public Set entrySet() {
    return entrySet;
  }
  
  public Object put(Object key, Object value) {
    for(Iterator i = entrySet.iterator(); i.hasNext(); ) {
      PropertyEntry pe = (PropertyEntry) i.next();
      if(pe.getKey().equals(key)) {
        return pe.setValue(value);
      }
    }
    throw new IllegalArgumentException("BeanAsMap does not support key: " + key);
  }
  
  private class PropertySet
    extends
      AbstractSet
  {
    private PropertyEntry[] entries;
    
    public PropertySet() {
      entries = new PropertyEntry[descriptors.length];
      for(int i = 0; i < entries.length; i++) {
        entries[i] = new PropertyEntry(descriptors[i]);
      }
    }
    
    public int size() {
      return entries.length;
    }
    
    public Iterator iterator() {
      return new Iterator() {
        int i = 0;
        
        public boolean hasNext() {
          return i < entries.length;
        }
        
        public Object next() {
          return entries[i++];
        }
        
        public void remove() {
          throw new UnsupportedOperationException();
        }
      };
    }
  }
  
  private class PropertyEntry implements Map.Entry {
    private final PropertyDescriptor pd;
    
    public PropertyEntry(PropertyDescriptor pd) {
      this.pd = pd;
    }
    
    public Object getKey() {
      return pd.getName();
    }
    
    public Object getValue() {
      try {
        return pd.getReadMethod().invoke(bean, NO_PARAMS);
      } catch (IllegalAccessException iae) {
        throw new AssertionFailure("Could not set property", iae);
      } catch (InvocationTargetException ite) {
        throw new AssertionFailure("Could not invoke property", ite);
      }
    }
    
    public Object setValue(Object value) {
      try {
        Object old = getValue();
        pd.getWriteMethod().invoke(bean, new Object[] { value });
        
        return old;
      } catch (IllegalAccessException iae) {
        throw new AssertionFailure("Could not access property", iae);
      } catch (InvocationTargetException ite) {
        throw new AssertionFailure("Could not invoke property", ite);
      }
    }
  }
}
