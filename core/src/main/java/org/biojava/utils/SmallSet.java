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

import java.util.AbstractSet;
import java.util.Collection;
import java.util.ConcurrentModificationException;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

/**
 * Lightweight implementation of Set which uses little memory to store a small
 * number of items, at the expense of scalability. Not recomended for more than
 * 20-30 items.
 *
 * <p>
 * This implementation has the useful property that the iteration order is the
 * same as the order in which the items are added.
 * </p>
 *
 * @author Matthew Pocock
 */

public class SmallSet extends AbstractSet {
  private int vid = 0;
  private Object[] items;
  private int numItems;
  
  public SmallSet() {
    this(2);
  }
  
  public SmallSet(int size) {
    items = new Object[size];
    numItems = 0;
  }
  
  public SmallSet(Collection col) {
      if (! (col instanceof Set)) {
          col = new HashSet(col);
      }
      numItems = col.size();
      items = new Object[numItems];
      items = col.toArray(items);
  }
  
  public boolean contains(Object o) {
    for(int i = 0; i < numItems; i++) {
      if(items[i].equals(o)) {
        return true;
      }
    }
    
    return false;
  }
  
  public int size() {
    return numItems;
  }
  
  public boolean add(Object o) {
    if(this.contains(o)) {
      return false;
    }
    
    if(numItems == items.length) {
      Object[] tmp = new Object[items.length * items.length];
      System.arraycopy(items, 0, tmp, 0, items.length);
      this.items = tmp;
    }
    
    this.items[this.numItems++] = o;
    vid++;
    
    return true;
  }
  
  // remove by index
  private boolean remove(int i) {
    System.arraycopy(items, i, items, i-1, numItems - i);
    numItems--;
    vid++;
    
    return true;
  }
  
  public Iterator iterator() {
    return new Iterator() {
      int vid = SmallSet.this.vid;
      
      int i = 0;
      
      public boolean hasNext() {
        validate();
        return i < numItems;
      }
      
      public Object next() {
        validate();
        return items[i++];
      }
      
      public void remove() {
        validate();
        SmallSet.this.remove(i--);
        this.vid = SmallSet.this.vid;
      }
      
      private void validate() {
        if(this.vid != SmallSet.this.vid) {
          throw new ConcurrentModificationException();
        }
      }
    };
  }
}
