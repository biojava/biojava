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

package org.biojavax.ga.util;

import java.util.AbstractSet;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.NoSuchElementException;


/**
 * <p>Inspred by the BioJava Distribution objects the WeightedSet is a map from
 * a Key to a Weight. Unlike Distributions the Keys do not have to be Symbols.
 * In the GA package the WeightedMap is useful for sampling Organisms according
 * to their fitness.</p>
 *
 * <p>When Symbols are added or their weights are set then the weights are internally
 * normalized to 1</p>
 *
 * @author Mark Schreiber
 * @version 1.0
 * @since 1.5
 */

public class WeightedSet extends AbstractSet implements java.io.Serializable{
  private HashMap key2Weight;
  double totalWeight;

  public WeightedSet() {
    key2Weight = new HashMap();
  }

  /**
   * Converts the Set to a map from key <code>Objects</code> to <code>Double</code>
   * weights.
   * @return a Map with all the key-weight mappings. Weights are not normalized in this map.
   */
  public Map asMap(){
    return key2Weight;
  }

  /**
   * Randomly samples an <code>Object</code> from the <code>Set</code> according
   * to its weight.
   * @return the Object sampled.
   */
  public Object sample(){
    double p = Math.random();
    for (Iterator i = this.iterator(); i.hasNext(); ) {
      Object o = i.next();
      double weight = getWeight(o);

      p -= weight;
      if(p <= 0.0){
        return o;
      }
    }
    throw new org.biojava.bio.BioError("Cannot sample an object, does this set contain any objects?");
  }

  /**
   * Determines the normalized weight for <code>o</code>
   * @param o the <code>Object</code> you want to know the weight of
   * @return the normalized weight
   * @throws NoSuchElementException if <code>o</code> is not found in this set
   */
  public double getWeight(Object o) throws NoSuchElementException{
    if(!( key2Weight.containsKey(o)))
      throw new NoSuchElementException(o+" not found in this WeightedSet");

    Double d = (Double)key2Weight.get(o);
    if(totalWeight == 0.0)
      return 0.0;


    return d.doubleValue() / totalWeight;
  }

  /**
   * The total weight that has been added to this Set.
   * @return the total weight (the value that can be used for normalizing)
   */
  protected double getTotalWeight(){
    return totalWeight;
  }

  /**
   * Sets the weight of an <code>Object</code>. If the <code>Object</code> is
   * not in this <code>Set</code> then it is added.
   * @param o the <code>Object</code>
   * @param w the weight.
   * @throws IllegalArgumentException if <code>w</code> is < 0.0
   */
  public void setWeight(Object o, double w){
    if(w < 0.0){
      throw new IllegalArgumentException("Weight must be >= 0.0");
    }
    if(key2Weight.containsKey(o)){
      remove(o);
    }
    totalWeight += w;
    key2Weight.put(o, Double.valueOf(w));
  }

  public boolean contains(Object o) {
    return key2Weight.containsKey(o);
  }

  public boolean remove(Object o) {
    if(key2Weight.containsKey(o)){
      totalWeight -= ((Double)key2Weight.get(o)).doubleValue();
      key2Weight.remove(o);
      return true;
    }
    return false;
  }

  public boolean isEmpty() {
    return key2Weight.isEmpty();
  }
  public boolean retainAll(Collection c) {
    boolean b = false;
    Collection toRemove = new ArrayList();

    for (Iterator i = iterator(); i.hasNext(); ) {
      Object item = i.next();
      if(c.contains(item) == false){
        b = true;
        toRemove.add(item);
      }
    }

    removeAll(toRemove);

    return b;
  }

  /**
   * Adds a new <code>Object</code> with a weight of zero. Equivalent to
   * setWeight(o, 0.0);
   * @param o the object to add.
   * @return true if this Object has not been added before.
   */
  public boolean add(Object o) {
    boolean b = !(key2Weight.containsKey(o));
    setWeight(o, 0.0);
    return b;
  }
  public int size() {
    return key2Weight.size();
  }

  public boolean containsAll(Collection c) {
    if(size() == 0)
      return false;

    for (Iterator i = iterator(); i.hasNext(); ) {
      Object item = i.next();
      if(!(key2Weight.containsKey(item))){
        return false;
      }
    }
    return true;
  }
  public Object[] toArray() {
    Object[] o = new Object[size()];
    int j = 0;
    for (Iterator i = iterator(); i.hasNext(); ) {
      o[j++] = i.next();
    }

    return o;
  }

  public void clear() {
    key2Weight = new HashMap();
    totalWeight = 0.0;
  }

  /**
   * Returns an unmodifiable iterator over the keys of the set.
   * @return an Iterator
   */
  public Iterator iterator() {
    return Collections.unmodifiableSet(key2Weight.keySet()).iterator();
  }

  public boolean addAll(Collection c) {
    boolean b = false;

    for (Iterator i = c.iterator(); i.hasNext(); ) {

      Object item = i.next();
      if(!(key2Weight.containsKey(item)))
         b = true;

      add(item);
    }
    return b;
  }
}