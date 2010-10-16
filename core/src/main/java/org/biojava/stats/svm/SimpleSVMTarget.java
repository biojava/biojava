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


package org.biojava.stats.svm;

import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

/**
 * No-frills implementation of SVMTarget.
 *
 * @author Matthew Pocock
 */
public class SimpleSVMTarget implements SVMTarget {
  private Set itemTargetSet;
  private Map itemToItemTarget;

  public SimpleSVMTarget() {
    itemTargetSet = new HashSet();
    itemToItemTarget = new HashMap();
  }
  
  public SimpleSVMTarget(Collection items) {
    this();
    for(Iterator i = items.iterator(); i.hasNext(); ) {
      addItem(i.next());
    }
  }
  
  public Set items() {
    return itemToItemTarget.keySet();
  }
  
  public Set itemTargets() {
    return itemTargetSet;
  }

  public double getTarget(Object item) {
    return ((ItemValue) itemToItemTarget.get(item)).getValue();
  }
  
  public void setTarget(Object item, double target) {
    ItemValue iv = (ItemValue) itemToItemTarget.get(item);
    iv.setValue(target);
  }
  
  public void addItem(Object item) {
    ItemValue iv = new SimpleItemValue(item, 0.0);
    itemToItemTarget.put(item, iv);
    itemTargetSet.add(iv);
  }
  
  public void addItemTarget(Object item, double target) {
    ItemValue iv = new SimpleItemValue(item, target);
    itemToItemTarget.put(item, iv);
    itemTargetSet.add(iv);
  }
  
  public void removeItem(Object item) {
    itemToItemTarget.remove(item);
    itemTargetSet.remove(item);
  }
  
  public void clear() {
    itemToItemTarget.clear();
    itemTargetSet.clear();
  }
}
