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
 * A no-frills implementation of an SVM classifier model.
 *
 * @author Matthew Pocock
 */
public class SimpleSVMClassifierModel extends AbstractSVMClassifierModel {
  private SVMKernel kernel;
  private double threshold;
  
  private Set itemAlphaSet;
  private Map itemToItemAlpha;
  
  public SimpleSVMClassifierModel(SVMKernel kernel) {
    this.kernel = kernel;
    itemAlphaSet = new HashSet();
    itemToItemAlpha = new HashMap();
  }
  
  public SimpleSVMClassifierModel(SVMKernel kernel, Collection items) {
    this(kernel);
    for(Iterator i = items.iterator(); i.hasNext(); ) {
      addItem(i.next());
    }
  }    

  public SimpleSVMClassifierModel(SVMKernel kernel, SVMTarget target) {
    this(kernel, target.items());
  }
  
  public SVMKernel getKernel() {
    return kernel;
  }
  
  public void setThreshold(double threshold) {
    this.threshold = threshold;
  }
  
  public double getThreshold() {
    return threshold;
  }
  
  public Set items() {
    return itemToItemAlpha.keySet();
  }
  
  public Set itemAlphas() {
    return itemAlphaSet;
  }
  
  public double getAlpha(Object item) {
    return ((ItemValue) itemToItemAlpha.get(item)).getValue();
  }
  
  public void setAlpha(Object item, double alpha) {
    ItemValue iv = (ItemValue) itemToItemAlpha.get(item);
    iv.setValue(alpha);
  }
  
  public void addItem(Object item) {
    ItemValue iv = new SimpleItemValue(item, 0.0);
    itemToItemAlpha.put(item, iv);
    itemAlphaSet.add(iv);
  }
  
  public void addItemAlpha(Object item, double alpha) {
    ItemValue iv = new SimpleItemValue(item, alpha);
    itemToItemAlpha.put(item, iv);
    itemAlphaSet.add(iv);
  }
  
  public void removeItem(Object item) {
    itemToItemAlpha.remove(item);
    itemAlphaSet.remove(item);
  }
  
  public void clear() {
    itemAlphaSet.clear();
    itemToItemAlpha.clear();
  }
}
