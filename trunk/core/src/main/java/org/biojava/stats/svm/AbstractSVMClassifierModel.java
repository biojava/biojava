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

import java.util.Iterator;

/**
 * <p>
 * Abstract implementation of SVMClassifierModel.
 * </p>
 *
 * <p>
 * To implement a read-only implementation, you need only implement
 * getThreshold and getAlpha.
 * </p>
 *
 * @author Matthew Pocock
 */
public abstract class AbstractSVMClassifierModel implements SVMClassifierModel {
  public void setThreshold()
  throws UnsupportedOperationException {
    throw new UnsupportedOperationException(
      "setThreshold not supported by " + getClass()
    );
  }
  
  public void setAlpha(Object item, double alpha)
  throws UnsupportedOperationException {
    throw new UnsupportedOperationException(
      "setAlpha not supported by " + getClass()
    );
  }
  
  public void addItem(Object item)
  throws UnsupportedOperationException {
    throw new UnsupportedOperationException(
      "addItem not supported by " + getClass()
    );
  }
  
  public void addItemAlpha(Object item, double alpha)
  throws UnsupportedOperationException {
    throw new UnsupportedOperationException(
      "addItemAlpha not supported by " + getClass()
    );
  }
  
  public void removeItem(Object item)
  throws UnsupportedOperationException {
    throw new UnsupportedOperationException(
      "removeItem not supported by " + getClass()
     );
  }

  
  public void clear()
  throws UnsupportedOperationException {
    throw new UnsupportedOperationException(
      "clear not supported by " + getClass()
    );
  }
  
  public double classify(Object item) {
    double delta=0;
    SVMKernel kernel = getKernel();
    double threshold = getThreshold();
    for(Iterator i = itemAlphas().iterator(); i.hasNext(); ) {
      ItemValue itemValue = (ItemValue) i.next();
      double alpha = itemValue.getValue();
	    if (alpha != 0)
      delta += alpha * kernel.evaluate(itemValue.getItem(), item);
    }
    return delta - threshold;
  }
}
