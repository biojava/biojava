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


/**
 * <p>
 * An abstract implementation of an SVMModel.
 * </p>
 *
 * <p>
 * You only need implement items, itemTargets and getTarget to make a
 * read-only implementation.
 * </p>
 *
 * @author Matthew Pocock
 */
public abstract class AbstractSVMTarget implements SVMTarget {
  public void setTarget(Object item, double target)
  throws UnsupportedOperationException {
    throw new UnsupportedOperationException();
  }
  
  public void addItem(Object item)
  throws UnsupportedOperationException {
    throw new UnsupportedOperationException();
  }
  
  public void addItemTarget(Object item, double target)
  throws UnsupportedOperationException {
    throw new UnsupportedOperationException();
  }
  
  public void removeItem(Object item)
  throws UnsupportedOperationException {
    throw new UnsupportedOperationException();
  }
  
  public void clear()
  throws UnsupportedOperationException {
    throw new UnsupportedOperationException();
  }
}
