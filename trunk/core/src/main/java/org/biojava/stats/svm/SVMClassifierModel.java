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

import java.util.Set;

/**
 * An SVM classifier model.
 * <p>
 * This is the interface for objects that contain the model for a binary
 * classification task.
 *
 * @author Matthew Pocock
 */
public interface SVMClassifierModel {
  public SVMKernel getKernel();
  
  public void setThreshold(double threshold)
  throws UnsupportedOperationException;
  
  public double getThreshold();
  
  public Set items();
  
  public Set itemAlphas();

  public double getAlpha(Object item);
  
  public void setAlpha(Object item, double alpha)
  throws UnsupportedOperationException;
  
  public void addItem(Object item)
  throws UnsupportedOperationException;
  
  public void addItemAlpha(Object item, double alpha)
  throws UnsupportedOperationException;
  
  public void removeItem(Object item)
  throws UnsupportedOperationException;
  
  public void clear()
  throws UnsupportedOperationException;
  
  public double classify(Object item);
}
