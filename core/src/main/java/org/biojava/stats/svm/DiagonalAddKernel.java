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

import java.util.HashSet;
import java.util.Set;

/**
 * Adds a class specific constant to k(x, x).
 *
 * @author Matthew Pocock
 */
public class DiagonalAddKernel extends NestedKernel {
  private Set posClass;
  private Set negClass;

  {
    posClass = new HashSet();
    negClass = new HashSet();
  }
  
  public void addPos(Object o) {
    posClass.add(o);
  }
  
  public void addNeg(Object o) {
    negClass.add(o);
  }
  
  /**
   * The scale vactor.
   */
  private double lambda = 1.0;
  
  /**
   * Set the scale factor.
   *
   * @param l  the new scale factor
   */
  public void setLambda(double l) {
    this.lambda = l;
  }
  
  /**
   * Retrieve the scale factor.
   *
   * @return the current scale factor
   */
  public double getLambda() {
    return lambda;
  }
  
  /**
   * Return the dot product of a, b.
   * <p>
   * This is equal to
   * <code>k(a, b) + d(a, b) * ||class(a)|| / (||class||)</code>
   * where d(a, b) is zero if a != b, and 1 if a == b. class(a) is the set of all
   * items in the same class as a. class is all items with a classification.
   */
  public double evaluate(Object a, Object b) {
    double dot = getNestedKernel().evaluate(a, b);
    if(a == b) {
      int size = 0;
      if(posClass.contains(a)) {
        size = posClass.size();
      } else if(negClass.contains(a)) {
        size = negClass.size();
      }
      dot += getLambda() * size / (posClass.size() + negClass.size());
    }
    return dot;
  }
  
  public String toString() {
    return
     "DiagonalAdd K(a, b | l, s+, s-, k) = k(a, b) + d[a, b]; d[a, b] = " +
     "{ a != b, 0; a == b, l * {class(a == +), s+; class(a == -), s-} }; k = " +
     getNestedKernel().toString();
  }
}
