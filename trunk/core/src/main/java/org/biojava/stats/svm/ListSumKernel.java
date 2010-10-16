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
import java.util.List;

/**
 * This kernel computes the sum of the dot products between items of two lists
 * at corresponding indexes. We define k(x, null) to be zero for when list
 * elements are null.
 *
 * @author Matthew Pocock
 */
public class ListSumKernel extends NestedKernel {
    public double evaluate(Object a, Object b) {
      List l1 = (List) a;
      List l2 = (List) b;
      
      double dot = 0.0;
      SVMKernel k = getNestedKernel();
      
      Iterator i1 = l1.iterator();
      Iterator i2 = l2.iterator();
      
      while(i1.hasNext() && i2.hasNext()) {
        Object o1 = i1.next();
        Object o2 = i2.next();
        if(o1 != null && o2 != null) {
          dot += k.evaluate(i1.next(), i2.next());
        }
      }
      
      return dot;
    }

    public String toString() {
      return "List kernel K(x, y, k) = sum_i(k(x_i, y_i))"
        + "; k = " + getNestedKernel().toString();
    }
}
