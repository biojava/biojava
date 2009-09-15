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
/*
 * @(#)CachingKernel.java      0.1 00/01/20
 *
 * By Thomas Down <td2@sanger.ac.uk>
 */

package org.biojava.stats.svm;

import java.util.HashMap;
import java.util.Map;

/**
 * <p>
 * Caches the results of a nested kernel so that k(a, b) need only be calculated
 * once.
 * </p>
 *
 * <p>
 * This kernel is thread-safe. However, care must be taken when setting the
 * nested kernel that no other thread is retrieving values at the same time.
 * This would cause a race condition in which the newly flushed cache may
 * contain a value from the previous kernel.
 * </p>
 *
 * @author Thomas Down
 * @author Matthew Pocock
 */
public class CachingKernel extends NestedKernel {
    transient private Map cache;
    
    public CachingKernel() {
    }
    
    public CachingKernel(SVMKernel k) {
      super(k);
    }

    public void setNestedKernel(SVMKernel k) {
      super.setNestedKernel(k);
      if(cache == null) {
        cache = new HashMap();
      }
      cache.clear();
    }

    public double evaluate(Object x, Object y) {
      ObjectPair op = new ObjectPair(x, y);
      Double d = null;
      synchronized (cache) {
        d = (Double) cache.get(op);
      }
      if (d == null) {
        d = new Double(getNestedKernel().evaluate(x, y));
        synchronized (cache) {
          cache.put(op, d);
        }
      }
      return d.doubleValue();
    }

    public String toString() {
      return getNestedKernel().toString();
    }
    
    private static class ObjectPair {
      Object a;
      Object b;

      public ObjectPair(Object a, Object b) {
        this.a = a;
        this.b = b;
      }

      public boolean equals(Object x) {
        if (! (x instanceof ObjectPair))
          return false;
        ObjectPair op = (ObjectPair) x;
        return ((op.a == a && op.b == b) || 
		            (op.a == b && op.b == a));
      }

      public int hashCode() {
        return a.hashCode() + b.hashCode();
      }
    }
}
