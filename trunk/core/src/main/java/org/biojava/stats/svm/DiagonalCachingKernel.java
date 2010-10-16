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

import java.util.HashMap;
import java.util.Map;

/**
 * Caches the leading diagonal of a kernel matrix.
 * <p>
 * Several kernels need to repeatedly access k(x,x) to do things like
 * normalization, or to calculate distances. This kernel wraps k so that these
 * leading diagonal elements do not need to be calculated each time.
 * <p>
 * This kernel is thread-safe. However, care must be taken when setting the
 * nested kernel that no other thread is retrieving values at the same time.
 * This would cause a race condition in which the newly flushed cache may
 * contain a value from the previous kernel.
 *
 * @author Matthew Pocock
 * @author Thomas Down
 */
public class DiagonalCachingKernel extends NestedKernel {
  /**
   * The cache of values.
   */
  transient private Map cache;

  {
    cache = new HashMap();
  }
  
  /**
   * Create a new CachingKernel.
   */
  public DiagonalCachingKernel() {
  }

  /**
   * Creates a new DiagonalCachingKernel that nests k.
   */
  public DiagonalCachingKernel(SVMKernel k) {
    setNestedKernel(k);
  }
  
  /**
   * <p>
   * Set the kernel to nest.
   * </p>
   *
   * <p>
   * This will flush the cache.
   * </p>
   *
   * @param k  the kernel to nest.
   */
  public void setNestedKernel(SVMKernel k) {
    super.setNestedKernel(k);
    synchronized(cache) {
      cache.clear();
    }
  }

  /**
   * <p>
   * Returns the kernel product of two Objects.
   * </p>
   *
   * <p>
   * This returns <code>getNestedKernel.evaluate(x, y)</code>. If
   * <code>x.equals(y)</code> then it will cache the result first time, and do
   * a hash table look up to retrieve the value in subsequent calls.
   * </p>
   */
  public double evaluate(Object x, Object y) {
    if(x.equals(y)) {
      Double d = null;
      synchronized(cache) {
        d = (Double) cache.get(x);
      }
      if (d == null) {
        d = new Double(getNestedKernel().evaluate(x, x));
        synchronized(cache) {
          cache.put(x, d);
        }
      }
      return d.doubleValue();
    } else {
      return getNestedKernel().evaluate(x, y);
    }
  }
  
  public String toString() {
    return getNestedKernel().toString();
  }
}
