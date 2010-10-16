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

import java.io.Serializable;

/**
 * Encapsulates a kernel that wraps another kernel up.
 *
 * @author Matthew Pocock
 */
public abstract class NestedKernel implements SVMKernel, Serializable {
  /**
   * The <span class="type">SVMKernel</span> being wrapped.
   */
  private SVMKernel nested;

  /**
   * Create a new <span class="type">NestedKernel</span>.
   */
  public NestedKernel() {}

  /**
   * Create a new <span class="type">NestedKernel</span> that wraps
   * <span class="arg">k</span>.
   *
   * @param k  the <span class="type">SVMKernel</span> to wrap
   */
  public NestedKernel(SVMKernel k) {
    setNestedKernel(k);
  }

  /**
   * Set the <span class="type">SVMKernel</span> to nest to
   * <span class="arg">k</span>.
   *
   * @param k  the <span class="type">SVMKernel</span> to nest.
   */
  public void setNestedKernel(SVMKernel k) {
    nested = k;
  }

  /**
   * Retrieve the currently nested <span class="type">SVMKernel</span>.
   *
   * @return the nested <span class="type">SVMKernel</span>
   */
  public SVMKernel getNestedKernel() {
    return nested;
  }
}
