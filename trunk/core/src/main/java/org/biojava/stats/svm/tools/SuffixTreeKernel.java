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

package org.biojava.stats.svm.tools;

import java.io.Serializable;
import java.util.BitSet;

import org.biojava.bio.symbol.SuffixTree;
import org.biojava.stats.svm.SVMKernel;

/**
 * Computes the dot-product of two suffix-trees as the sum of the products
 * of the counts of all nodes they have in common.
 * <p>
 * This implementation allows you to scale the sub-space for each word length
 * independently.
 *
 * @author Matthew Pocock
 */
public class SuffixTreeKernel implements SVMKernel, Serializable {
  /**
   * The <span class="type">DepthScaler</span> that will scale each sub-space.
   * This defaults to <span class="type">UniformScaler</type>.
   */
  private DepthScaler depthScaler = new UniformScaler();
  
  /**
   * Retrieve the current <span class="type">DepthScaler</span>.
   *
   * @return the current <span class="type">DepthScaler</span>
   */
  public DepthScaler getDepthScaler() {
    return depthScaler;
  }
  
  /**
   * Change the current <span class="type">DepthScaler</span> to
   * <span class="arg">depthScaler</span>.
   *
   * @param depthScaler  the new <span class="type">DepthScaler</span> to use
   */
  public void setDepthScaler(DepthScaler depthScaler) {
    this.depthScaler = depthScaler;
  }
  
  /**
   * Calculate the dot product between the
   * <span class="type">SuffixTree</span>s <span class="arg">a</span> and
   * <span class="arg">b</span>.
   * <p>
   * This is the sum of the dot products of each subspace for a given word
   * length. Each subspace is scaled using the
   * <span class="type">DepthScaler</span> returned by
   * <span class="method">getDepthScaler</span>.
   *
   * @param a  the first <span class="type">Object</span>
   * @param b  the second <span class="type">Object</span>
   * @return <span class="arg">a</span>.<span class="arg">b</span>
   * @throws <span class="type">ClassCastException</span> if either
   *         <span class="arg">a</span> or <span class="arg">b</span> are not
   *         castable to <span class="type">SuffixTree</span>
   */
  public double evaluate(Object a, Object b) {
    SuffixTree st1 = (SuffixTree) a;
    SuffixTree st2 = (SuffixTree) b;
    SuffixTree.SuffixNode n1 = st1.getRoot();
    SuffixTree.SuffixNode n2 = st2.getRoot();
      
    return dot(st1, n1, st2, n2, st1.getAlphabet().size(), 0);
  }
  
  /**
   * Recursive method to compute the dot product of the
   * <span class="type">SuffixTree.SuffixNode</span>s
   * <span class="arg">n1</span> and <span class="arg">n2</span>.
   * <p>
   * This scales <span class="arg">n1</span>.
   * <span class="method">getNumber</span><code>()</code> *
   * <span class="arg">n2</span>.
   * <span class="method">getNumber</span><code>()</code>
   * by <span class="const">this</span>.<span class="method">getDepthScaler</span>
   * (<span class="arg">depth</span>), and then returns the sum of this and the
   * dot products for all children of the suffix nodes.
   */
  private double dot(SuffixTree st1,
		     SuffixTree.SuffixNode n1,
		     SuffixTree st2,
                     SuffixTree.SuffixNode n2,
		     int size,
		     int depth)
  {
    double scale = getDepthScaler().getScale(depth);
    double dot = n1.getNumber() * n2.getNumber() * scale * scale;
    for(int i = 0; i < size; i++) {
      if(n1.hasChild(i) && n2.hasChild(i)) {
        dot += dot(st1, st1.getChild(n1, i), st2, st2.getChild(n2, i), size, depth+1);
      }
    }
    return dot;
  }
    
  public String toString() {
    return new String("Suffix tree kernel");
  }
  
  /**
   * Encapsulates the scale factor to apply at a given depth.
   *
   * @author Matthew Pocock
   */
  public interface DepthScaler {
    /**
     * Retrieve the scaling factor at a given depth
     *
     * @param depth  word length
     * @return the scaling factor for the subspace at that length
     */
    double getScale(int depth);
  }
  
  /**
   * Scales by 4^depth - equivalent to dividing by a probablistic flatt prior
   * null model
   *
   * @author Matthew Pocock
   */
  public static class NullModelScaler implements DepthScaler, Serializable {
    public double getScale(int depth) {
      return Math.pow(4.0, (double) depth);
    }
  }
  
  /**
   * Scale all depths by 1.0
   *
   * @author Matthew Pocock
   */
  public static class UniformScaler implements DepthScaler, Serializable {
    public double getScale(int depth) {
      return 1.0;
    }
  }
  
  /**
   * Scale using a <span class="type">BitSet</span> to allow/disallow depths.
   *
   * @author Matthew Pocock
   */
  public static class SelectionScalar implements DepthScaler, Serializable {
    private BitSet bSet;
    
    /**
     * Make a new <span class="type">SelectionScalar</span> that masks in different
     * depths.
     *
     * @param bSet  the mask for which depths to allow
     */
    public SelectionScalar(BitSet bSet) {
      this.bSet = new BitSet();
      this.bSet.or(bSet);
    }
    
    /**
     * @return 1.0 or 0.0 depending on whether the bit at
     *         <span class="arg">depth</span> is set or not
     */
    public double getScale(int depth) {
      if(bSet.get(depth)) {
        return 1.0;
      } else {
        return 0.0;
      }
    }
  }
  
  /**
   * Scale using a multiple of two <span class="type">DepthScaler</span>s.
   *
   * @author Matthew Pocock
   */
  public static class MultipleScalar implements DepthScaler, Serializable {
    private DepthScaler a;
    private DepthScaler b;
    
    public MultipleScalar(DepthScaler a, DepthScaler b) {
      this.a = a;
      this.b = b;
    }
    
    public double getScale(int depth) {
      return a.getScale(depth) * b.getScale(depth);
    }
  }
}
