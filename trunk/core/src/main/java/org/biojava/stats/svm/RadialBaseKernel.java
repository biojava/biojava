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
 * This kernel computes the radial base kernel that corresponds to a gausian
 * distribution. 
 * <p>
 * The formula for this is <code>exp( -||a - b|| / (2* width ^ 2))</code>. The
 * term a-b can be represented in an arbitrary feature space by using a nested
 * kernel k, and becomes <code>k(a, a) + k(b, b) - 2 * k(a, b)</code>.
 * <p>
 * As k(x, x) is required repeatedly, I suggest using a DiagonalCachingKernel as
 * the immediately nested kernel function.
 *
 * @author Matthew Pocock
 */
public class RadialBaseKernel extends NestedKernel {
    private double width;

    public RadialBaseKernel() {
      this(null, 1.0);
    }

    public RadialBaseKernel(SVMKernel nested, double width) {
      super(nested);
      this.width = width;
    }
    
    public double evaluate(Object a, Object b) {
      SVMKernel k = getNestedKernel();
      double w = getWidth();
      return Math.exp(-Math.abs(2.0 * k.evaluate(a, b) - k.evaluate(a, a) -
                                k.evaluate(b, b)
                               ) / ( 2.0 * w * w ));
    }

    public double getWidth() {
      return width;
    }

    public void setWidth(double width) {
      this.width = width;
    }
    
    public String toString() {
      return "Radial base kernel K(x, k) = exp(-abs(k(x,x) - k(y,y)) / (2*"
        + getWidth() + "^2)"
        + "; k = " + getNestedKernel().toString();
    }
}
