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
 * This kernel computes all possible products of order features in feature
 * space. This is done by computing (a.k(i,j) + c)^order for some other kernel k
 * that defines a dot product in some feature space.
 *
 * @author Thomas Down
 * @author Matthew Pocock
 */
public class PolynomialKernel extends NestedKernel {
    private double order;
    private double a;
    private double c;

    public PolynomialKernel() {
      this(null, 3.0, 1.0, 1.0);
    }

    public PolynomialKernel(SVMKernel nested, double order, double a, double c) {
      super(nested);
      this.order = order;
      this.a = a;
      this.c = c;
    }
    
    public double evaluate(Object a, Object b) {
      return Math.pow(getMultiplier()*getNestedKernel().evaluate(a, b)
                      + getConstant(),
                      getOrder());
    }

    public double getOrder() {
      return order;
    }

    public void setOrder(double o) {
      this.order = o;
    }

    public double getConstant() {
      return c;
    }

    public void setConstant(double c) {
      this.c = c;
    }
    
    public double getMultiplier() {
      return a;
    }
    
    public void setMultiplier(double m) {
      this.a = m;
    }

    public String toString() {
      return "Polynomial kernel K(x, y | k) = ("
        + getMultiplier() + " * k(x, y) + " + c + ")^" + order
        + ". k = " + getNestedKernel().toString();
    }
}
