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
 * This kernel implements a three layer neural net. This is calculated as:
 *   tanh(a*k(x,y)+c)
 *
 * @author Matthew Pocock
 */
public class SigmoidKernel implements SVMKernel {
    private double a;
    private double c;
    private SVMKernel kernel;

    public SigmoidKernel() {
      a = 1.0;
      c = 1.0;
      kernel = null;
    }

    public double evaluate(Object a, Object b) {
      return tanh(getMultiplier()*getWrappedKernel().evaluate(a, b)
                  + getConstant());
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

    public SVMKernel getWrappedKernel() {
      return kernel;
    }
    
    public void setWrappedKernel(SVMKernel kernel) {
      this.kernel = kernel;
    }
    
    public String toString() {
      return "Sigmoid kernel K(x, k) = tanh("
        + getMultiplier() + ".k(x) + " + c + ")"
        + ". k = " + getWrappedKernel().toString();
    }
    
    public double tanh(double a) {
      double x = Math.exp(a);
      double y = Math.exp(-a);
      
      return (x - y) / (x + y);
    }
}
