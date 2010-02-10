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
 * Classic linear kernel. It actualy just delegates to SparseVector.kernel.
 *
 * @deprecated Just use SparseVector.kernel instead...
 * @author Matthew Pocock
 */
public class LinearKernel implements SVMKernel, Serializable {
    /**
     * The linear kernel is equal to the dot product of a and b.
     */
    public double evaluate(Object a, Object b) {
      return SparseVector.kernel.evaluate(a, b);
    }
}
