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
 * Kernel for support vector machines and related methods.
 * <p>
 * It is hoped that all implementations of SVMKernel will be serializable,
 * as this will allow models to be stored and retrieved without inventing
 * complex data formats.
 *
 * @author Thomas Down.
 */
public interface SVMKernel {
    /**
     * Return the dot product of two vectors in an arbitrary
     * feature space.  In this implementation, the `vectors'
     * can actually be arbitrary objects.
     */
    public double evaluate(Object a, Object b);
}
