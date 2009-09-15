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

package org.biojava.utils.math;


/**
 * interface for classes that return a single
 * double precision value for a single double
 * precision argument.
 * Used to represent equations of type y = f(x) = 0.
 *
 * @author David Huen
 * @since 1.22
 */
public interface ComputeObject
{
    /**
     * workhorse method for this class.
     * computes f(x) for given x.
     */
    public double compute(double arg);
}
