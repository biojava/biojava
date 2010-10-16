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

import org.biojava.bio.BioException;

/**
 * solves y = f(x) = 0 by binary search.
 * Only really suitable for monotonic functions as
 * the method will check that the initial values
 * lie on opposite sides of the X=0 axis.
 *
 * @author David Huen
 */
public class BinarySearch
{

    /**
     * method that will attempt solving the equation.
     *
     * @param min lower bound of search space.
     * @param max upper bound of search space.
     * @param tolerance change in x required to continue iteration.
     * @param obj the class of ComputeObject class representing the equation to be solved.
     */
    public static double solve(double min, double max, double tolerance, ComputeObject obj)
        throws BioException
    {
        // compute initial values
        double x1 = min;
        double y1 = obj.compute(min);
        double x2 = max;
        double y2 = obj.compute(max);

        // validate that function standas some chance of monotonicity
        if ((y1 <  0.0) && (y2 < 0.0)) throw new BioException("Illegal initial range limits.");
        if ((y1 >  0.0) && (y2 > 0.0)) throw new BioException("Illegal initial range limits.");

        // iterate
        while (Math.abs(x1 - x2) > tolerance) {
            // compute a value midway within the current interval
            double newX = 0.5 * (x1 + x2);
            double newY = obj.compute(newX);

            // determine new interval
            if (newY >= 0.0) {
                if (y1 >= 0.0) {
                    y1 = newY;
                    x1 = newX;
                }
                else {
                    y2 = newY;
                    x2 = newX;
                }
            }
            else if (newY < 0.0) {
                if (y1 >= 0.0) {
                    y2 = newY;
                    x2 = newX;
                }
                else {
                    y1 = newY;
                    x1 = newX;
                }
            }
        }

        // converged: return midpoint of interval
        return 0.5 * (x1 + x2);
    }

}

