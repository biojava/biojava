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
package org.biojava.nbio.structure.symmetry.internal;

import static java.lang.Math.*;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.util.RotationAxis;
import org.biojava.nbio.structure.symmetry.internal.OrderDetector;
import org.biojava.nbio.structure.symmetry.internal.RefinerFailedException;

/**
 * Guesses an order of rotational symmetry from the angle.
 * <p>
 * Improves upon the AngleOrderDetector used in the paper by checking all
 * rotations, not just the base one.
 *
 * @author Spencer Bliven
 * @since 4.2.0
 *
 */
public class AngleOrderDetectorPlus implements OrderDetector {

	private int maxOrder;
	private final double error;
	private boolean normalizeError;

	/**
	 * @param error
	 *            maximum angular error, in radians
	 */
	public AngleOrderDetectorPlus(double angleError) {
		this(8, angleError, false);
	}

	public AngleOrderDetectorPlus(int maxOrder) {
		// PI is the max error for C1
		this(maxOrder, PI, false);
	}

	/**
	 *
	 * @param maxOrder
	 *            maximum order to consider
	 * @param error
	 *            maximum angular error, in radians
	 */
	public AngleOrderDetectorPlus(int maxOrder, double angleError) {
		this(maxOrder, angleError, false);
	}

	/**
	 * Determine order by finding the order (up to the maxOrder) which has the
	 * closest rotation angle to the observed rotation.
	 *
	 * If normalized is false, then the error is taken to be the absolute error
	 * from the closest ideal angle (in radians). If normalized is true, error
	 * is taken to be relative to the fundamental rotation for a given order.
	 * For instance, for an error of .25, C2 order would be accepted for angles
	 * within .25*pi radians of 0 or pi, while C3 order would be acceptable
	 * within .25*2pi/3 radians of 0, 2pi/3, or 4pi/3. In the normalized case,
	 * numbers between 0 and .5 are sensible for error.
	 *
	 * @param maxOrder
	 *            maximum order to consider
	 * @param error
	 *            maximum angular error
	 * @param normalize
	 *            indicates whether error should be normalized by the order
	 */
	public AngleOrderDetectorPlus(int maxOrder, double angleError,
			boolean normalize) {
		super();
		this.maxOrder = maxOrder;
		this.error = angleError;
		this.normalizeError = normalize;
	}

	@Override
	public int calculateOrder(AFPChain afpChain, Atom[] ca)
			throws RefinerFailedException {
		final double tol = 1e-6; // tolerance to floating point errors
		try {
			RotationAxis axis = new RotationAxis(afpChain);
			double theta = axis.getAngle();

			double bestDelta = error;
			int bestOrder = 1;
			for (int order = 1; order <= maxOrder; order++) {
				// Triangle wave starting at 0 with period 2pi/order
				double delta = abs(abs(theta * order / (2 * PI) - .5) % 1.0 - .5);
				// Triangle waves have amplitude 1, so need to un-normalize
				if (!normalizeError)
					delta *= 2 * PI / order;

				if (delta < bestDelta - tol) {
					bestOrder = order;
					bestDelta = delta;
				}
			}
			return bestOrder;
		} catch (Exception e) {
			throw new RefinerFailedException(e);
		}
	}

	@Override
	public String toString() {
		return getClass().getSimpleName() + "[maxOrder=" + maxOrder
				+ ", error=" + error + ", normalizeError=" + normalizeError
				+ "]";
	}

}
