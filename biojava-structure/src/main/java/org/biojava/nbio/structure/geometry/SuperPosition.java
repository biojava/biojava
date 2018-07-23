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
package org.biojava.nbio.structure.geometry;

import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;

/**
 * The SuperPosition interface defines and documents the required methods for
 * any superpostion algorithm implementation, so that the input and expected
 * output are uniform.
 * 
 * @author Aleix Lafita
 * @since 5.0.0
 *
 */
public interface SuperPosition {

	/**
	 * Obtain the superposition matrix that minimizes the RMSD between two
	 * arrays of equivalent points.
	 * <p>
	 * The two point arrays have to be of the same length and the order of
	 * points have to be the same, so that a specific position in the one array
	 * is equivalent to the same position in the other array.
	 * 
	 * @param fixed
	 *            point array as reference, onto which the other point array is
	 *            superposed. Original coordinates will not be modified.
	 * @param moved
	 *            point array to which the resulting transformation matrix is
	 *            applied. Original coordinates will not be modified.
	 * @return transformation matrix as a Matrix4d to superpose moved onto fixed
	 *         point arrays
	 */
	public Matrix4d superpose(Point3d[] fixed, Point3d[] moved);

	/**
	 * Transform an array of points so that the coordinates of its points
	 * minimize the RMSD to the other array of equivalent points, and return the
	 * transformation matrix applied.
	 * <p>
	 * The two point arrays have to be of the same length and the order of
	 * points have to be the same, so that a specific position in the one array
	 * is equivalent to the same position in the other array.
	 * 
	 * @param fixed
	 *            point array as reference, onto which the other point array is
	 *            superposed. Original coordinates will not be modified.
	 * @param moved
	 *            point array to which the resulting transformation matrix is
	 *            applied. Original coordinates will be transformed.
	 * @return transformation matrix as a Matrix4d to superpose moved onto fixed
	 *         point arrays
	 */
	public Matrix4d superposeAndTransform(Point3d[] fixed, Point3d[] moved);

	/**
	 * Calculate the RMSD between two arrays of equivalent points that are not
	 * superposed.
	 * <p>
	 * This is equivalent to first superposing the point arrays with
	 * {@link SuperPosition#superposeAndTransform(Point3d[], Point3d[])} and
	 * then calculating the RMSD of the superposed point arrays with
	 * {@link CalcPoint#rmsd(Point3d[], Point3d[])}, but it will be faster when
	 * the transformation matrix is not needed.
	 * <p>
	 * The two point arrays have to be of the same length and the order of
	 * points have to be the same, so that a specific position in the one array
	 * is equivalent to the same position in the other array.
	 * 
	 * @param x
	 *            an array of points. Original coordinates will not be modified.
	 * @param y
	 *            an array of points. Original coordinates will not be modified.
	 * @return the minimum RMSD between the equivalent point arrays (after
	 *         superposition)
	 */
	public double getRmsd(Point3d[] x, Point3d[] y);

}
