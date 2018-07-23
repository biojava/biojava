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

import javax.vecmath.AxisAngle4d;
import javax.vecmath.Point3d;
import javax.vecmath.Quat4d;

import org.biojava.nbio.structure.jama.EigenvalueDecomposition;
import org.biojava.nbio.structure.jama.Matrix;

/**
 * UnitQuaternions is a static Class that contains methods for calculating and
 * using unit quaternions. It assumes the use of {@link Quat4d} Class from
 * vecmath to represent the unit quaternions, and it also implements some of the
 * basic methods that the library is missing.
 * <p>
 * A Unit Quaternion is a four-dimensional vector used to describe a
 * three-dimensional attitude representation (axis and angle of rotation). By
 * definition, unit quaternions are always normalized, so their length is always
 * 1.
 * 
 * @author Aleix Lafita
 * @since 5.0.0
 *
 */
public class UnitQuaternions {

	/** Prevent instantiation */
	private UnitQuaternions() {
	}

	/**
	 * The orientation metric is obtained by comparing the quaternion
	 * orientations of the principal axes of each set of points in 3D.
	 * <p>
	 * First, the quaternion orientation of each set of points is calculated
	 * using their principal axes with {@link #orientation(Point3d[])}. Then,
	 * the two quaternions are compared using the method
	 * {@link #orientationMetric(Quat4d, Quat4d)}.
	 * <p>
	 * A requisite for this method to work properly is that both sets of points
	 * have to define the same shape (or very low RMSD), otherwise some of the
	 * principal axes might change or be inverted, resulting in an unreliable
	 * metric. For shapes with some deviations in their shape, use the metric
	 * {@link #orientationAngle(Point3d[], Point3d[])}.
	 * 
	 * @param a
	 *            array of Point3d
	 * @param b
	 *            array of Point3d
	 * @return the quaternion orientation metric
	 */
	public static double orientationMetric(Point3d[] a, Point3d[] b) {

		Quat4d qa = orientation(a);
		Quat4d qb = orientation(b);

		return orientationMetric(qa, qb);
	}

	/**
	 * The orientation metric is obtained by comparing two unit quaternion
	 * orientations.
	 * <p>
	 * The two quaternions are compared using the formula: d(q1,q2) =
	 * arccos(|q1*q2|). The range of the metric is [0, Pi/2], where 0 means the
	 * same orientation and Pi/2 means the opposite orientation.
	 * <p>
	 * The formula is taken from: Huynh, D. Q. (2009). Metrics for 3D rotations:
	 * comparison and analysis. Journal of Mathematical Imaging and Vision,
	 * 35(2), 155–164. http://doi.org/10.1007/s10851-009-0161-2
	 * 
	 * @param q1
	 *            quaternion as Quat4d object
	 * @param q2
	 *            quaternion as Quat4d object
	 * @return the quaternion orientation metric
	 */
	public static double orientationMetric(Quat4d q1, Quat4d q2) {
		return Math.acos(Math.abs(dotProduct(q1, q2)));
	}

	/**
	 * The orientation represents the rotation of the principal axes with
	 * respect to the axes of the coordinate system (unit vectors [1,0,0],
	 * [0,1,0] and [0,0,1]).
	 * <p>
	 * The orientation can be expressed as a unit quaternion.
	 * 
	 * @param points
	 *            array of Point3d
	 * @return the orientation of the point cloud as a unit quaternion
	 */
	public static Quat4d orientation(Point3d[] points) {

		MomentsOfInertia moi = new MomentsOfInertia();

		for (Point3d p : points)
			moi.addPoint(p, 1.0);

		// Convert rotation matrix to quaternion
		Quat4d quat = new Quat4d();
		quat.set(moi.getOrientationMatrix());

		return quat;
	}

	/**
	 * Calculate the rotation angle component of the input unit quaternion.
	 * 
	 * @param q
	 *            unit quaternion Quat4d
	 * @return the angle in radians of the input quaternion
	 */
	public static double angle(Quat4d q) {
		AxisAngle4d axis = new AxisAngle4d();
		axis.set(q);
		return axis.angle;
	}

	/**
	 * The angle of the relative orientation of the two sets of points in 3D.
	 * Equivalent to {@link #angle(Quat4d)} of the unit quaternion obtained by
	 * {@link #relativeOrientation(Point3d[], Point3d[])}.
	 * <p>
	 * The arrays of points need to be centered at the origin. To center the
	 * points use {@link CalcPoint#center(Point3d[])}.
	 * 
	 * @param fixed
	 *            array of Point3d, centered at origin. Original coordinates
	 *            will not be modified.
	 * @param moved
	 *            array of Point3d, centered at origin. Original coordinates
	 *            will not be modified.
	 * @return the angle in radians of the relative orientation of the points,
	 *         angle to rotate moved to bring it to the same orientation as
	 *         fixed.
	 */
	public static double orientationAngle(Point3d[] fixed, Point3d[] moved) {
		Quat4d q = relativeOrientation(fixed, moved);
		return angle(q);
	}

	/**
	 * The angle of the relative orientation of the two sets of points in 3D.
	 * Equivalent to {@link #angle(Quat4d)} of the unit quaternion obtained by
	 * {@link #relativeOrientation(Point3d[], Point3d[])}.
	 * 
	 * @param fixed
	 *            array of Point3d. Original coordinates will not be modified.
	 * @param moved
	 *            array of Point3d. Original coordinates will not be modified.
	 * @param centered
	 *            true if the points are already centered at the origin
	 * @return the angle in radians of the relative orientation of the points,
	 *         angle to rotate moved to bring it to the same orientation as
	 *         fixed.
	 */
	public static double orientationAngle(Point3d[] fixed, Point3d[] moved,
			boolean centered) {
		if (!centered) {
			fixed = CalcPoint.clonePoint3dArray(fixed);
			moved = CalcPoint.clonePoint3dArray(moved);
			CalcPoint.center(fixed);
			CalcPoint.center(moved);
		}
		return orientationAngle(fixed, moved);
	}

	/**
	 * Calculate the relative quaternion orientation of two arrays of points.
	 * 
	 * @param fixed
	 *            point array, coordinates will not be modified
	 * @param moved
	 *            point array, coordinates will not be modified
	 * @return a unit quaternion representing the relative orientation, to
	 *         rotate moved to bring it to the same orientation as fixed.
	 */
	public static Quat4d relativeOrientation(Point3d[] fixed, Point3d[] moved) {
		Matrix m = CalcPoint.formMatrix(moved, fixed); // inverse
		EigenvalueDecomposition eig = m.eig();
		double[][] v = eig.getV().getArray();
		Quat4d q = new Quat4d(v[1][3], v[2][3], v[3][3], v[0][3]);
		q.normalize();
		q.conjugate();
		return q;
	}

	/**
	 * Compute the dot (inner) product of two quaternions.
	 * 
	 * @param q1
	 *            quaternion as Quat4d object
	 * @param q2
	 *            quaternion as Quat4d object
	 * @return the value of the quaternion dot product
	 */
	public static double dotProduct(Quat4d q1, Quat4d q2) {
		return q1.x * q2.x + q1.y * q2.y + q1.z * q2.z + q1.w * q2.w;
	}

}
