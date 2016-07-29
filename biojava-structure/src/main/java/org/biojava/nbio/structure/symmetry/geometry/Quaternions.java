package org.biojava.nbio.structure.symmetry.geometry;

import javax.vecmath.Point3d;
import javax.vecmath.Quat4d;

/**
 * Quaternions is a static Class that contains methods for calculating and using
 * quaternions.
 * <p>
 * A Quaternion is a four-dimensional vector used to describe a
 * three-dimensional attitude representation (axis and angle of rotation).
 * 
 * @author Aleix Lafita
 * @since 5.0.0
 *
 */
public class Quaternions {

	/** Prevent instantiation */
	private Quaternions() {
	}

	/**
	 * The orientation metric is obtained by comparing the quaternion
	 * orientations of two superposed sets of points in 3D.
	 * <p>
	 * First, the quaternion orientation of each set of points is calculated
	 * using their principal axes. Then, the two quaternions are compared using
	 * the formula: d(q1,q2) = arccos(|q1*q2|)
	 * <p>
	 * The formula is taken from: Huynh, D. Q. (2009). Metrics for 3D rotations:
	 * comparison and analysis. Journal of Mathematical Imaging and Vision,
	 * 35(2), 155–164. http://doi.org/10.1007/s10851-009-0161-2
	 * 
	 * @return the quaternion orientation metric
	 */
	public static double orientationMetric(Point3d[] a, Point3d[] b) {

		Quat4d qa = orientation(a);
		Quat4d qb = orientation(b);

		qa.mul(qb);

		double score;

		return 0.0;
	}

	/**
	 * The orientation represents the rotation of the principal axes with
	 * respect to the axes of the coordinate system (unit vectors [1,0,0],
	 * [0,1,0] and [0,0,1]). The rotation can be expressed as a quaternion.
	 * 
	 * @param points
	 *            array of Point3d
	 * @return the orientation as a quaternion
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
	 * Return the length of the quaternion (the norm, the magnitude). The length
	 * of the quaternion is obtained by multiplying by its conjugate and taking
	 * the square root of the sum of terms.
	 * 
	 * @param q
	 *            quaternion as Quat4d object
	 * @return the length of the quaterion
	 */
	public static double length(Quat4d q) {
		return 0.0;
	}

}
