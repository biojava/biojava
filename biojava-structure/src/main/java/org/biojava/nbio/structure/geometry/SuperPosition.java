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

import javax.vecmath.*;

/**
 * The SuperPosition is a static Class that contains methods to superpose arrays
 * of Points in 3D.
 * 
 * @author Peter Rose
 * 
 */
public final class SuperPosition {

	/** Prevent instantiation */
	private SuperPosition() {
	}

	/**
	 * This method superposes x onto y, so it transforms the coordinates of x.
	 * The translation component of the matrix is wrong.
	 * 
	 * @param x
	 *            point array, point coordinates will be modified
	 * @param y
	 *            point array, point coordinates will not be modified
	 * @return the Matrix4d used for superposition, without translation
	 * @deprecated use {@link #superposeWithTranslation(Point3d[], Point3d[])}
	 */
	@Deprecated
	public static Matrix4d superpose(Point3d[] x, Point3d[] y) {
		// superpose x onto y
		Point3d[] ref = CalcPoint.clonePoint3dArray(y);

		Point3d ytrans = CalcPoint.centroid(ref);
		ytrans.negate();
		CalcPoint.translate(ytrans, ref);

		CalcPoint.center(x);

		// calculate quaternion from relative orientation
		Quat4d q = UnitQuaternions.relativeOrientation(x, ref);

		Matrix4d rotTrans = new Matrix4d();
		rotTrans.set(q);

		// set translational component of transformation matrix TODO wrong
		ytrans.negate();
		rotTrans.setTranslation(new Vector3d(ytrans));

		// tranform coordinates
		CalcPoint.transform(rotTrans, x);

		return rotTrans;
	}

	/**
	 * This method superposes x onto y, so it transforms the coordinates of x.
	 * 
	 * @param x
	 *            point array, point coordinates will be modified
	 * @param y
	 *            point array, point coordinates will not be modified
	 * @return the Matrix4d used for superposition
	 */
	public static Matrix4d superposeWithTranslation(Point3d[] x, Point3d[] y) {

		// translate to origin
		Point3d[] xref = CalcPoint.clonePoint3dArray(x);
		Point3d xtrans = CalcPoint.centroid(xref);
		xtrans.negate();
		CalcPoint.translate(xtrans, xref);

		Point3d[] yref = CalcPoint.clonePoint3dArray(y);
		Point3d ytrans = CalcPoint.centroid(yref);
		ytrans.negate();
		CalcPoint.translate(ytrans, yref);

		// calculate rotational component (rotation around origin)
		Quat4d q = UnitQuaternions.relativeOrientation(xref, yref);
		Matrix4d rotTrans = new Matrix4d();
		rotTrans.set(q);

		// combine with x -> origin translation
		Matrix4d trans = new Matrix4d();
		trans.setIdentity();
		trans.setTranslation(new Vector3d(xtrans));
		rotTrans.mul(rotTrans, trans);

		// combine with origin -> y translation
		ytrans.negate();
		Matrix4d transInverse = new Matrix4d();
		transInverse.setIdentity();
		transInverse.setTranslation(new Vector3d(ytrans));
		rotTrans.mul(transInverse, rotTrans);

		// transform x coordinates onto y coordinate frame
		CalcPoint.transform(rotTrans, x);

		return rotTrans;
	}

	/**
	 * This method superposes x onto y, so it transforms the coordinates of x.
	 * It requires the points to be centered (use
	 * {@link CalcPoint#center(Point3d[])}, because the transformation matrix
	 * will not have translation component (faster).
	 * 
	 * @param x
	 *            point array, point coordinates will be modified
	 * @param y
	 *            point array, point coordinates will not be modified
	 * @return the Matrix4d used for superposition, without translation
	 */
	public static Matrix4d superposeAtOrigin(Point3d[] x, Point3d[] y) {
		Quat4d q = UnitQuaternions.relativeOrientation(x, y);

		Matrix4d rotTrans = new Matrix4d();
		rotTrans.set(q);

		CalcPoint.transform(rotTrans, x);

		return rotTrans;
	}

	/**
	 * This method superposes x onto y, so it transforms the coordinates of x.
	 * It requires the points to be centered (use
	 * {@link CalcPoint#center(Point3d[])}, because the transformation matrix
	 * will not have translation component (faster).
	 * 
	 * @param x
	 *            point array, point coordinates will be modified
	 * @param y
	 *            point array, point coordinates will not be modified
	 * @param axisAngle
	 *            to be filled with rotation axis
	 * @return the Matrix4d used for superposition, without translation
	 */
	public static Matrix4d superposeAtOrigin(Point3d[] x, Point3d[] y,
			AxisAngle4d axisAngle) {
		Quat4d q = UnitQuaternions.relativeOrientation(x, y);
		Matrix4d rotTrans = new Matrix4d();
		rotTrans.setIdentity();
		rotTrans.set(q);
		axisAngle.set(q);
		Vector3d axis = new Vector3d(axisAngle.x, axisAngle.y, axisAngle.z);
		if (axis.lengthSquared() < 1.0E-6) {
			// System.err.println("Error: SuperPosition.superposeAtOrigin: axis vector undefined!");
			axisAngle.x = 0;
			axisAngle.y = 0;
			axisAngle.z = 1;
			axisAngle.angle = 0;
		} else {
			axis.normalize();
			axisAngle.x = axis.x;
			axisAngle.y = axis.y;
			axisAngle.z = axis.z;
		}
		CalcPoint.transform(rotTrans, x);

		return rotTrans;
	}

	/**
	 * Calculate the RMSD of two point arrays, already superposed.
	 * 
	 * @param x
	 *            array of points superposed to y
	 * @param y
	 *            array of points superposed to x
	 * @return RMSD
	 */
	public static double rmsd(Point3d[] x, Point3d[] y) {
		double sum = 0.0;
		for (int i = 0; i < x.length; i++) {
			sum += x[i].distanceSquared(y[i]);
		}
		return Math.sqrt(sum / x.length);
	}

	/*
	 * Needs documentation!
	 * 
	 * @param x
	 * 
	 * @param y
	 * 
	 * @return
	 */
	public static double rmsdMin(Point3d[] x, Point3d[] y) {
		double sum = 0.0;
		for (int i = 0; i < x.length; i++) {
			double minDist = Double.MAX_VALUE;
			for (int j = 0; j < y.length; j++) {
				minDist = Math.min(minDist, x[i].distanceSquared(y[j]));
			}
			sum += minDist;
		}
		return Math.sqrt(sum / x.length);
	}

	/**
	 * Returns the TM-Score for two superimposed sets of coordinates Yang Zhang
	 * and Jeffrey Skolnick, PROTEINS: Structure, Function, and Bioinformatics
	 * 57:702–710 (2004)
	 * 
	 * @param x
	 *            coordinate set 1
	 * @param y
	 *            coordinate set 2
	 * @param lengthNative
	 *            total length of native sequence
	 * @return
	 */
	public static double TMScore(Point3d[] x, Point3d[] y, int lengthNative) {
		double d0 = 1.24 * Math.cbrt(x.length - 15.0) - 1.8;
		double d0Sq = d0 * d0;

		double sum = 0;
		for (int i = 0; i < x.length; i++) {
			sum += 1.0 / (1.0 + x[i].distanceSquared(y[i]) / d0Sq);
		}

		return sum / lengthNative;
	}

	/*
	 * Needs documentation!
	 * 
	 * @param x
	 * @param y
	 * @return
	 */
	public static double GTSlikeScore(Point3d[] x, Point3d[] y) {
		int contacts = 0;

		for (Point3d px : x) {
			double minDist = Double.MAX_VALUE;

			for (Point3d py : y) {
				minDist = Math.min(minDist, px.distanceSquared(py));
			}

			if (minDist > 64)
				continue;
			contacts++;

			if (minDist > 16)
				continue;
			contacts++;

			if (minDist > 4)
				continue;
			contacts++;

			if (minDist > 1)
				continue;
			contacts++;
		}

		return contacts * 25.0 / x.length;
	}

	/*
	 * Needs documentation!
	 * @param x
	 * @param y
	 * @param maxDistance
	 * @return
	 */
	public static int contacts(Point3d[] x, Point3d[] y, double maxDistance) {
		int contacts = 0;
		for (int i = 0; i < x.length; i++) {
			double minDist = Double.MAX_VALUE;
			for (int j = 0; j < y.length; j++) {
				minDist = Math.min(minDist, x[i].distanceSquared(y[j]));
			}
			if (minDist < maxDistance * maxDistance) {
				contacts++;
			}
		}
		return contacts;
	}

}
