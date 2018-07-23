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
import javax.vecmath.Vector3d;

import org.biojava.nbio.structure.jama.Matrix;

/**
 * Utility operations on Point3d.
 * 
 * @author Aleix Lafita
 * @since 5.0.0
 *
 */
public class CalcPoint {

	/** Prevent instantiation */
	private CalcPoint() {
	}

	/**
	 * Center a cloud of points. This means subtracting the {@lin
	 * #centroid(Point3d[])} of the cloud to each point.
	 * 
	 * @param x
	 *            array of points. Point objects will be modified
	 */
	public static void center(Point3d[] x) {
		Point3d center = centroid(x);
		center.negate();
		translate(new Vector3d(center), x);
	}

	/**
	 * Calculate the centroid of the point cloud.
	 * 
	 * @param x
	 *            array of points. Point objects will not be modified
	 * @return centroid as Point3d
	 */
	public static Point3d centroid(Point3d[] x) {
		Point3d center = new Point3d();
		for (Point3d p : x) {
			center.add(p);
		}
		center.scale(1.0 / x.length);
		return center;
	}

	/**
	 * Transform all points with a 4x4 transformation matrix.
	 * 
	 * @param rotTrans
	 *            4x4 transformation matrix
	 * @param x
	 *            array of points. Point objects will be modified
	 */
	public static void transform(Matrix4d rotTrans, Point3d[] x) {
		for (Point3d p : x) {
			rotTrans.transform(p);
		}
	}

	/**
	 * Translate all points with a translation vector.
	 * 
	 * @param trans
	 *            the translation vector to apply
	 * @param x
	 *            array of points. Point objects will be modified
	 */
	public static void translate(Vector3d trans, Point3d[] x) {
		for (Point3d p : x) {
			p.add(trans);
		}
	}

	/**
	 * Clone an array of points.
	 * 
	 * @param x
	 *            original array of points. Point objects will not be modified
	 * @return new array of points, identical clone of x
	 */
	public static Point3d[] clonePoint3dArray(Point3d[] x) {
		Point3d[] clone = new Point3d[x.length];
		for (int i = 0; i < x.length; i++) {
			clone[i] = new Point3d(x[i]);
		}
		return clone;
	}

	/*
	 * Peter can you document this method? TODO
	 * 
	 * @param moved
	 * @param fixed
	 * @return
	 */
	public static Matrix formMatrix(Point3d[] a, Point3d[] b) {
		double xx = 0.0, xy = 0.0, xz = 0.0;
		double yx = 0.0, yy = 0.0, yz = 0.0;
		double zx = 0.0, zy = 0.0, zz = 0.0;

		for (int i = 0; i < a.length; i++) {
			xx += a[i].x * b[i].x;
			xy += a[i].x * b[i].y;
			xz += a[i].x * b[i].z;
			yx += a[i].y * b[i].x;
			yy += a[i].y * b[i].y;
			yz += a[i].y * b[i].z;
			zx += a[i].z * b[i].x;
			zy += a[i].z * b[i].y;
			zz += a[i].z * b[i].z;
		}

		double[][] f = new double[4][4];
		f[0][0] = xx + yy + zz;
		f[0][1] = zy - yz;
		f[1][0] = f[0][1];
		f[1][1] = xx - yy - zz;
		f[0][2] = xz - zx;
		f[2][0] = f[0][2];
		f[1][2] = xy + yx;
		f[2][1] = f[1][2];
		f[2][2] = yy - zz - xx;
		f[0][3] = yx - xy;
		f[3][0] = f[0][3];
		f[1][3] = zx + xz;
		f[3][1] = f[1][3];
		f[2][3] = yz + zy;
		f[3][2] = f[2][3];
		f[3][3] = zz - xx - yy;

		return new Matrix(f);
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
		
		if (x.length != y.length) {
			throw new IllegalArgumentException(
					"Point arrays are not of the same length.");
		}
		
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
	 * 
	 * @param y
	 * 
	 * @return
	 */
	public static double GTSlikeScore(Point3d[] x, Point3d[] y) {
		
		if (x.length != y.length) {
			throw new IllegalArgumentException(
					"Point arrays are not of the same length.");
		}
		
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
		
		if (x.length != y.length) {
			throw new IllegalArgumentException(
					"Point arrays are not of the same length.");
		}
		
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
		
		if (x.length != y.length) {
			throw new IllegalArgumentException(
					"Point arrays are not of the same length.");
		}
		
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

	/*
	 * Needs documentation!
	 * 
	 * @param x
	 * 
	 * @param y
	 * 
	 * @param maxDistance
	 * 
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
