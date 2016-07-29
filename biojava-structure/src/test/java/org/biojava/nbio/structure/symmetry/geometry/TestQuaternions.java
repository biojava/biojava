package org.biojava.nbio.structure.symmetry.geometry;

import static org.junit.Assert.*;

import javax.vecmath.AxisAngle4d;
import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;
import javax.vecmath.Quat4d;
import javax.vecmath.Vector3d;

import org.junit.Test;

/**
 * Test the methods in the {@link Quaternions} class.
 * 
 * @author Aleix Lafita
 * @since 5.0.0
 *
 */
public class TestQuaternions {

	/**
	 * Test {@link Quaternions#orientation(javax.vecmath.Point3d[])}.
	 * <p>
	 * Tests the identity orientation, orientation around one coordinate axis
	 * and orientation around a non-coordinate axis.
	 */
	@Test
	public void testOrientation() {

		// Cloud of points colinear with the coordinate axes centroid at orig
		// Put more mass to x, then y, then z - give correct order of PCA axes
		Point3d[] cloud = new Point3d[180];
		for (int p = 0; p < 60; p++) {

			// Two points, one positive one negative
			double[] coords = { 0, 0, 0 };
			double[] coords_neg = { 0, 0, 0 };

			int value = 60 - p;
			int index = p / 20;

			coords[index] = 2 * value;
			coords_neg[index] = -value;

			cloud[3 * p] = new Point3d(coords);
			cloud[3 * p + 1] = new Point3d(coords_neg);
			cloud[3 * p + 2] = new Point3d(coords_neg);
		}

		Quat4d orientation = Quaternions.orientation(cloud);
		orientation.normalize();

		// No rotation is equivalent to a quaternion with scalar 1 and rest 0
		assertEquals(orientation.x, 0.0, 0.01);
		assertEquals(orientation.y, 0.0, 0.01);
		assertEquals(orientation.z, 0.0, 0.01);
		assertEquals(Math.abs(orientation.w), 1.0, 0.01);

		// Rotate cloud 90 degrees through x axis and recover orientation
		AxisAngle4d axis90x = new AxisAngle4d(new Vector3d(1, 0, 0), 1.57079633);
		Matrix4d mat90x = new Matrix4d();
		mat90x.set(axis90x);

		Point3d[] cloud2 = SuperPosition.clonePoint3dArray(cloud);
		SuperPosition.transform(mat90x, cloud2);

		orientation = Quaternions.orientation(cloud2);
		orientation.normalize();
		AxisAngle4d orientaxis = new AxisAngle4d();
		orientaxis.set(orientation);

		// Angle of rotation 90 degrees around x axis
		assertEquals(Math.abs(orientaxis.x), 1.0, 0.01);
		assertEquals(orientaxis.y, 0.0, 0.01);
		assertEquals(orientaxis.z, 0.0, 0.01);
		assertEquals(orientaxis.angle, axis90x.angle, 0.01);

		// Now try a rotation through a non-coordinate axis
		Quat4d quat = new Quat4d(1.0, 1.0, 1.0, 1.0);

		Matrix4d mat = new Matrix4d();
		mat.set(quat);

		SuperPosition.transform(mat, cloud);

		orientation = Quaternions.orientation(cloud);
		orientation.normalize();

		// Test recovering the quaternion (q and -q same rotation)
		assertEquals(Math.abs(orientation.x), quat.x, 0.01);
		assertEquals(Math.abs(orientation.y), quat.y, 0.01);
		assertEquals(Math.abs(orientation.z), quat.z, 0.01);
		assertEquals(Math.abs(orientation.w), quat.w, 0.01);
	}

	/**
	 * Test {@link Quaternions#orientationMetric(Point3d[], Point3d[])}.
	 * <p>
	 * Tests the range of values of the metric with a perfect correlation,
	 * perfect anticorrelation and intermediate values.
	 */
	@Test
	public void testOrientationMetric() {

		// no rotation quaternion
		Quat4d qa = new Quat4d(0, 0, 0, 1);
		Quat4d qb = new Quat4d(qa);
		
		assertEquals(Quaternions.orientationMetric(qa, qb), 0, 0.01);
		
		// rotate 180 degrees to get opposite orientation
		qb = new Quat4d(1, 0, 0, 0);
		
		assertEquals(Quaternions.orientationMetric(qa, qb), Math.PI / 2, 0.01);
		
		// rotate 90 degrees to get half-way orientation
		qb = new Quat4d(0.706, 0, 0, 0.706);
		
		assertEquals(Quaternions.orientationMetric(qa, qb), Math.PI / 4, 0.01);
	}
}
