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

import static org.junit.Assert.*;

import java.io.IOException;

import javax.vecmath.AxisAngle4d;
import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;
import javax.vecmath.Quat4d;
import javax.vecmath.Vector3d;

import org.biojava.nbio.structure.Calc;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.geometry.UnitQuaternions;
import org.junit.Test;

/**
 * Test the methods in the {@link UnitQuaternions} class.
 * 
 * @author Aleix Lafita
 * @since 5.0.0
 *
 */
public class TestUnitQuaternions {

	/**
	 * Test {@link UnitQuaternions#orientation(javax.vecmath.Point3d[])}.
	 * <p>
	 * Tests the identity orientation, orientation around one coordinate axis
	 * and orientation around a non-coordinate axis.
	 * 
	 * @throws StructureException
	 * @throws IOException
	 */
	@Test
	public void testOrientation() throws IOException, StructureException {

		// Get points from a structure. It is difficult to generate points
		// with no bias in their distribution (too uniform, ie).
		Structure pdb = StructureIO.getStructure("4hhb.A");
		Point3d[] cloud = Calc.atomsToPoints(StructureTools
				.getRepresentativeAtomArray(pdb));

		// Center the cloud at the origin
		CalcPoint.center(cloud);

		// Orient its principal axes to the coordinate axis
		Quat4d orientation = UnitQuaternions.orientation(cloud);
		Matrix4d transform = new Matrix4d();
		transform.set(orientation);
		transform.invert();
		CalcPoint.transform(transform, cloud);

		// The orientation found now should be 0 (it has been re-oriented)
		orientation = UnitQuaternions.orientation(cloud);
		AxisAngle4d axis = new AxisAngle4d();
		axis.set(orientation);

		// No significant rotation
		assertEquals(orientation.x, 0.0, 0.01);
		assertEquals(orientation.y, 0.0, 0.01);
		assertEquals(orientation.z, 0.0, 0.01);
		assertEquals(axis.angle, 0.0, 0.01);

		// Now try to recover an orientation
		Quat4d quat = new Quat4d(0.418, 0.606, 0.303, 0.606);

		Matrix4d mat = new Matrix4d();
		mat.set(quat);

		CalcPoint.transform(mat, cloud);

		orientation = UnitQuaternions.orientation(cloud);

		// Test recovering the quaternion (q and -q same rotation)
		assertEquals(Math.abs(orientation.x), quat.x, 0.01);
		assertEquals(Math.abs(orientation.y), quat.y, 0.01);
		assertEquals(Math.abs(orientation.z), quat.z, 0.01);
		assertEquals(Math.abs(orientation.w), quat.w, 0.01);
	}

	/**
	 * Test {@link UnitQuaternions#orientationMetric(Point3d[], Point3d[])}.
	 * <p>
	 * Tests the range of values of the metric with a perfect correlation,
	 * perfect anticorrelation and intermediate values.
	 */
	@Test
	public void testOrientationMetricRange() {

		// no rotation quaternion
		Quat4d qa = new Quat4d(0, 0, 0, 1);
		Quat4d qb = new Quat4d(qa);

		// Two equal quaternions produce the minimum score of 0
		assertEquals(UnitQuaternions.orientationMetric(qa, qb), 0, 0.01);

		// 90 degrees rotation over x
		qa = new Quat4d(0.707, 0, 0, 0.707);

		// 270 degrees rotation over x
		qb = new Quat4d(0.707, 0, 0, -0.707);

		// two quaternions with 180 degree axis produce the max score Pi / 2
		assertEquals(UnitQuaternions.orientationMetric(qa, qb), Math.PI / 2,
				0.01);

		// 90 degrees rotation over y
		qb = new Quat4d(0, 0.707, 0, 0.707);

		// two quaternions with 90 degree axis produce the score Pi / 4
		assertEquals(UnitQuaternions.orientationMetric(qa, qb), Math.PI / 3,
				0.01);

		// two quaternions with 45 degree axis produce the score Pi / 8
		qb = new Quat4d(0.383, 0, 0, 0.924);

		assertEquals(UnitQuaternions.orientationMetric(qa, qb), Math.PI / 8,
				0.01);

		// 90 degrees rotation over x in negative
		qb = new Quat4d(0, -0.707, 0, -0.707);

		// assert no negative angles are returned
		assertEquals(UnitQuaternions.orientationMetric(qa, qb), Math.PI / 3,
				0.01);

	}

	/**
	 * Test {@link UnitQuaternions#orientationMetric(Point3d[], Point3d[])} on a
	 * real structure, which will be deviating a little bit every time.
	 * 
	 * @throws StructureException
	 * @throws IOException
	 */
	@Test
	public void testOrientationMetricIncrement() throws IOException,
			StructureException {

		// The rotation increment will be Pi/10, Pi/15 and Pi/12 degrees in X,Y
		// and Z
		Matrix4d transform = new Matrix4d();
		transform.rotX(Math.PI / 10);
		transform.rotY(Math.PI / 12);
		transform.rotZ(Math.PI / 15);

		// Get points from a structure.
		Structure pdb = StructureIO.getStructure("4hhb.A");
		Point3d[] cloud = Calc.atomsToPoints(StructureTools
				.getRepresentativeAtomArray(pdb));
		Point3d[] cloud2 = CalcPoint.clonePoint3dArray(cloud);

		// Center the clouds at the origin
		CalcPoint.center(cloud);
		CalcPoint.center(cloud2);

		// Their orientation is equal at this stage
		double m0 = UnitQuaternions.orientationMetric(cloud, cloud2);
		assertEquals(m0, 0.0, 0.001);

		// Assert it keeps incrementing every time transform is applied
		CalcPoint.transform(transform, cloud2);
		double m1 = UnitQuaternions.orientationMetric(cloud, cloud2);
		assertTrue(m1 > m0);

		CalcPoint.transform(transform, cloud2);
		double m2 = UnitQuaternions.orientationMetric(cloud, cloud2);
		assertTrue(m2 > m1);

		CalcPoint.transform(transform, cloud2);
		double m3 = UnitQuaternions.orientationMetric(cloud, cloud2);
		assertTrue(m3 > m2);

		CalcPoint.transform(transform, cloud2);
		double m4 = UnitQuaternions.orientationMetric(cloud, cloud2);
		assertTrue(m4 > m3);

		CalcPoint.transform(transform, cloud2);
		double m5 = UnitQuaternions.orientationMetric(cloud, cloud2);
		assertTrue(m5 > m4);

		CalcPoint.transform(transform, cloud2);
		double m6 = UnitQuaternions.orientationMetric(cloud, cloud2);
		assertTrue(m6 > m5);

		CalcPoint.transform(transform, cloud2);
		double m7 = UnitQuaternions.orientationMetric(cloud, cloud2);
		assertTrue(m7 > m6);

		CalcPoint.transform(transform, cloud2);
		double m8 = UnitQuaternions.orientationMetric(cloud, cloud2);
		assertTrue(m8 > m7);

		CalcPoint.transform(transform, cloud2);
		double m9 = UnitQuaternions.orientationMetric(cloud, cloud2);
		assertTrue(m9 > m8);
	}

	/**
	 * Test {@link UnitQuaternions#relativeOrientation(Point3d[], Point3d[])} on
	 * a real structure. Test recovering of the angle applied.
	 * 
	 * @throws StructureException
	 * @throws IOException
	 */
	@Test
	public void testRelativeOrientation() throws IOException,
			StructureException {

		// Get points from a structure.
		Structure pdb = StructureIO.getStructure("4hhb.A");
		Point3d[] cloud = Calc.atomsToPoints(StructureTools
				.getRepresentativeAtomArray(pdb));
		Point3d[] cloud2 = CalcPoint.clonePoint3dArray(cloud);
		
		// Test orientation angle equal to 0 at this point
		double angle = UnitQuaternions.orientationAngle(cloud, cloud2, false);
		assertEquals(angle, 0, 0.001);
		
		// Apply a 30 degree rotation to cloud
		AxisAngle4d axis = new AxisAngle4d(new Vector3d(1,1,1), Math.PI / 6);
		Matrix4d transform = new Matrix4d();
		transform.set(axis);
		
		CalcPoint.transform(transform, cloud);
		angle = UnitQuaternions.orientationAngle(cloud, cloud2, false);
		angle = Math.min(Math.abs(2 * Math.PI - angle), angle);
		
		// Test that angle was recovered
		assertEquals(angle, Math.PI / 6, 0.001);
	}

}
