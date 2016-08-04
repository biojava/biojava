package org.biojava.nbio.structure.geometry;

import static org.junit.Assert.*;

import java.io.IOException;

import javax.vecmath.AxisAngle4d;
import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;
import javax.vecmath.Quat4d;

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
	public void testOrientationMetric() {

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
	
	
}
