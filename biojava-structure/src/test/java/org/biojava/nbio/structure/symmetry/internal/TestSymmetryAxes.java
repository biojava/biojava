package org.biojava.nbio.structure.symmetry.internal;

import static org.junit.Assert.*;

import java.util.Arrays;
import java.util.List;

import javax.vecmath.AxisAngle4d;
import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

import org.biojava.nbio.structure.symmetry.internal.CESymmParameters.SymmetryType;
import org.junit.Test;


public class TestSymmetryAxes {

	@Test
	public void testClosedCase() {
		// D4 case
		SymmetryAxes axes = new SymmetryAxes();

		// Level 1 is C4 along Z
		Matrix4d r90 = new Matrix4d();
		r90.set(new AxisAngle4d(0, 0, 1, -Math.PI/2));
		axes.addAxis(r90, 4, SymmetryType.CLOSED);
		
		// Level 2 is C2 along X
		Matrix4d r180 = new Matrix4d();
		r180.set(new AxisAngle4d(1, 0, 0, Math.PI));
		axes.addAxis(r180, 2, SymmetryType.CLOSED);

		assertEquals(2,axes.getElementaryAxes().size());

		Matrix4d expectedEven = new Matrix4d();
		expectedEven.set(1);
		Matrix4d expectedOdd = new Matrix4d(r180);
		assertEquals(expectedEven, axes.getRepeatTransform(0));
		assertEquals(expectedOdd, axes.getRepeatTransform(1));
		expectedEven.mul(r90);
		expectedOdd.mul(r180,expectedEven);
		assertEquals(expectedEven, axes.getRepeatTransform(2));
		assertEquals(expectedOdd, axes.getRepeatTransform(3));
		expectedEven.mul(r90);
		expectedOdd.mul(r180,expectedEven);
		assertEquals(expectedEven, axes.getRepeatTransform(4));
		assertEquals(expectedOdd, axes.getRepeatTransform(5));
		expectedEven.mul(r90);
		expectedOdd.mul(r180,expectedEven);
		assertEquals(expectedEven, axes.getRepeatTransform(6));
		assertEquals(expectedOdd, axes.getRepeatTransform(7));

		List<List<Integer>> relation = Arrays.asList(
				Arrays.asList(0,1,2,3,4,5,6,7),
				Arrays.asList(2,3,4,5,6,7,0,1)
				);
		assertEquals(relation,axes.getRepeatRelation(0));
		relation = Arrays.asList(
				Arrays.asList(0,1),
				Arrays.asList(1,0)
				);
		assertEquals(relation,axes.getRepeatRelation(1));
		try {
			axes.getRepeatRelation(2);
			fail("Invalid level");
		} catch(IndexOutOfBoundsException e) {}
		
		
		// Expected location of each repeat
		Point3d[] repeats = new Point3d[] {
				new Point3d(1,1,1),
				new Point3d(1,-1,-1),
				new Point3d(-1,1,1),
				new Point3d(1,1,-1),
				new Point3d(-1,-1,1),
				new Point3d(-1,1,-1),
				new Point3d(1,-1,1),
				new Point3d(-1,-1,-1)
		};
		// Inverse of (1,1,1) should give above points
		for(int i=0;i<8;i++) {
			Matrix4d m = axes.getRepeatTransform(i);
			m.invert();
			Point3d x = new Point3d(repeats[0]);
			m.transform(x);
			assertTrue("Transformation "+i+"^-1 of "+repeats[0]+ "="+x+" not "+repeats[i],x.epsilonEquals(repeats[i], 1e-5));
		}
		// Forward should map the above points onto the first one
		for(int i=0;i<8;i++) {
			Matrix4d m = axes.getRepeatTransform(i);
			Point3d x = new Point3d(repeats[i]);
			m.transform(x);
			assertTrue("Transformation "+i+" of "+repeats[i]+ "="+x+" not 1,1,1",x.epsilonEquals(repeats[0], 1e-5));
		}
		
		Point3d x;
		
		List<Matrix4d> symmetryAxes = axes.getSymmetryAxes();
		assertEquals(5,symmetryAxes.size());
		int axisNum = 0;
		// Repeat 2 -> 0 (90 deg around z)
		x = new Point3d(repeats[2]);
		symmetryAxes.get(axisNum).transform(x);
		assertTrue(String.format("SymmAxis %d of %s=%s not %s",axisNum,round(repeats[2]),round(x),round(repeats[0])),x.epsilonEquals(repeats[0], 1e-5));
		axisNum++;
		// Repeat 1 -> 0 (180 deg around x)
		x = new Point3d(repeats[1]);
		symmetryAxes.get(axisNum).transform(x);
		assertTrue(String.format("SymmAxis %d of %s=%s not %s",axisNum,round(repeats[1]),round(x),round(repeats[0])),x.epsilonEquals(repeats[0], 1e-5));
		axisNum++;
		// Repeat 3 -> 2 (180 deg around y)
		x = new Point3d(repeats[3]);
		symmetryAxes.get(axisNum).transform(x);
		assertTrue(String.format("SymmAxis %d of %s=%s not %s",axisNum,round(repeats[3]),round(x),round(repeats[2])),x.epsilonEquals(repeats[2], 1e-5));
		axisNum++;
		// Repeat 5 -> 4 (180 deg around x)
		x = new Point3d(repeats[5]);
		symmetryAxes.get(axisNum).transform(x);
		assertTrue(String.format("SymmAxis %d of %s=%s not %s",axisNum,round(repeats[5]),round(x),round(repeats[4])),x.epsilonEquals(repeats[4], 1e-5));
		axisNum++;
		// Repeat 7 -> 6 (180 deg around y)
		x = new Point3d(repeats[7]);
		symmetryAxes.get(axisNum).transform(x);
		assertTrue(String.format("SymmAxis %d of %s=%s not %s",axisNum,round(repeats[7]),round(x),round(repeats[6])),x.epsilonEquals(repeats[6], 1e-5));
		axisNum++;
	}
	private static Point3d round(Point3d p) {
		return new Point3d(Math.round(p.x*100)/100.,Math.round(p.y*100)/100.,Math.round(p.z*100)/100.);
	}

	@Test
	public void testOpenCase() {
		// D4 case
		SymmetryAxes axes = new SymmetryAxes();

		// Level 1 is R4 along X
		Matrix4d t10 = new Matrix4d();
		t10.set(1,new Vector3d(-10,0,0));
		axes.addAxis(t10, 4, SymmetryType.OPEN);

		// Level 2 is C2 along X
		Matrix4d r180 = new Matrix4d();
		r180.set(new AxisAngle4d(1, 0, 0, Math.PI));
		axes.addAxis(r180, 2, SymmetryType.CLOSED);

		assertEquals(2,axes.getElementaryAxes().size());

		Matrix4d expectedEven = new Matrix4d();
		expectedEven.set(1);
		Matrix4d expectedOdd = new Matrix4d(r180);
		assertEquals(expectedEven, axes.getRepeatTransform(0));
		assertEquals(expectedOdd, axes.getRepeatTransform(1));
		expectedEven.mul(t10);
		expectedOdd.mul(r180,expectedEven);
		assertEquals(expectedEven, axes.getRepeatTransform(2));
		assertEquals(expectedOdd, axes.getRepeatTransform(3));
		expectedEven.mul(t10);
		expectedOdd.mul(r180,expectedEven);
		assertEquals(expectedEven, axes.getRepeatTransform(4));
		assertEquals(expectedOdd, axes.getRepeatTransform(5));
		expectedEven.mul(t10);
		expectedOdd.mul(r180,expectedEven);
		assertEquals(expectedEven, axes.getRepeatTransform(6));
		assertEquals(expectedOdd, axes.getRepeatTransform(7));

		List<List<Integer>> relation = Arrays.asList(
				Arrays.asList(0,1,2,3,4,5),
				Arrays.asList(2,3,4,5,6,7)
				);
		assertEquals(relation,axes.getRepeatRelation(0));
		relation = Arrays.asList(
				Arrays.asList(0,1),
				Arrays.asList(1,0)
				);
		assertEquals(relation,axes.getRepeatRelation(1));
		try {
			axes.getRepeatRelation(2);
			fail("Invalid level");
		} catch(IndexOutOfBoundsException e) {}

		// Expected location of each repeat
		Point3d[] repeats = new Point3d[] {
				new Point3d(-15,1,1),
				new Point3d(-15,-1,-1),
				new Point3d( -5,1,1),
				new Point3d( -5,-1,-1),
				new Point3d(  5,1,1),
				new Point3d(  5,-1,-1),
				new Point3d( 15,1,1),
				new Point3d( 15,-1,-1),

		};
		// Inverse of first point should give above points
		for(int i=0;i<8;i++) {
			Matrix4d m = axes.getRepeatTransform(i);
			m.invert();
			Point3d x = new Point3d(repeats[0]);
			m.transform(x);
			assertTrue("Transformation "+i+"^-1 of "+repeats[0]+ "="+x+" not "+repeats[i],x.epsilonEquals(repeats[i], 1e-5));
		}
		// Forward should map the above points onto the first one
		for(int i=0;i<8;i++) {
			Matrix4d m = axes.getRepeatTransform(i);
			Point3d x = new Point3d(repeats[i]);
			m.transform(x);
			assertTrue("Transformation "+i+" of "+repeats[i]+ "="+x+" not "+repeats[0],x.epsilonEquals(repeats[0], 1e-5));
		}
		
		Point3d x;
		
		List<Matrix4d> symmetryAxes = axes.getSymmetryAxes();
		assertEquals(5,symmetryAxes.size());
		int axisNum = 0;
		// Repeat 2 -> 0 (shift 1)
		x = new Point3d(repeats[2]);
		symmetryAxes.get(axisNum).transform(x);
		assertTrue(String.format("SymmAxis %d of %s=%s not %s",axisNum,round(repeats[2]),round(x),round(repeats[0])),x.epsilonEquals(repeats[0], 1e-5));
		axisNum++;
		// All of these are actually equivalent
		// Repeat 1 -> 0 (180 deg around x)
		x = new Point3d(repeats[1]);
		symmetryAxes.get(axisNum).transform(x);
		assertTrue(String.format("SymmAxis %d of %s=%s not %s",axisNum,round(repeats[1]),round(x),round(repeats[0])),x.epsilonEquals(repeats[0], 1e-5));
		axisNum++;
		// Repeat 3 -> 2 (180 deg around x)
		x = new Point3d(repeats[3]);
		symmetryAxes.get(axisNum).transform(x);
		assertTrue(String.format("SymmAxis %d of %s=%s not %s",axisNum,round(repeats[3]),round(x),round(repeats[2])),x.epsilonEquals(repeats[2], 1e-5));
		axisNum++;
		// Repeat 5 -> 4 (180 deg around x)
		x = new Point3d(repeats[5]);
		symmetryAxes.get(axisNum).transform(x);
		assertTrue(String.format("SymmAxis %d of %s=%s not %s",axisNum,round(repeats[5]),round(x),round(repeats[4])),x.epsilonEquals(repeats[4], 1e-5));
		axisNum++;
		// Repeat 7 -> 6 (180 deg around x)
		x = new Point3d(repeats[7]);
		symmetryAxes.get(axisNum).transform(x);
		assertTrue(String.format("SymmAxis %d of %s=%s not %s",axisNum,round(repeats[7]),round(x),round(repeats[6])),x.epsilonEquals(repeats[6], 1e-5));
		axisNum++;
	}
	/**
	 * Test that the deprecated addAxis still works
	 */
	@SuppressWarnings("deprecation")
	@Test
	public void testOpenCaseOld() {
		// D4 case
		SymmetryAxes axes = new SymmetryAxes();

		// Level 1 is R4 along X
		Matrix4d t10 = new Matrix4d();
		t10.set(1,new Vector3d(1,0,0));
		List<Integer> repeats = Arrays.asList(0,0,1,1,2,2,3,3);
		List<List<Integer>> superposition = Arrays.asList(
				Arrays.asList(0,1,2,3,4,5),
				Arrays.asList(2,3,4,5,6,7));
		axes.addAxis(t10, superposition, repeats, 4);

		// Level 2 is C2 along X
		Matrix4d r180 = new Matrix4d();
		r180.set(new AxisAngle4d(1, 0, 0, Math.PI));
		repeats = Arrays.asList(0,1,0,1,0,1,0,1);
		superposition = Arrays.asList(
				Arrays.asList(0,1),
				Arrays.asList(1,0));
		axes.addAxis(r180, superposition, repeats, 2);


		assertEquals(2,axes.getElementaryAxes().size());

		Matrix4d expectedEven = new Matrix4d();
		expectedEven.set(1);
		Matrix4d expectedOdd = new Matrix4d(r180);
		assertEquals(expectedEven, axes.getRepeatTransform(0));
		assertEquals(expectedOdd, axes.getRepeatTransform(1));
		expectedEven.mul(t10);
		expectedOdd.mul(r180,expectedEven);
		assertEquals(expectedEven, axes.getRepeatTransform(2));
		assertEquals(expectedOdd, axes.getRepeatTransform(3));
		expectedEven.mul(t10);
		expectedOdd.mul(r180,expectedEven);
		assertEquals(expectedEven, axes.getRepeatTransform(4));
		assertEquals(expectedOdd, axes.getRepeatTransform(5));
		expectedEven.mul(t10);
		expectedOdd.mul(r180,expectedEven);
		assertEquals(expectedEven, axes.getRepeatTransform(6));
		assertEquals(expectedOdd, axes.getRepeatTransform(7));

		List<List<Integer>> relation = Arrays.asList(
				Arrays.asList(0,1,2,3,4,5),
				Arrays.asList(2,3,4,5,6,7)
				);
		assertEquals(relation,axes.getRepeatRelation(0));
		relation = Arrays.asList(
				Arrays.asList(0,1),
				Arrays.asList(1,0)
				);
		assertEquals(relation,axes.getRepeatRelation(1));
		try {
			axes.getRepeatRelation(2);
			fail("Invalid level");
		} catch(IndexOutOfBoundsException e) {}

	}
}
