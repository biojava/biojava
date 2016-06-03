package org.biojava.nbio.structure.symmetry.internal;

import static org.junit.Assert.*;

import java.util.Arrays;
import java.util.List;

import javax.vecmath.AxisAngle4d;
import javax.vecmath.Matrix4d;
import javax.vecmath.Vector3d;

import org.junit.Test;


public class TestSymmetryAxes {

	@Test
	public void testClosedCase() {
		// D4 case
		SymmetryAxes axes = new SymmetryAxes();
		
		// Level 1 is C4 along Z
		Matrix4d r90 = new Matrix4d();
		r90.set(new AxisAngle4d(0, 0, 1, Math.PI/2));
		// number of times to apply this op for each repeat
		List<Integer> repeats = Arrays.asList(0,0,1,1,2,2,3,3);
		// aligned repeats
		List<List<Integer>> superposition = Arrays.asList(
				Arrays.asList(0,1,2,3,4,5,6,7),
				Arrays.asList(2,3,4,5,6,7,0,1));
		axes.addAxis(r90, superposition, repeats, 4);
		
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
		assertNull(axes.getRepeatRelation(2));

	}
	
	@Test
	public void testOpenCase() {
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
		assertNull(axes.getRepeatRelation(2));

	}
}
