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
package org.biojava.nbio.structure;

import static org.junit.Assert.*;

import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;

import org.biojava.nbio.structure.jama.Matrix;
import org.junit.Test;

public class TestCalc {

	@Test
	public void testAngle() {
		Atom a = getAtom(1,1,1);
		
		Atom b = getAtom(2,2,2);
		
		// 2 parallel vectors should return 0
		assertEquals(0.0, Calc.angle(a, b),0.00001);
		
		// range should be [0,180]
		Atom ref = getAtom(1,0,0);
		Atom a1 = getAtom(1,1,0);
		Atom a2 = getAtom(0,1,0);
		Atom a3 = getAtom(-1,1,0);
		Atom a4 = getAtom(-1,0,0);
		Atom a5 = getAtom(-1,-1,0);
		Atom a6 = getAtom(0,-1,0);
		Atom a7 = getAtom(1,-1,0);
		
		assertEquals(180.0,Calc.angle(ref, a4),0.00001);
		
		assertTrue(Calc.angle(ref, a1)>=0 && Calc.angle(ref,a1)<=180.0);
		assertTrue(Calc.angle(ref, a2)>=0 && Calc.angle(ref,a2)<=180.0);
		assertTrue(Calc.angle(ref, a3)>=0 && Calc.angle(ref,a3)<=180.0);
		assertTrue(Calc.angle(ref, a4)>=0 && Calc.angle(ref,a4)<=180.0);
		assertTrue(Calc.angle(ref, a5)>=0 && Calc.angle(ref,a5)<=180.0);
		assertTrue(Calc.angle(ref, a6)>=0 && Calc.angle(ref,a6)<=180.0);
		assertTrue(Calc.angle(ref, a7)>=0 && Calc.angle(ref,a7)<=180.0);
		
		Atom c = getAtom(0,0,0);
		Atom d = getAtom(0,0,0);
		
		assertEquals(Double.NaN, Calc.angle(a,c),0.00001);
		assertEquals(Double.NaN, Calc.angle(c,d),0.00001);
	}
	
	@Test
	public void testTorsionAngle() {
		Atom a = getAtom(0,0,0);
		
		Atom b = getAtom(1,0,0);
		
		Atom c = getAtom(2,0,0);
		
		Atom d = getAtom(3,0,0);
		
		// all 4 points colinear
		
		assertEquals(Double.NaN, Calc.torsionAngle(a, b, c, d),0.00001);
		
		// first 3 colinear
		d = getAtom(3,1,0);
		assertEquals(Double.NaN, Calc.torsionAngle(a, b, c, d),0.00001);
		
		// second 3 colinear
		d = getAtom(3,0,0);
		a = getAtom(1,1,0);
		assertEquals(Double.NaN, Calc.torsionAngle(a, b, c, d),0.00001);
		
		// coplanar vectors
		a = getAtom(0,0,0);
		b = getAtom(1,0,0);
		c = getAtom(1, 1, 0);
		d = getAtom(2, 1, 0);
		
		assertEquals(180, Calc.torsionAngle(a, b, c, d),0.00001);
		
		c = getAtom(-1,-1, 0);
		d = getAtom(-2,-1, 0);
		
		assertEquals(0, Calc.torsionAngle(a, b, c, d),0.00001);
	}
	
	@Test
	public void testJamaTransformation() {
		
		Atom atom = getAtom(1.0, 1.0, 1.0);
		
		//Identity transform
		Matrix identR = Matrix.identity(3, 3);
		Atom identT = getAtom(0, 0, 0);
		Calc.rotate(atom, identR);
		Calc.shift(atom, identT);
		
		Point3d expected = new Point3d(1.0, 1.0, 1.0);
		Point3d actual = new Point3d(atom.getCoords());
		
		assertEquals(expected, actual);
		
		//Sample transform: calc transposes automatically the matrix
		//because it is a pre-multiplication rotation matrix
		Matrix sampleR = Calc.getRotationMatrix(getSampleTransform());
		Atom sampleT = Calc.getTranslationVector(getSampleTransform());
		Calc.rotate(atom, sampleR);
		Calc.shift(atom, sampleT);
		
		expected = new Point3d(2.0, 7.0, -1.3);
		actual = new Point3d(atom.getCoords());
		
		assertEquals(expected, actual);
	}
	
	@Test
	public void testVecmathTransformation() {
		
		Atom atom = getAtom(1.0, 1.0, 1.0);
		
		//Identity transform
		Matrix4d ident = new Matrix4d();
		ident.setIdentity();
		Calc.transform(atom, ident);
		
		Point3d expected = new Point3d(1.0, 1.0, 1.0);
		Point3d actual = new Point3d(atom.getCoords());
		
		assertEquals(expected, actual);
		
		//Sample transform
		Matrix4d sample = getSampleTransform();
		Calc.transform(atom, sample);
		
		expected = new Point3d(2.0, 7.0, -1.3);
		actual = new Point3d(atom.getCoords());
		
		assertEquals(expected, actual);
	}
	
	private static Atom getAtom(double x, double y, double z) {
		Atom a = new AtomImpl();
		a.setX(x);
		a.setY(y);
		a.setZ(z);
		return a;
	}
	
	private static Matrix4d getSampleTransform(){
		
		Matrix4d sample = new Matrix4d(new double[] {1.0,1.5,0.5,-1.0,
		                                0.5,0.5,1.0,5.0,
		                                0.3,0.4,0.5,-2.5,
		                                0.0,0.0,0.0,1.0});
		return sample;
	}

}
