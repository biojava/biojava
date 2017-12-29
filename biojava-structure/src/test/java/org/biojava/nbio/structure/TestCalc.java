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
import javax.vecmath.Vector3d;

import org.biojava.nbio.structure.geometry.Matrices;
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
		Point3d actual = atom.getCoordsAsPoint3d();

		assertEquals(expected, actual);

		//Sample transform: calc transposes automatically the matrix
		//because it is a pre-multiplication rotation matrix
		Matrix sampleR = Matrices.getRotationJAMA(getSampleTransform());
		Atom sampleT = Calc.getTranslationVector(getSampleTransform());
		Calc.rotate(atom, sampleR);
		Calc.shift(atom, sampleT);

		expected = new Point3d(2.0, 7.0, -1.3);
		actual = atom.getCoordsAsPoint3d();

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
		Point3d actual = atom.getCoordsAsPoint3d();

		assertEquals(expected, actual);

		//Sample transform
		Matrix4d sample = getSampleTransform();
		Calc.transform(atom, sample);

		expected = new Point3d(2.0, 7.0, -1.3);
		actual = atom.getCoordsAsPoint3d();

		assertEquals(expected, actual);
	}
	
	/**
	 * Issue https://github.com/biojava/biojava/issues/715
	 */
	@Test
	public void testChainTransform() {
		
		Chain c = createDummyChain();
		
		Matrix4d m = new Matrix4d(1,0,0,1, 0,1,0,0, 0,0,1,0, 0,0,0,1); // shift of 1 in x axis
 		Calc.transform(c, m);
 		
 		Group thegroup = c.getAtomGroup(0);
 		Group thealtlocgroup = thegroup.getAltLocs().get(0);
 		
 		Atom atom1 = thegroup.getAtom("CA");
 		Atom atom2 = thealtlocgroup.getAtom("CA");
 		
 		// x should be shifted by 1
 		assertEquals(2, atom1.getX(), 0.00001);
 		assertEquals(1, atom1.getY(), 0.00001);
 		assertEquals(1, atom1.getZ(), 0.00001);
 		
 		// x should be shifted by 1
 		assertEquals(3, atom2.getX(), 0.00001);
 		assertEquals(2, atom2.getY(), 0.00001);
 		assertEquals(2, atom2.getZ(), 0.00001);

 		
	}

	/**
	 * Issue https://github.com/biojava/biojava/issues/715
	 */
	@Test
	public void testStructureTransform() {
		
		Structure s = createDummyStructure();
		
		Matrix4d m = new Matrix4d(1,0,0,1, 0,1,0,0, 0,0,1,0, 0,0,0,1); // shift of 1 in x axis
 		Calc.transform(s, m);
 		
 		// testing 1st chain
 		Group thegroup = s.getChain("A").getAtomGroup(0);
 		Group thealtlocgroup = thegroup.getAltLocs().get(0);
 		
 		Atom atom1 = thegroup.getAtom("CA");
 		Atom atom2 = thealtlocgroup.getAtom("CA");
 		
 		// x should be shitfted by 1
 		assertEquals(2, atom1.getX(), 0.00001);
 		assertEquals(1, atom1.getY(), 0.00001);
 		assertEquals(1, atom1.getZ(), 0.00001);
 		
 		// x should be shitfted by 1
 		assertEquals(3, atom2.getX(), 0.00001);
 		assertEquals(2, atom2.getY(), 0.00001);
 		assertEquals(2, atom2.getZ(), 0.00001);

 		// testing 2nd chain
 		thegroup = s.getChain("B").getAtomGroup(0);
 		thealtlocgroup = thegroup.getAltLocs().get(0);
 		
 		atom1 = thegroup.getAtom("CA");
 		atom2 = thealtlocgroup.getAtom("CA");
 		
 		// x should be shitfted by 1
 		assertEquals(4, atom1.getX(), 0.00001);
 		assertEquals(3, atom1.getY(), 0.00001);
 		assertEquals(3, atom1.getZ(), 0.00001);
 		
 		// x should be shitfted by 1
 		assertEquals(5, atom2.getX(), 0.00001);
 		assertEquals(4, atom2.getY(), 0.00001);
 		assertEquals(4, atom2.getZ(), 0.00001);

 		
	}
	
	@Test
	public void testChainTranslate() {
		Chain c = createDummyChain();

		Vector3d translation = new Vector3d(1, 0, 0);
 		Calc.translate(c, translation);
 		
 		Group thegroup = c.getAtomGroup(0);
 		Group thealtlocgroup = thegroup.getAltLocs().get(0);
 		
 		Atom atom1 = thegroup.getAtom("CA");
 		Atom atom2 = thealtlocgroup.getAtom("CA");
 		
 		// x should be shifted by 1
 		assertEquals(2, atom1.getX(), 0.00001);
 		assertEquals(1, atom1.getY(), 0.00001);
 		assertEquals(1, atom1.getZ(), 0.00001);
 		
 		// x should be shifted by 1
 		assertEquals(3, atom2.getX(), 0.00001);
 		assertEquals(2, atom2.getY(), 0.00001);
 		assertEquals(2, atom2.getZ(), 0.00001);
	}
	
	@Test
	public void testStructureTranslate() {
		Structure s = createDummyStructure();
		
		Vector3d translation = new Vector3d(1, 0, 0);
 		Calc.translate(s, translation);
 		
 		// testing 1st chain
 		Group thegroup = s.getChain("A").getAtomGroup(0);
 		Group thealtlocgroup = thegroup.getAltLocs().get(0);
 		
 		Atom atom1 = thegroup.getAtom("CA");
 		Atom atom2 = thealtlocgroup.getAtom("CA");
 		
 		// x should be shitfted by 1
 		assertEquals(2, atom1.getX(), 0.00001);
 		assertEquals(1, atom1.getY(), 0.00001);
 		assertEquals(1, atom1.getZ(), 0.00001);
 		
 		// x should be shitfted by 1
 		assertEquals(3, atom2.getX(), 0.00001);
 		assertEquals(2, atom2.getY(), 0.00001);
 		assertEquals(2, atom2.getZ(), 0.00001);

 		// testing 2nd chain
 		thegroup = s.getChain("B").getAtomGroup(0);
 		thealtlocgroup = thegroup.getAltLocs().get(0);
 		
 		atom1 = thegroup.getAtom("CA");
 		atom2 = thealtlocgroup.getAtom("CA");
 		
 		// x should be shitfted by 1
 		assertEquals(4, atom1.getX(), 0.00001);
 		assertEquals(3, atom1.getY(), 0.00001);
 		assertEquals(3, atom1.getZ(), 0.00001);
 		
 		// x should be shitfted by 1
 		assertEquals(5, atom2.getX(), 0.00001);
 		assertEquals(4, atom2.getY(), 0.00001);
 		assertEquals(4, atom2.getZ(), 0.00001);
	}

	private static Atom getAtom(String name, double x, double y, double z) {
		Atom a = new AtomImpl();
		a.setX(x);
		a.setY(y);
		a.setZ(z);
		a.setName(name);
		return a;
	}
	
	private static Atom getAtom(double x, double y, double z) {
		return getAtom(null, x, y, z);
	}

	private static Matrix4d getSampleTransform(){

		Matrix4d sample = new Matrix4d(new double[] {1.0,1.5,0.5,-1.0,
		                                0.5,0.5,1.0,5.0,
		                                0.3,0.4,0.5,-2.5,
		                                0.0,0.0,0.0,1.0});
		return sample;
	}
	
	private static Chain createDummyChain() {
		Group g = new AminoAcidImpl();
		Atom a = getAtom("CA", 1, 1, 1);
		g.addAtom(a);
		Group altLocG = new AminoAcidImpl();
		Atom a2 = getAtom("CA", 2, 2, 2);
		altLocG.addAtom(a2);
		
		g.addAltLoc(altLocG);
		
		Chain c = new ChainImpl();
		c.addGroup(g);
		return c;
	}
	
	private static Structure createDummyStructure() {
		Group g = new AminoAcidImpl();
		Atom a = getAtom("CA", 1, 1, 1);
		g.addAtom(a);
		Group altLocG = new AminoAcidImpl();
		Atom a2 = getAtom("CA", 2, 2, 2);
		altLocG.addAtom(a2);
		
		g.addAltLoc(altLocG);
		
		Chain c1 = new ChainImpl();
		c1.addGroup(g);
		c1.setId("A");
		
		Group gc2 = new AminoAcidImpl();
		Atom ac2 = getAtom("CA", 3, 3, 3);
		gc2.addAtom(ac2);
		Group altLocGc2 = new AminoAcidImpl();
		Atom ac22 = getAtom("CA", 4, 4, 4);
		altLocGc2.addAtom(ac22);
		
		gc2.addAltLoc(altLocGc2);
		
		Chain c2 = new ChainImpl();
		c2.addGroup(gc2);
		c2.setId("B");
		
		Structure s = new StructureImpl();
		s.addChain(c1);
		s.addChain(c2);
		return s;
	}

}
