package org.biojava.nbio.structure;

import static org.junit.Assert.*;

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
	
	
	
	private static Atom getAtom(double x, double y, double z) {
		Atom a = new AtomImpl();
		a.setX(x);
		a.setY(y);
		a.setZ(z);
		return a;
	}

	
}
