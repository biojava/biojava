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
package org.biojava.bio.structure.xtal;

import org.junit.Before;
import org.junit.Test;

import javax.vecmath.Point3d;
import javax.vecmath.Point3i;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

public class TestCrystalCell {

	@Before
	public void setUp() throws Exception {
	}

	@Test
	public void testGetCellIndices() {
		CrystalCell cell = new CrystalCell(100, 100, 100, 90, 90, 45);

		Point3i result,expected;
		Point3d query;

		query = new Point3d(0,0,0);
		expected = new Point3i(0,0,0);
		result = cell.getCellIndices(query);
		assertEquals("Wrong index for "+query,expected,result);

		query = new Point3d(99.9,0,0);
		expected = new Point3i(0,0,0);
		result = cell.getCellIndices(query);
		assertEquals("Wrong index for "+query,expected,result);

		query = new Point3d(100,0,0);
		expected = new Point3i(1,0,0);
		result = cell.getCellIndices(query);
		assertEquals("Wrong index for "+query,expected,result);

		query = new Point3d(0,50,0);
		expected = new Point3i(-1,0,0);
		result = cell.getCellIndices(query);
		assertEquals("Wrong index for "+query,expected,result);

		query = new Point3d(51,50,0);
		expected = new Point3i(0,0,0);
		result = cell.getCellIndices(query);
		assertEquals("Wrong index for "+query,expected,result);

		query = new Point3d(72,71,0);
		expected = new Point3i(0,1,0);
		result = cell.getCellIndices(query);
		assertEquals("Wrong index for "+query,expected,result);

		query = new Point3d(500,0,0);
		expected = new Point3i(5,0,0);
		result = cell.getCellIndices(query);
		assertEquals("Wrong index for "+query,expected,result);

		query = new Point3d(-500,0,0);
		expected = new Point3i(-5,0,0);
		result = cell.getCellIndices(query);
		assertEquals("Wrong index for "+query,expected,result);

		query = new Point3d(-550,0,0);
		expected = new Point3i(-6,0,0);
		result = cell.getCellIndices(query);
		assertEquals("Wrong index for "+query,expected,result);

		query = new Point3d(2,1,500);
		expected = new Point3i(0,0,5);
		result = cell.getCellIndices(query);
		assertEquals("Wrong index for "+query,expected,result);
	}

	@Test
	public void testTransformToOrigin() {
		CrystalCell cell = new CrystalCell(100, 100, 100, 90, 90, 45);
		double h = 100/Math.sqrt(2);
		double tol = 1e-6;

		Point3d query,expected;


		// Note that the 0 boundaries are instable
		// If tests break, it's ok to move into the cell (e.g. <2,1,1> is unambiguous)
		query = new Point3d(0,0,0);
		expected = new Point3d(0,0,0);
		cell.transfToOriginCell(query);
		assertTrue("Error transforming to origin. Expected:"+expected+" but was:"+query, expected.epsilonEquals(query, tol));

		query = new Point3d(99.9,0,0);
		expected = new Point3d(99.9,0,0);
		cell.transfToOriginCell(query);
		assertTrue("Error transforming to origin. Expected:"+expected+" but was:"+query, expected.epsilonEquals(query, tol));

		query = new Point3d(100,0,0);
		expected = new Point3d(0,0,0);
		cell.transfToOriginCell(query);
		assertTrue("Error transforming to origin. Expected:"+expected+" but was:"+query, expected.epsilonEquals(query, tol));

		query = new Point3d(0,50,0);
		expected = new Point3d(100,50,0);
		cell.transfToOriginCell(query);
		assertTrue("Error transforming to origin. Expected:"+expected+" but was:"+query, expected.epsilonEquals(query, tol));

		query = new Point3d(51,50,0);
		expected = new Point3d(51,50,0);
		cell.transfToOriginCell(query);
		assertTrue("Error transforming to origin. Expected:"+expected+" but was:"+query, expected.epsilonEquals(query, tol));

		query = new Point3d(h+2,h+1,0);
		expected = new Point3d(2,1,0);
		cell.transfToOriginCell(query);
		assertTrue("Error transforming to origin. Expected:"+expected+" but was:"+query, expected.epsilonEquals(query, tol));

		query = new Point3d(500,0,0);
		expected = new Point3d(0,0,0);
		cell.transfToOriginCell(query);
		assertTrue("Error transforming to origin. Expected:"+expected+" but was:"+query, expected.epsilonEquals(query, tol));

		query = new Point3d(-500,0,0);
		expected = new Point3d(0,0,0);
		cell.transfToOriginCell(query);
		assertTrue("Error transforming to origin. Expected:"+expected+" but was:"+query, expected.epsilonEquals(query, tol));

		query = new Point3d(2,1,500);
		expected = new Point3d(2,1,0);
		cell.transfToOriginCell(query);
		assertTrue("Error transforming to origin. Expected:"+expected+" but was:"+query, expected.epsilonEquals(query, tol));
	}

	@Test
	public void testTransformToOriginArray() {
		CrystalCell cell = new CrystalCell(100, 100, 100, 90, 90, 45);
		double h = 100/Math.sqrt(2);
		double tol = 1e-6;

		Point3d query,expected;


		// Note that the 0 boundaries are instable
		// If tests break, it's ok to move into the cell (e.g. <2,1,1> is unambiguous)
		query = new Point3d(0,0,0);
		expected = new Point3d(0,0,0);
		cell.transfToOriginCell(new Point3d[] {query}, query);
		assertTrue("Error transforming to origin. Expected:"+expected+" but was:"+query, expected.epsilonEquals(query, tol));

		query = new Point3d(99.9,0,0);
		expected = new Point3d(99.9,0,0);
		cell.transfToOriginCell(new Point3d[] {query}, query);
		assertTrue("Error transforming to origin. Expected:"+expected+" but was:"+query, expected.epsilonEquals(query, tol));

		query = new Point3d(100,0,0);
		expected = new Point3d(0,0,0);
		cell.transfToOriginCell(new Point3d[] {query}, query);
		assertTrue("Error transforming to origin. Expected:"+expected+" but was:"+query, expected.epsilonEquals(query, tol));

		query = new Point3d(0,50,0);
		expected = new Point3d(100,50,0);
		cell.transfToOriginCell(new Point3d[] {query}, query);
		assertTrue("Error transforming to origin. Expected:"+expected+" but was:"+query, expected.epsilonEquals(query, tol));

		query = new Point3d(51,50,0);
		expected = new Point3d(51,50,0);
		cell.transfToOriginCell(new Point3d[] {query}, query);
		assertTrue("Error transforming to origin. Expected:"+expected+" but was:"+query, expected.epsilonEquals(query, tol));

		query = new Point3d(h+2,h+1,0);
		expected = new Point3d(2,1,0);
		cell.transfToOriginCell(new Point3d[] {query}, query);
		assertTrue("Error transforming to origin. Expected:"+expected+" but was:"+query, expected.epsilonEquals(query, tol));

		query = new Point3d(500,0,0);
		expected = new Point3d(0,0,0);
		cell.transfToOriginCell(new Point3d[] {query}, query);
		assertTrue("Error transforming to origin. Expected:"+expected+" but was:"+query, expected.epsilonEquals(query, tol));

		query = new Point3d(-500,0,0);
		expected = new Point3d(0,0,0);
		cell.transfToOriginCell(new Point3d[] {query}, query);
		assertTrue("Error transforming to origin. Expected:"+expected+" but was:"+query, expected.epsilonEquals(query, tol));

		query = new Point3d(2,1,500);
		expected = new Point3d(2,1,0);
		cell.transfToOriginCell(new Point3d[] {query}, query);
		assertTrue("Error transforming to origin. Expected:"+expected+" but was:"+query, expected.epsilonEquals(query, tol));
	}

}
