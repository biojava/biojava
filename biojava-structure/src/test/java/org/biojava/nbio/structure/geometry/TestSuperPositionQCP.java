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

import java.util.Random;

import javax.vecmath.AxisAngle4d;
import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

import org.biojava.nbio.structure.geometry.SuperPositionQuat;
import org.biojava.nbio.structure.geometry.SuperPositionQCP;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Test the Quaternion-Based Characteristic Polynomial {@link SuperPositionQCP}
 * algorithm for RMSD and Superposition calculations.
 * 
 * @author Aleix Lafita
 * @since 5.0.0
 *
 */
public class TestSuperPositionQCP {

	private static final Logger LOGGER = LoggerFactory.getLogger(TestSuperPosition.class);

	/**
	 * Test case proposed by Peter Rose from his observations about quaternary
	 * symmetry artifacts with the QCP algorithm.
	 */
	@Test
	public void testSymmetryQCP() {

		// Generate an array of points with symmetry
		Point3d[] set1 = new Point3d[16];
		set1[0] = new Point3d(14.065934, 47.068832, -32.895836);
		set1[1] = new Point3d(-14.065934, -47.068832, -32.895836);
		set1[2] = new Point3d(-47.068832, 14.065934, -32.895836);
		set1[3] = new Point3d(47.068832, -14.065934, -32.895836);
		set1[4] = new Point3d(-14.065934, 47.068832, 32.895836);
		set1[5] = new Point3d(14.065934, -47.068832, 32.895836);
		set1[6] = new Point3d(47.068832, 14.065934, 32.895836);
		set1[7] = new Point3d(-47.068832, -14.065934, 32.895836);
		set1[8] = new Point3d(43.813946, 22.748293, -32.14434);
		set1[9] = new Point3d(-43.813946, -22.748293, -32.14434);
		set1[10] = new Point3d(-22.748293, 43.813946, -32.14434);
		set1[11] = new Point3d(22.748293, -43.813946, -32.14434);
		set1[12] = new Point3d(-43.813946, 22.748293, 32.14434);
		set1[13] = new Point3d(43.813946, -22.748293, 32.14434);
		set1[14] = new Point3d(22.748293, 43.813946, 32.14434);
		set1[15] = new Point3d(-22.748293, -43.813946, 32.14434);

		Point3d[] set2 = CalcPoint.clonePoint3dArray(set1);

		// Use a random transformation to set2
		AxisAngle4d rotAxis = new AxisAngle4d(0.440, 0.302, 0.845, 1.570);
		Vector3d translation = new Vector3d(0.345, 2.453, 5.324);
		Matrix4d transform = new Matrix4d();
		transform.set(rotAxis);
		transform.setTranslation(translation);
		CalcPoint.transform(transform, set2);

		// Use Quaternion superposition to obtain the RMSD
		SuperPosition algorithm = new SuperPositionQuat(false);
		long quatStart = System.nanoTime();
		double quatrmsd = algorithm.getRmsd(set1, set2);
		long quatTime = (System.nanoTime() - quatStart) / 1000;

		// Use QCP algorithm to get the RMSD
		algorithm = new SuperPositionQCP(false);
		long qcpStart = System.nanoTime();
		double qcprmsd = algorithm.getRmsd(set1, set2);
		long qcpTime = (System.nanoTime() - qcpStart) / 1000;

		LOGGER.info(String.format("RMSD Symmetry: Quat time: %d us" + ", QCP time: %d us", quatTime, qcpTime));

		// Check that the returned RMSDs are equal
		assertEquals(quatrmsd, qcprmsd, 0.001);

	}

	/**
	 * Test case proposed by Peter Rose to check the alternative use of QCP,
	 * where first the RMSD is checked before obtaining the transformation
	 * matrix, in order to speed up large-scale calculations.
	 */
	@Test
	public void testAlternativeUsageQCP() {

		// Transformation applied to cloud points 1 that needs to be recovered
		// by the superposition method
		AxisAngle4d rotAxis = new AxisAngle4d(0.440, 0.302, 0.845, 1.570);
		Vector3d translation = new Vector3d(0.345, 2.453, 5.324);

		Matrix4d transform = new Matrix4d();
		transform.set(rotAxis);
		transform.setTranslation(translation);

		// Generate a random artificial array of points
		Random rnd = new Random(0);

		transform = new Matrix4d();
		transform.set(rotAxis);
		transform.setTranslation(translation);

		Point3d[] c1 = new Point3d[500];
		Point3d[] c2 = new Point3d[500];

		for (int p = 0; p < 500; p++) {

			Point3d a = new Point3d(rnd.nextInt(100), rnd.nextInt(50), rnd.nextInt(150));
			c1[p] = a;

			// Add some noise to the second point
			Point3d b = new Point3d(a.x + rnd.nextDouble(), a.y + rnd.nextDouble(), a.z + rnd.nextDouble());
			c2[p] = b;

		}

		CalcPoint.transform(transform, c1);
		
		SuperPositionQCP qcp = new SuperPositionQCP(false);

		// Step 1 calculate RMSD
		long start = System.nanoTime() / 1000;
		qcp.getRmsd(c1, c2);
		long rmsdTime = (System.nanoTime() / 1000 - start);

		// Step 2 Obtain the matrix after RMSD
		Matrix4d trans1 = qcp.superposeAfterRmsd();
		long trans1time = (System.nanoTime() / 1000 - start) - rmsdTime;

		// Now obtain the matrix from scratch
		Matrix4d trans2 = qcp.superpose(c1, c2);
		long trans2time = (System.nanoTime() / 1000 - start) - trans1time;

		LOGGER.info(String.format(
				"Time for RMSD: %d us, superposition after RMSD: %d us, and superposition from scratch: %d us",
				rmsdTime, trans1time, trans2time));

		// Check the results are the same
		assertTrue(trans1.epsilonEquals(trans2, 0.05));

	}

}
