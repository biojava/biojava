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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

import javax.vecmath.AxisAngle4d;
import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.geometry.SuperPositionQuat;
import org.biojava.nbio.structure.geometry.SuperPositionQCP;
import org.junit.Before;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Test and compare the different superposition methods implemented in BioJava.
 * 
 * @author Aleix Lafita
 * @author Jose Duarte
 * @since 5.0.0
 *
 */
public class TestSuperPosition {

	private static final Logger LOGGER = LoggerFactory.getLogger(TestSuperPosition.class);

	private static List<Point3d[]> cloud1;
	private static List<Point3d[]> cloud2;

	// the transformation to apply to cloud points 1 that needs to be recovered by the superposition code
	private static final AxisAngle4d rotAxis = new AxisAngle4d(0.440, 0.302, 0.845, 1.570);
	private static final Vector3d translation = new Vector3d(0.345, 2.453, 5.324);;
	private static Matrix4d transform;
		
	// a translation to apply to cloud point 2 for the rmsd test only
	private static final Vector3d translation2 = new Vector3d(1.32, -1.03, 6.28);

	/**
	 * Generate two clouds of random points of different sizes to test
	 * correctness and performance of superposition algorithms.
	 * 
	 * @throws StructureException
	 */
	@Before
	public void setUp() throws StructureException {

		cloud1 = new ArrayList<Point3d[]>(5);
		cloud2 = new ArrayList<Point3d[]>(5);

		Random rnd = new Random(0);

		transform = new Matrix4d();
		transform.set(rotAxis);
		transform.setTranslation(translation);

		List<Integer> sizes = Arrays.asList(5, 50, 500, 5000, 50000, 500000);

		for (Integer size : sizes) {

			Point3d[] c1 = new Point3d[size];
			Point3d[] c2 = new Point3d[size];

			for (int p = 0; p < size; p++) {

				Point3d a = new Point3d(rnd.nextInt(100), rnd.nextInt(50),
						rnd.nextInt(150));
				c1[p] = a;

				// Add some noise
				Point3d b = new Point3d(a.x + rnd.nextDouble(), a.y
						+ rnd.nextDouble(), a.z + rnd.nextDouble());
				c2[p] = b;
			}

			CalcPoint.center(c1);
			CalcPoint.center(c2);

			CalcPoint.transform(transform, c1);			

			cloud1.add(c1);
			cloud2.add(c2);
					
			Point3d centroid1 = CalcPoint. centroid(c1);
			Point3d centroid2 = CalcPoint. centroid(c2);
			LOGGER.debug("Centroid c1 (size %d): (%.2f, %.2f, %.2f)\n", size, centroid1.x, centroid1.y, centroid1.z);
			LOGGER.debug("Centroid c2 (size %d): (%.2f, %.2f, %.2f)\n", size, centroid2.x, centroid2.y, centroid2.z);
		}

	}

	/**
	 * Test method to obtain the transformation matrix from superpositions.
	 */
	@Test
	public void testSuperposition() {

		for (int c = 0; c < cloud1.size(); c++) {

			// Use SVD superposition
			SuperPosition svd = new SuperPositionSVD(false);
			long svdStart = System.nanoTime();
			Matrix4d svdTransform = svd.superpose(cloud1.get(c), cloud2.get(c));
			long svdTime = (System.nanoTime() - svdStart) / 1000;

			// Use quaternion superposition
			SuperPosition quat = new SuperPositionQuat(false);
			long quatStart = System.nanoTime();
			Matrix4d quatTransform = quat.superpose(cloud1.get(c), cloud2.get(c));
			long quatTime = (System.nanoTime() - quatStart) / 1000;

			// Use QCP algorithm
			SuperPosition qcp = new SuperPositionQCP(false);
			long qcpStart = System.nanoTime();
			Matrix4d qcpTransform = qcp.superpose(cloud1.get(c), cloud2.get(c));
			long qcpTime = (System.nanoTime() - qcpStart) / 1000;

			LOGGER.info(String.format("Transformation Matrix %d points: "
					+ "SVD time %d us, SP time: %d us, QCP time: %d us",
					cloud1.get(c).length, svdTime, quatTime, qcpTime));

			// Check that the transformation matrix was recovered
			assertTrue(transform.epsilonEquals(svdTransform, 0.05));
			assertTrue(transform.epsilonEquals(quatTransform, 0.05));
			assertTrue(transform.epsilonEquals(qcpTransform, 0.05));
		}

	}

	/**
	 * Test method to obtain the RMSD of a superposition.
	 */
	@Test
	public void testRMSD() {
		
		// for the rmsd test we first make sure that both cloud points are not centered in origin so that the centering is tested too
		// first cloud points are already centered, we translate cloud2 only
		for (int c=0; c<cloud2.size(); c++) {
			CalcPoint.translate(translation2, cloud2.get(c));
			
			Point3d centroid2 = CalcPoint. centroid(cloud2.get(c));			
			LOGGER.debug("Centroid c2 (index %d): (%.2f, %.2f, %.2f)\n", c, centroid2.x, centroid2.y, centroid2.z);
		}
		

		for (int c = 0; c < cloud1.size(); c++) {

			// Use SVD superposition
			SuperPosition svd = new SuperPositionSVD(false);
			long svdStart = System.nanoTime();
			double svdrmsd = svd.getRmsd(cloud1.get(c), cloud2.get(c));
			long svdTime = (System.nanoTime() - svdStart) / 1000;

			// Use quaternion superposition
			SuperPosition quat = new SuperPositionQuat(false);
			long quatStart = System.nanoTime();
			double quatrmsd = quat.getRmsd(cloud1.get(c), cloud2.get(c));
			long quatTime = (System.nanoTime() - quatStart) / 1000;

			// Use QCP algorithm
			SuperPosition qcp = new SuperPositionQCP(false);
			long qcpStart = System.nanoTime();
			double qcprmsd = qcp.getRmsd(cloud1.get(c), cloud2.get(c));
			long qcpTime = (System.nanoTime() - qcpStart) / 1000;

			LOGGER.info(String.format("RMSD %d points: SVD time %d us, "
					+ "Quat time: %d us, QCP time: %d us", cloud1.get(c).length,
					svdTime, quatTime, qcpTime));

			// Check that the returned RMSDs are equal
			assertEquals(quatrmsd, qcprmsd, 0.001);
			assertEquals(svdrmsd, quatrmsd, 0.001);
			assertEquals(svdrmsd, qcprmsd, 0.001);
		}
	}

}
