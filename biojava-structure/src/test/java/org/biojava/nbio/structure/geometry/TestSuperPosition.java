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
 * Test the Quaternion-Based Characteristic Polynomial {@link SuperPositionQCP}
 * algorithm for RMSD and Superposition calculations.
 * 
 * @author Aleix Lafita
 * @since 5.0.0
 *
 */
public class TestSuperPosition {

	private static final Logger logger = LoggerFactory
			.getLogger(TestSuperPosition.class);

	private List<Point3d[]> cloud1;
	private List<Point3d[]> cloud2;

	private AxisAngle4d rotAxis;
	private Vector3d translation;
	private Matrix4d transform;

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

		rotAxis = new AxisAngle4d(0.440, 0.302, 0.845, 1.570);
		translation = new Vector3d(0.345, 2.453, 5.324);
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

			logger.error(String.format("Transformation Matrix %d points: "
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

			logger.error(String.format("RMSD %d points: SVD time %d us, "
					+ "Quat time: %d us, QCP time: %d us", cloud1.get(c).length,
					svdTime, quatTime, qcpTime));

			// Check that the returned RMSDs are equal
			assertEquals(svdrmsd, quatrmsd, 0.001);
			assertEquals(svdrmsd, qcprmsd, 0.001);
		}
	}

}
