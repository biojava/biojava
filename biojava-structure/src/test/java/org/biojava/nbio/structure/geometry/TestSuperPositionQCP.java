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

import org.biojava.nbio.structure.SVDSuperimposer;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.geometry.SuperPosition;
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
public class TestSuperPositionQCP {

	private static final Logger logger = LoggerFactory
			.getLogger(TestSuperPositionQCP.class);

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
	 * Test method to obtain the transformation matrix from superposition
	 * {@link SuperPositionQCP#getTransformationMatrix()},
	 * {@link SuperPositionQCP#getRmsd()}.
	 * 
	 * @throws StructureException
	 */
	@Test
	public void testTransformationMatrix() throws StructureException {

		for (int c = 0; c < cloud1.size(); c++) {
			// Use SVD superposition to obtain the optimal transformation matrix
			long svdStart = System.nanoTime();
			SVDSuperimposer svd = new SVDSuperimposer(cloud1.get(c),
					cloud2.get(c));
			Matrix4d svdTransform = svd.getTransformation();
			long svdTime = (System.nanoTime() - svdStart) / 1000;

			Point3d[] c2c = CalcPoint.clonePoint3dArray(cloud2.get(c));

			// Use SuperPosition to obtain the optimal transformation matrix
			long spStart = System.nanoTime();
			Matrix4d spTransform = SuperPosition.superposeWithTranslation(c2c,
					cloud1.get(c));
			long spTime = (System.nanoTime() - spStart) / 1000;

			// Use QCP algorithm to get the optimal transformation matrix
			SuperPositionQCP qcp = new SuperPositionQCP();
			qcp.set(cloud2.get(c), cloud1.get(c));
			long qcpStart = System.nanoTime();
			Matrix4d qcpTransform = qcp.getTransformationMatrix();
			long qcpTime = (System.nanoTime() - qcpStart) / 1000;

			logger.info(String.format("Transformation Matrix %d points: "
					+ "SVD time %d us, SP time: %d us, QCP time: %d us",
					cloud1.get(c).length, svdTime, spTime, qcpTime));

			// Check that the transformation matrix was recovered
			assertTrue(transform.epsilonEquals(svdTransform, 0.01));
			assertTrue(transform.epsilonEquals(spTransform, 0.01));
			assertTrue(transform.epsilonEquals(qcpTransform, 0.01));
		}

	}

	/**
	 * Test method to obtain the RMSD of a superposition
	 * {@link SuperPositionQCP#getRmsd()}.
	 * 
	 * @throws StructureException
	 */
	@Test
	public void testRMSD() throws StructureException {

		for (int c = 0; c < cloud1.size(); c++) {
			Point3d[] c2c = CalcPoint.clonePoint3dArray(cloud2.get(c));

			// Use SVD superposition to obtain the RMSD
			long svdStart = System.nanoTime();
			SVDSuperimposer svd = new SVDSuperimposer(cloud1.get(c),
					cloud2.get(c));
			Matrix4d svdTransform = svd.getTransformation();
			CalcPoint.transform(svdTransform, c2c);
			double svdrmsd = SuperPosition.rmsd(cloud1.get(c), c2c);
			long svdTime = (System.nanoTime() - svdStart) / 1000;

			c2c = CalcPoint.clonePoint3dArray(cloud2.get(c));

			// Use SVD superposition to obtain the RMSD
			long spStart = System.nanoTime();
			SuperPosition.superposeWithTranslation(c2c, cloud1.get(c));
			double sprmsd = SuperPosition.rmsd(c2c, cloud1.get(c));
			long spTime = (System.nanoTime() - spStart) / 1000;

			// Use QCP algorithm to obtain the RMSD
			SuperPositionQCP qcp = new SuperPositionQCP();
			qcp.set(cloud2.get(c), cloud1.get(c));
			long qcpStart = System.nanoTime();
			double qcprmsd = qcp.getRmsd();
			long qcpTime = (System.nanoTime() - qcpStart) / 1000;

			logger.info(String.format("RMSD %d points: SVD time %d us, "
					+ "SP time: %d us, QCP time: %d us", cloud1.get(c).length,
					svdTime, spTime, qcpTime));

			// Check that the returned RMSDs are equal
			assertEquals(svdrmsd, sprmsd, 0.001);
			assertEquals(svdrmsd, qcprmsd, 0.001);
		}
	}

	@Test
	public void testSymmetryQCP() {

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

		// Use SP superposition to obtain the RMSD
		long spStart = System.nanoTime();
		SuperPosition.superposeWithTranslation(set1, set2);
		double sprmsd = SuperPosition.rmsd(set1, set2);
		long spTime = (System.nanoTime() - spStart) / 1000;

		set2 = CalcPoint.clonePoint3dArray(set1);

		// Use QCP algorithm to get the RMSD
		SuperPositionQCP qcp = new SuperPositionQCP();
		qcp.set(set1, set2);
		long qcpStart = System.nanoTime();
		double qcprmsd = qcp.getRmsd();
		long qcpTime = (System.nanoTime() - qcpStart) / 1000;

		logger.info(String.format("RMSD Symmetry: SP time: %d us"
				+ ", QCP time: %d us", spTime, qcpTime));

		// Check that the returned RMSDs are equal
		assertEquals(sprmsd, qcprmsd, 0.001);

	}

}
