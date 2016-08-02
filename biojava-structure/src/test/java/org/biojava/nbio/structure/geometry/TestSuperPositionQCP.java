package org.biojava.nbio.structure.geometry;

import static org.junit.Assert.*;

import java.util.Random;

import javax.vecmath.AxisAngle4d;
import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.AtomImpl;
import org.biojava.nbio.structure.Calc;
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

	private Atom[] cloud1;
	private Atom[] cloud2;
	private Atom[] cloud2clone;

	private AxisAngle4d axis;
	private Matrix4d transform;

	@Before
	public void setUp() {

		// Generate two big random clouds of points
		Random rnd = new Random();

		cloud1 = new Atom[5000];
		cloud2 = new Atom[5000];
		cloud2clone = new Atom[5000];

		axis = new AxisAngle4d(0.440, 0.302, 0.845, 1.570);
		transform = new Matrix4d();
		transform.set(axis);

		for (int p = 0; p < 5000; p++) {

			Atom a = new AtomImpl();
			a.setCoords(new double[] { rnd.nextInt(100), rnd.nextInt(50),
					rnd.nextInt(150) });
			cloud1[p] = a;

			Atom b = new AtomImpl();
			// Add some noise
			b.setCoords(new double[] { a.getX() + rnd.nextDouble(),
					a.getY() + rnd.nextDouble(), a.getZ() + rnd.nextDouble() });
			Calc.transform(b, transform);
			cloud2[p] = b;

			Atom c = (Atom) b.clone();
			cloud2clone[p] = c;
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

		// Use SVD superposition to obtain the optimal transformation matrix
		long svdStart = System.nanoTime();
		SVDSuperimposer svd = new SVDSuperimposer(cloud1, cloud2);
		Matrix4d svdTransform = svd.getTransformation();
		long svdTime = (System.nanoTime() - svdStart) / 1000;
		AxisAngle4d svdAxis = new AxisAngle4d();
		svdAxis.set(svdTransform);

		// Use SuperPosition to obtain the optimal transformation matrix
		Point3d[] cloud1p = Calc.atomsToPoints(cloud1);
		Point3d[] cloud2p = Calc.atomsToPoints(cloud2);
		long spStart = System.nanoTime();
		Matrix4d spTransform = SuperPosition.superpose(cloud1p, cloud2p);
		long spTime = (System.nanoTime() - spStart) / 1000;
		AxisAngle4d spAxis = new AxisAngle4d();
		spAxis.set(spTransform);

		// Use QCP algorithm to get the optimal transformation matrix
		SuperPositionQCP qcp = new SuperPositionQCP();
		qcp.set(Calc.atomsToPoints(cloud1), Calc.atomsToPoints(cloud2));
		long qcpStart = System.nanoTime();
		Matrix4d qcpTransform = qcp.getTransformationMatrix();
		long qcpTime = (System.nanoTime() - qcpStart) / 1000;
		AxisAngle4d qcpAxis = new AxisAngle4d();
		qcpAxis.set(qcpTransform);

		logger.info(String
				.format("Transformation Matrix: SVD time %d us"
						+ ", SP time: %d us, QCP time: %d us", svdTime, spTime,
						qcpTime));

		// Check that the rotation angle was recovered in both
		assertEquals(svdAxis.angle, axis.angle, 0.01);
		assertEquals(spAxis.angle, axis.angle, 0.01);
		assertEquals(qcpAxis.angle, axis.angle, 0.01);

		// Check that the axis vector was recovered
		AxisAngle4d negAxis = new AxisAngle4d(-axis.x, -axis.y, -axis.z,
				axis.angle);

		assertTrue(axis.epsilonEquals(svdAxis, 0.01)
				|| negAxis.epsilonEquals(svdAxis, 0.01));

		assertTrue(axis.epsilonEquals(spAxis, 0.01)
				|| negAxis.epsilonEquals(spAxis, 0.01));
		
		assertTrue(axis.epsilonEquals(qcpAxis, 0.01)
				|| negAxis.epsilonEquals(qcpAxis, 0.01));

	}

	/**
	 * Test method to obtain the RMSD of a superposition
	 * {@link SuperPositionQCP#getRmsd()}.
	 * 
	 * @throws StructureException
	 */
	@Test
	public void testRMSD() throws StructureException {

		// Use SVD superposition to obtain the RMSD
		long svdStart = System.nanoTime();
		SVDSuperimposer svd = new SVDSuperimposer(cloud1, cloud2);
		Matrix4d svdTransform = svd.getTransformation();
		Calc.transform(cloud2clone, svdTransform);
		double svdrmsd = SVDSuperimposer.getRMS(cloud1, cloud2clone);
		long svdTime = (System.nanoTime() - svdStart) / 1000;

		// Use QCP algorithm to get the optimal transformation matrix
		SuperPositionQCP qcp = new SuperPositionQCP();
		qcp.set(Calc.atomsToPoints(cloud1), Calc.atomsToPoints(cloud2));
		long qcpStart = System.nanoTime();
		double qcprmsd = qcp.getRmsd();
		long qcpTime = (System.nanoTime() - qcpStart) / 1000;

		logger.info(String.format("RMSD: SVD time %d us" + ", QCP time: %d us",
				svdTime, qcpTime));

		// Check that the returned RMSDs are equal
		assertEquals(svdrmsd, qcprmsd, 0.01);
	}

}
