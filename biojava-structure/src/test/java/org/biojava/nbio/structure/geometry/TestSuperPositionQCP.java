package org.biojava.nbio.structure.geometry;

import static org.junit.Assert.*;

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

	private static final Logger logger = LoggerFactory
			.getLogger(TestSuperPositionQCP.class);

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

		logger.info(String.format("RMSD Symmetry: Quat time: %d us"
				+ ", QCP time: %d us", quatTime, qcpTime));

		// Check that the returned RMSDs are equal
		assertEquals(quatrmsd, qcprmsd, 0.001);

	}

}
