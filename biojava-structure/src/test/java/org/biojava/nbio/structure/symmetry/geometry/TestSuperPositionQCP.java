package org.biojava.nbio.structure.symmetry.geometry;

import static org.junit.Assert.*;

import java.util.Random;

import javax.vecmath.AxisAngle4d;
import javax.vecmath.Matrix4d;
import javax.vecmath.Vector3d;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.AtomImpl;
import org.biojava.nbio.structure.Calc;
import org.biojava.nbio.structure.SVDSuperimposer;
import org.biojava.nbio.structure.StructureException;
import org.junit.Test;

/**
 * Test the Quaternion-Based Characteristic Polynomial {@link SuperPositionQCP}
 * algorithm for RMSD and Superposition calculations.
 * 
 * @author Aleix Lafita
 * @since 5.0.0
 *
 */
public class TestSuperPositionQCP {

	/**
	 * Test all the methods for superposition
	 * {@link SuperPositionQCP#getTransformationMatrix()}.
	 * 
	 * @throws StructureException
	 */
	@Test
	public void testSuperposition() throws StructureException {

		// Generate two random clouds of points
		Random rnd = new Random();

		Atom[] cloud1 = new Atom[50];
		Atom[] cloud2 = new Atom[50];
		Atom[] cloud2clone = new Atom[50];

		Matrix4d transform = new Matrix4d();
		transform.setRotation(new AxisAngle4d(1, 0, 0, 0.78));
		transform.setTranslation(new Vector3d(0, 1, 0));

		for (int p = 0; p < 50; p++) {

			Atom a = new AtomImpl();
			a.setCoords(new double[] { rnd.nextInt(100), rnd.nextInt(100),
					rnd.nextInt(100) });
			cloud1[p] = a;

			// Add some noise
			Atom b = new AtomImpl();
			b.setCoords(new double[] { a.getX() + rnd.nextDouble(),
					a.getY() + rnd.nextDouble(), a.getZ() + rnd.nextDouble() });
			Calc.transform(b, transform);
			cloud2[p] = b;

			Atom c = (Atom) b.clone();
			cloud2clone[p] = c;
		}

		// Use SVD superposition to obtain the optimal transformation matrix
		SVDSuperimposer svd = new SVDSuperimposer(cloud1, cloud2);
		Matrix4d svdTransform = svd.getTransformation();
		AxisAngle4d svdAxis = new AxisAngle4d();
		svdAxis.set(svdTransform);
		//Calc.transform(cloud2clone, svdTransform);
		double svdrmsd = SVDSuperimposer.getRMS(cloud1, cloud2clone);

		// Use QCP algorithm to get the optimal transformation matrix
		SuperPositionQCP qcp = new SuperPositionQCP();
		qcp.set(Calc.atomsToPoints(cloud1), Calc.atomsToPoints(cloud2));
		Matrix4d qcpTransform = qcp.getTransformationMatrix();
		AxisAngle4d qcpAxis = new AxisAngle4d();
		qcpAxis.set(qcpTransform);
		double qcprmsd = qcp.getRmsd();

		//assertTrue(svdTransform.epsilonEquals(qcpTransform, 0.01));

	}
}
