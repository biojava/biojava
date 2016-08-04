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

import org.biojava.nbio.structure.jama.EigenvalueDecomposition;
import org.biojava.nbio.structure.jama.Matrix;

import javax.vecmath.Matrix3d;
import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

import java.util.ArrayList;
import java.util.List;

/**
 * The moment of inertia, otherwise known as the angular mass or rotational
 * inertia, of a rigid body determines the torque needed for a desired angular
 * acceleration about a rotational axis. It depends on the body's mass
 * distribution and the axis chosen, with larger moments requiring more torque
 * to change the body's rotation.
 * <p>
 * More in https://en.wikipedia.org/wiki/Moment_of_inertia.
 * 
 * @author Peter Rose
 * @author Aleix Lafita
 * 
 */
public class MomentsOfInertia {

	private List<Point3d> points = new ArrayList<Point3d>();
	private List<Double> masses = new ArrayList<Double>();

	private boolean modified = true;

	private double[] principalMomentsOfInertia = new double[3];
	private Vector3d[] principalAxes = new Vector3d[3];

	private Matrix3d orientation = new Matrix3d();

	public enum SymmetryClass {
		LINEAR, PROLATE, OBLATE, SYMMETRIC, ASYMMETRIC
	};

	/** Creates a new empty instance of MomentsOfInertia */
	public MomentsOfInertia() {
	}

	public void addPoint(Point3d point, double mass) {
		points.add(point);
		masses.add(mass);
		modified = true;
	}

	public Point3d getCenterOfMass() {

		if (points.size() == 0) {
			throw new IllegalStateException(
					"MomentsOfInertia: no points defined");
		}

		Point3d center = new Point3d();
		double totalMass = 0.0;

		for (int i = 0, n = points.size(); i < n; i++) {
			double mass = masses.get(i);
			totalMass += mass;
			center.scaleAdd(mass, points.get(i), center);
		}
		center.scale(1.0 / totalMass);
		return center;
	}

	public double[] getPrincipalMomentsOfInertia() {

		if (modified) {
			diagonalizeTensor();
			modified = false;
		}
		return principalMomentsOfInertia;
	}

	/**
	 * The principal axes of intertia
	 * 
	 * @return
	 */
	public Vector3d[] getPrincipalAxes() {

		if (modified) {
			diagonalizeTensor();
			modified = false;
		}
		return principalAxes;
	}

	/**
	 * The orientation Matrix is a 3x3 Matrix with a column for each principal
	 * axis. It represents the orientation (rotation) of the principal axes with
	 * respect to the axes of the coordinate system (unit vectors [1,0,0],
	 * [0,1,0] and [0,0,1]).
	 * <p>
	 * The orientation matrix indicates the rotation to bring the coordinate
	 * axes to the principal axes, in this direction.
	 * 
	 * @return the orientation Matrix as a Matrix3d object
	 */
	public Matrix3d getOrientationMatrix() {

		if (modified) {
			diagonalizeTensor();
			modified = false;
		}
		return orientation;
	}

	/**
	 * The effective value of this distance for a certain body is known as its
	 * radius of / gyration with respect to the given axis. The radius of
	 * gyration corresponding to Ijj / is defined as /
	 * http://www.eng.auburn.edu/~marghitu/MECH2110/C_4.pdf / radius of gyration
	 * k(j) = sqrt(I(j)/m)
	 */
	public double[] getElipsisRadii() {
		if (modified) {
			diagonalizeTensor();
			modified = false;
		}
		double m = 0;
		for (int i = 0, n = points.size(); i < n; i++) {
			m += masses.get(i);
		}
		double[] r = new double[3];
		for (int i = 0; i < 3; i++) {
			r[i] = Math.sqrt(principalMomentsOfInertia[i] / m);
		}
		return r;
	}

	public double getRadiusOfGyration() {
		Point3d c = getCenterOfMass();
		Point3d t = new Point3d();
		double sum = 0;
		for (int i = 0, n = points.size(); i < n; i++) {
			t.set(points.get(i));
			sum += t.distanceSquared(c);
		}
		sum /= points.size();
		return Math.sqrt(sum);
	}

	public SymmetryClass getSymmetryClass(double threshold) {
		if (modified) {
			diagonalizeTensor();
			modified = false;
		}
		double ia = principalMomentsOfInertia[0];
		double ib = principalMomentsOfInertia[1];
		double ic = principalMomentsOfInertia[2];
		boolean c1 = (ib - ia) / (ib + ia) < threshold;
		boolean c2 = (ic - ib) / (ic + ib) < threshold;

		if (c1 && c2) {
			return SymmetryClass.SYMMETRIC;
		}
		if (c1) {
			return SymmetryClass.OBLATE;
		}
		if (c2) {
			return SymmetryClass.PROLATE;
		}
		return SymmetryClass.ASYMMETRIC;
	}

	public double symmetryCoefficient() {
		if (modified) {
			diagonalizeTensor();
			modified = false;
		}
		double ia = principalMomentsOfInertia[0];
		double ib = principalMomentsOfInertia[1];
		double ic = principalMomentsOfInertia[2];
		double c1 = 1.0f - (ib - ia) / (ib + ia);
		double c2 = 1.0f - (ic - ib) / (ic + ib);
		return Math.max(c1, c2);
	}

	public double getAsymmetryParameter(double threshold) {
		if (modified) {
			diagonalizeTensor();
			modified = false;
		}
		if (getSymmetryClass(threshold).equals(SymmetryClass.SYMMETRIC)) {
			return 0.0;
		}
		double a = 1.0 / principalMomentsOfInertia[0];
		double b = 1.0 / principalMomentsOfInertia[1];
		double c = 1.0 / principalMomentsOfInertia[2];
		return (2 * b - a - c) / (a - c);
	}

	public double[][] getInertiaTensor() {

		Point3d p = new Point3d();
		double[][] tensor = new double[3][3];

		// calculate the inertia tensor at center of mass
		Point3d com = getCenterOfMass();

		for (int i = 0, n = points.size(); i < n; i++) {
			double mass = masses.get(i);
			p.sub(points.get(i), com);

			tensor[0][0] += mass * (p.y * p.y + p.z * p.z);
			tensor[1][1] += mass * (p.x * p.x + p.z * p.z);
			tensor[2][2] += mass * (p.x * p.x + p.y * p.y);

			tensor[0][1] -= mass * p.x * p.y;
			tensor[0][2] -= mass * p.x * p.z;
			tensor[1][2] -= mass * p.y * p.z;
		}

		tensor[1][0] = tensor[0][1];
		tensor[2][0] = tensor[0][2];
		tensor[2][1] = tensor[1][2];

		return tensor;
	}

	private void diagonalizeTensor() {

		Matrix m = new Matrix(getInertiaTensor());
		EigenvalueDecomposition eig = m.eig();

		principalMomentsOfInertia = eig.getRealEigenvalues();
		double[][] eigenVectors = eig.getV().getArray();

		// Get the principal axes from the eigenVectors
		principalAxes[0] = new Vector3d(eigenVectors[0][0], eigenVectors[1][0],
				eigenVectors[2][0]);
		principalAxes[1] = new Vector3d(eigenVectors[0][1], eigenVectors[1][1],
				eigenVectors[2][1]);
		principalAxes[2] = new Vector3d(eigenVectors[0][2], eigenVectors[1][2],
				eigenVectors[2][2]);

		// Convert the principal axes into a rotation matrix
		for (int i = 0; i < 3; i++)
			orientation.setColumn(i, principalAxes[i]);
		orientation.negate();

	}
}
