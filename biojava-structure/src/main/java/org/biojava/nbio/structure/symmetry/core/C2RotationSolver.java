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

package org.biojava.nbio.structure.symmetry.core;

import org.biojava.nbio.structure.geometry.CalcPoint;
import org.biojava.nbio.structure.geometry.UnitQuaternions;

import javax.vecmath.AxisAngle4d;
import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;
import javax.vecmath.Quat4d;
import javax.vecmath.Vector3d;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 *
 * @author Peter
 */
public class C2RotationSolver implements QuatSymmetrySolver {
	private QuatSymmetrySubunits subunits = null;
	private QuatSymmetryParameters parameters = null;
	private Vector3d centroid = new Vector3d();
	private Matrix4d centroidInverse = new Matrix4d();

	private RotationGroup rotations = new RotationGroup();


	public C2RotationSolver(QuatSymmetrySubunits subunits, QuatSymmetryParameters parameters) {
		if (subunits.getSubunitCount() != 2) {
			throw new IllegalArgumentException("C2RotationSolver can only be applied to cases with 2 centers");
		}
		this.subunits = subunits;
		this.parameters = parameters;
	}

	@Override
	public RotationGroup getSymmetryOperations() {
		if (rotations.getOrder() == 0) {
			solve();
		}
		return rotations;
	}

	private void solve() {
		initialize();
		Vector3d trans = new Vector3d(subunits.getCentroid());
		trans.negate();
		List<Point3d[]> traces = subunits.getTraces();

//		Point3d[] x = SuperPosition.clonePoint3dArray(traces.get(0));
//		SuperPosition.center(x);
//		Point3d[] y = SuperPosition.clonePoint3dArray(traces.get(1));
//		SuperPosition.center(y);

		Point3d[] x = CalcPoint.clonePoint3dArray(traces.get(0));
		CalcPoint.translate(trans, x);
		Point3d[] y = CalcPoint.clonePoint3dArray(traces.get(1));
		CalcPoint.translate(trans, y);

		// TODO implement this piece of code using at origin superposition
		Quat4d quat = UnitQuaternions.relativeOrientation(
				x, y);
		AxisAngle4d axisAngle = new AxisAngle4d();
		Matrix4d transformation = new Matrix4d();
		
		transformation.set(quat);
		axisAngle.set(quat);
		
		Vector3d axis = new Vector3d(axisAngle.x, axisAngle.y, axisAngle.z);
		if (axis.lengthSquared() < 1.0E-6) {
			axisAngle.x = 0;
			axisAngle.y = 0;
			axisAngle.z = 1;
			axisAngle.angle = 0;
		} else {
			axis.normalize();
			axisAngle.x = axis.x;
			axisAngle.y = axis.y;
			axisAngle.z = axis.z;
		}
		
		CalcPoint.transform(transformation, y);

		// if rmsd or angle deviation is above threshold, stop
		double angleThresholdRadians = Math.toRadians(parameters.getAngleThreshold());
		double deltaAngle = Math.abs(Math.PI-axisAngle.angle);

		if (deltaAngle > angleThresholdRadians) {
			rotations.setC1(subunits.getSubunitCount());
			return;
		}

		// add unit operation
		addEOperation();

		// add C2 operation
		int fold = 2;
		combineWithTranslation(transformation);
		List<Integer> permutation = Arrays.asList(1,0);
		QuatSymmetryScores scores = QuatSuperpositionScorer.calcScores(subunits, transformation, permutation);
		scores.setRmsdCenters(0.0); // rmsd for superposition of two subunits centers is zero by definition

		if (scores.getRmsd() > parameters.getRmsdThreshold() || deltaAngle > angleThresholdRadians) {
			rotations.setC1(subunits.getSubunitCount());
			return;
		}

		Rotation symmetryOperation = createSymmetryOperation(permutation, transformation, axisAngle, fold, scores);
		rotations.addRotation(symmetryOperation);
	}

	private void addEOperation() {
		List<Integer> permutation = Arrays.asList(new Integer[]{0,1});
		Matrix4d transformation = new Matrix4d();
		transformation.setIdentity();
		combineWithTranslation(transformation);
		AxisAngle4d axisAngle = new AxisAngle4d();
		QuatSymmetryScores scores = new QuatSymmetryScores();
		int fold = 1; // ??
		Rotation rotation = createSymmetryOperation(permutation, transformation, axisAngle, fold, scores);
		rotations.addRotation(rotation);
	}

	/**
	 * Adds translational component to rotation matrix
	 * @param rotTrans
	 * @param rotation
	 * @return
	 */
	private void combineWithTranslation(Matrix4d rotation) {
		rotation.setTranslation(centroid);
		rotation.mul(rotation, centroidInverse);
	}

	private Rotation createSymmetryOperation(List<Integer> permutation, Matrix4d transformation, AxisAngle4d axisAngle, int fold, QuatSymmetryScores scores) {
		Rotation s = new Rotation();
		s.setPermutation(new ArrayList<Integer>(permutation));
		s.setTransformation(new Matrix4d(transformation));
		s.setAxisAngle(new AxisAngle4d(axisAngle));
		s.setFold(fold);
		s.setScores(scores);
		return s;
	}

	private void initialize() {
		// translation to centered coordinate system
		centroid = new Vector3d(subunits.getCentroid());
	   // translation back to original coordinate system
		Vector3d reverse = new Vector3d(centroid);
		reverse.negate();
		centroidInverse.set(reverse);
//        // On LINUX there seems to be a bug with vecmath, and element m33 is zero. Here we make sure it's 1.
		centroidInverse.setElement(3, 3, 1.0);
	}

}
