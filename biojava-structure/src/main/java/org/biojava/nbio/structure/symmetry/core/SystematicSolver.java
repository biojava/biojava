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
import org.biojava.nbio.structure.symmetry.utils.PermutationGenerator;

import javax.vecmath.AxisAngle4d;
import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;
import javax.vecmath.Quat4d;
import javax.vecmath.Vector3d;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;


/**
 *
 * @author Peter
 */
public class SystematicSolver implements QuatSymmetrySolver {
	private QuatSymmetrySubunits subunits = null;
	private QuatSymmetryParameters parameters = null;

	private Point3d[] originalCoords = null;
	private Point3d[] transformedCoords = null;
	private RotationGroup rotations = new RotationGroup();
	private Vector3d centroid = new Vector3d();
	private Matrix4d centroidInverse = new Matrix4d();
	private Set<List<Integer>> hashCodes = new HashSet<List<Integer>>();

	public SystematicSolver(QuatSymmetrySubunits subunits, QuatSymmetryParameters parameters) {
		if (subunits.getSubunitCount()== 2) {
			throw new IllegalArgumentException("SystematicSolver cannot be applied to subunits with 2 centers");
		}
		this.subunits = subunits;
		this.parameters = parameters;
	}

	@Override
	public RotationGroup getSymmetryOperations() {
		if (rotations.getOrder() == 0) {
			solve();
			rotations.complete();
		}
		return rotations;
	}

	private void solve() {
		initialize();
		int n = subunits.getSubunitCount();
		PermutationGenerator g = new PermutationGenerator(n);

		// loop over all permutations
		while (g.hasMore()) {
			int[] perm = g.getNext();
			List<Integer> permutation = new ArrayList<Integer>(perm.length);
			for (int j = 0; j < n; j++) {
				permutation.add(perm[j]);
			}

			if (! isValidPermutation(permutation)) {
				continue;
			}

			boolean newPermutation = evaluatePermutation(permutation);
			if (newPermutation) {
				completeRotationGroup();
			}

			if (rotations.getOrder() >= subunits.getSubunitCount()) {
				return;
			}
		}
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

	private void completeRotationGroup() {
		PermutationGroup g = new PermutationGroup();
		for (int i = 0; i < rotations.getOrder(); i++) {
			Rotation s = rotations.getRotation(i);
			g.addPermutation(s.getPermutation());
		}
		g.completeGroup();

//   	System.out.println("Completing rotation group from: " +symmetryOperations.getSymmetryOperationCount() + " to " + g.getPermutationCount());

		// the group is complete, nothing to do
		if (g.getOrder() == rotations.getOrder()) {
			return;
		}

//  	System.out.println("complete group: " +  rotations.getOrder() +"/" + g.getOrder());
		// try to complete the group
		for (int i = 0; i < g.getOrder(); i++) {
			List<Integer> permutation = g.getPermutation(i);
			if (isValidPermutation(permutation)) {
				  // perform permutation of subunits
				evaluatePermutation(permutation);
			}
		}
	}

	private boolean isValidPermutation(List<Integer> permutation) {
		if (permutation.size() == 0) {
			return false;
		}

		// if this permutation is a duplicate, return false
		if (hashCodes.contains(permutation)) {
			return false;
		}

		// check if permutation is pseudosymmetric
		if (! isAllowedPermuation(permutation)) {
			return false;
		}

		// get fold and make sure there is only one E (fold=1) permutation
		int fold = PermutationGroup.getOrder(permutation);
		if (rotations.getOrder() > 1 && fold == 1) {
			return false;
		}
		if (fold == 0 || subunits.getSubunitCount() % fold != 0) {
			return false;
		}

		// if this permutation is a duplicate, returns false
		return hashCodes.add(permutation);
	}

	private boolean isAllowedPermuation(List<Integer> permutation) {
		List<Integer> seqClusterId = subunits.getClusterIds();
		for (int i = 0; i < permutation.size(); i++) {
			int j = permutation.get(i);
			if (seqClusterId.get(i) != seqClusterId.get(j)) {
				return false;
			}
		}
		return true;
	}

	private boolean evaluatePermutation(List<Integer> permutation) {
		// permutate subunits
		for (int j = 0, n = subunits.getSubunitCount(); j < n; j++) {
			transformedCoords[j].set(originalCoords[permutation.get(j)]);
		}

		int fold = PermutationGroup.getOrder(permutation);
		
		// TODO implement this piece of code using at origin superposition
		Quat4d quat = UnitQuaternions.relativeOrientation(
				originalCoords, transformedCoords);
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
		
		CalcPoint.transform(transformation, transformedCoords);
		
		double subunitRmsd = CalcPoint.rmsd(transformedCoords, originalCoords);

		if (subunitRmsd <parameters.getRmsdThreshold()) {
			// transform to original coordinate system
			combineWithTranslation(transformation);
			QuatSymmetryScores scores = QuatSuperpositionScorer.calcScores(subunits, transformation, permutation);
			if (scores.getRmsd() < 0.0 || scores.getRmsd() > parameters.getRmsdThreshold()) {
				return false;
			}

			scores.setRmsdCenters(subunitRmsd);
			Rotation symmetryOperation = createSymmetryOperation(permutation, transformation, axisAngle, fold, scores);
			rotations.addRotation(symmetryOperation);
			return true;
		}
		return false;
	}

	private void initialize() {
		// translation to centered coordinate system
		centroid = new Vector3d(subunits.getCentroid());

		// translation back to original coordinate system
		Vector3d reverse = new Vector3d(centroid);
		reverse.negate();
		centroidInverse.set(reverse);
		// Make sure matrix element m33 is 1.0. An old version vecmath did not set this element.
		centroidInverse.setElement(3, 3, 1.0);

		List<Point3d> centers = subunits.getCenters();
		int n = subunits.getSubunitCount();

		originalCoords = new Point3d[n];
		transformedCoords = new Point3d[n];

		for (int i = 0; i < n; i++) {
			originalCoords[i] = centers.get(i);
			transformedCoords[i] = new Point3d();
		}
	}
}
