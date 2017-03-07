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

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.vecmath.AxisAngle4d;
import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;
import javax.vecmath.Quat4d;
import javax.vecmath.Vector3d;

import org.biojava.nbio.structure.geometry.CalcPoint;
import org.biojava.nbio.structure.geometry.MomentsOfInertia;
import org.biojava.nbio.structure.geometry.UnitQuaternions;
import org.biojava.nbio.structure.symmetry.geometry.DistanceBox;
import org.biojava.nbio.structure.symmetry.geometry.SphereSampler;


/**
 *
 * @author Peter
 */
public class RotationSolver implements QuatSymmetrySolver {
	private QuatSymmetrySubunits subunits = null;
	private QuatSymmetryParameters parameters = null;

	private double distanceThreshold = 0.0f;
	private DistanceBox<Integer> box = null;
	private Vector3d centroid = new Vector3d();
	private Matrix4d centroidInverse = new Matrix4d();
	private Point3d[] originalCoords = null;
	private Point3d[] transformedCoords = null;
	// Cache whether a permutation is invalid (null) vs has been added to rotations
	private Map<List<Integer>,Rotation> evaluatedPermutations = new HashMap<>();

	private RotationGroup rotations = new RotationGroup();

	public RotationSolver(QuatSymmetrySubunits subunits, QuatSymmetryParameters parameters) {
		if (subunits.getSubunitCount()== 2) {
			throw new IllegalArgumentException("RotationSolver cannot be applied to subunits with 2 centers");
		}
		this.subunits = subunits;
		this.parameters = parameters;
	}

	@Override
	public RotationGroup getSymmetryOperations() {
		if (rotations.getOrder() == 0) {
			solve();
			completeRotationGroup();
			rotations.complete();
		}

		return rotations;
	}

	private void solve() {
		initialize();

		int maxSymOps = subunits.getSubunitCount();

		boolean isSpherical = isSpherical();

		// for cases with icosahedral symmetry n cannot be higher than 60, should check for spherical symmetry here
		// isSpherical check added 08-04-11
		if (maxSymOps % 60 == 0 && isSpherical) {
			maxSymOps = 60;
		 }

		AxisAngle4d sphereAngle = new AxisAngle4d();
		Matrix4d transformation = new Matrix4d();

		int n = subunits.getSubunitCount();

		int sphereCount = SphereSampler.getSphereCount();

		List<Double> angles = getAngles();

		for (int i = 0; i < sphereCount; i++) {
			// Sampled orientation
			//TODO The SphereSampler samples 4D orientation space. We really
			// only need to sample 3D unit vectors, since we use limited
			// angles. -SB
			SphereSampler.getAxisAngle(i, sphereAngle);

			// Each valid rotation angle
			for (double angle : angles) {
				// apply rotation
				sphereAngle.angle = angle;
				transformation.set(sphereAngle);
				// Make sure matrix element m33 is 1.0. It's not on Linux.
				transformation.setElement(3, 3, 1.0);
				for (int j = 0; j < n; j++) {
					transformedCoords[j].set(originalCoords[j]);
					transformation.transform(transformedCoords[j]);
				}

				// get permutation of subunits and check validity/uniqueness
				List<Integer> permutation = getPermutation();
	//              System.out.println("Rotation Solver: permutation: " + i + ": " + permutation);

				// check if novel
				if ( evaluatedPermutations.containsKey(permutation)) {
					continue; //either invalid or already added
				}

				Rotation newPermutation = isValidPermutation(permutation);
				if (newPermutation != null) {
					completeRotationGroup(newPermutation);
				}

				// check if all symmetry operations have been found.
				if (rotations.getOrder() >= maxSymOps) {
					return;
				}
			}
		}
	}

	/**
	 * Combine current rotations to make all possible permutations.
	 * If these are all valid, add them to the rotations
	 * @param additionalRots Additional rotations we are considering adding to this.rotations
	 * @return whether the rotations were valid and added
	 */
	private boolean completeRotationGroup(Rotation... additionalRots) {
		PermutationGroup g = new PermutationGroup();
		for (Rotation s : rotations) {
			g.addPermutation(s.getPermutation());
		}
		for( Rotation s : additionalRots) {
			g.addPermutation(s.getPermutation());
			// inputs should not have been added already
			assert evaluatedPermutations.get(s.getPermutation()) == null;
		}
		g.completeGroup();

		// the group is complete, nothing to do
		if (g.getOrder() == rotations.getOrder()+additionalRots.length) {
			for( Rotation s : additionalRots) {
				addRotation(s);
			}
			return true;
		}

		// try to complete the group
		List<Rotation> newRots = new ArrayList<>(g.getOrder());
		// First, quick check for whether they're allowed
		for (List<Integer> permutation : g) {
			if( evaluatedPermutations.containsKey(permutation)) {
				Rotation rot = evaluatedPermutations.get(permutation);
				if( rot == null ) {
					return false;
				}
			} else {
				if( ! isAllowedPermutation(permutation)) {
					return false;
				}
			}
		}
		// Slower check including the superpositions
		for (List<Integer> permutation : g) {
			Rotation rot;
			if( evaluatedPermutations.containsKey(permutation)) {
				rot = evaluatedPermutations.get(permutation);
			} else {
				rot = isValidPermutation(permutation);
			}

			if( rot == null ) {
				// if any induced rotation is invalid, abort
				return false;
			}
			if(!evaluatedPermutations.containsKey(permutation)){ //novel
				newRots.add(rot);
			}
		}
		// Add rotations
		for( Rotation rot : newRots) {
			addRotation(rot);
		}
		return true;
	}

	private void addRotation(Rotation rot) {
		evaluatedPermutations.put(rot.getPermutation(),rot);
		rotations.addRotation(rot);
	}

	/**
	 * Superimpose subunits based on the given permutation. Then check whether
	 * the superposition passes RMSD thresholds and create a Rotation to
	 * represent it if so.
	 * @param permutation A list specifying which subunits should be aligned by the current transformation
	 * @return A Rotation representing the permutation, or null if the superposition did not meet thresholds.
	 */
	private Rotation superimposePermutation(List<Integer> permutation) {
		// permutate subunits
		for (int j = 0, n = subunits.getSubunitCount(); j < n; j++) {
			transformedCoords[j].set(originalCoords[permutation.get(j)]);
		}

		int fold = PermutationGroup.getOrder(permutation);

		// get optimal transformation and axisangle by subunit superposition
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

		if (subunitRmsd < parameters.getRmsdThreshold()) {
			combineWithTranslation(transformation);

			// evaluate superposition of CA traces
			QuatSymmetryScores scores = QuatSuperpositionScorer.calcScores(subunits, transformation, permutation);
			if (scores.getRmsd() < 0.0 || scores.getRmsd() > parameters.getRmsdThreshold()) {
				return null;
			}

			scores.setRmsdCenters(subunitRmsd);
			Rotation symmetryOperation = createSymmetryOperation(permutation, transformation, axisAngle, fold, scores);
			return symmetryOperation;
		}
		return null;
	}

	/**
	 * Get valid rotation angles given the number of subunits
	 * @return The rotation angle corresponding to each fold of {@link Subunits#getFolds()}
	 */
	private List<Double> getAngles() {
		int n = subunits.getSubunitCount();
		// for spherical symmetric cases, n cannot be higher than 60
		if (n % 60 == 0 && isSpherical()) {
			 n = 60;
		}
		List<Integer> folds = subunits.getFolds();
		List<Double> angles = new ArrayList<Double>(folds.size()-1);

		// note this loop starts at 1, we do ignore 1-fold symmetry, which is the first entry
		for (int fold: folds) {
			if (fold > 0 && fold <= n) {
				angles.add(2* Math.PI/fold);
			}
		}
		return angles;
	}

	private boolean isSpherical() {
		MomentsOfInertia m = subunits.getMomentsOfInertia();
		return m.getSymmetryClass(0.05) == MomentsOfInertia.SymmetryClass.SYMMETRIC;
	}

	/**
	 * Checks if a particular permutation is allowed and superimposes well.
	 * Caches results.
	 * @param permutation
	 * @return null if invalid, or a rotation if valid
	 */
	private Rotation isValidPermutation(List<Integer> permutation) {
		if (permutation.size() == 0) {
			return null;
		}

		// cached value
		if (evaluatedPermutations.containsKey(permutation)) {
			return evaluatedPermutations.get(permutation);
		}

		// check if permutation is allowed
		if (! isAllowedPermutation(permutation)) {
			evaluatedPermutations.put(permutation, null);
			return null;
		}

		// check if superimposes
		Rotation rot = superimposePermutation(permutation);
		return rot;
	}

	/**
	 * The permutation must map all subunits onto an equivalent subunit
	 * and no subunit onto itself
	 * @param permutation
	 * @return
	 */
	private boolean isAllowedPermutation(List<Integer> permutation) {
		List<Integer> seqClusterId = subunits.getClusterIds();
		int selfaligned = 0;
		for (int i = 0; i < permutation.size(); i++) {
			int j = permutation.get(i);
			if ( seqClusterId.get(i) != seqClusterId.get(j)) {
				return false;
			}
			if(i == j ) {
				selfaligned++;
			}
		}
		// either identity (all self aligned) or all unique
		return selfaligned == 0 || selfaligned == permutation.size();
	}
	/**
	 * Adds translational component to rotation matrix
	 * @param rotation
	 * @return
	 */
	private void combineWithTranslation(Matrix4d rotation) {
		rotation.setTranslation(centroid);
		rotation.mul(rotation, centroidInverse);
	}

	private static Rotation createSymmetryOperation(List<Integer> permutation, Matrix4d transformation, AxisAngle4d axisAngle, int fold, QuatSymmetryScores scores) {
		Rotation s = new Rotation();
		s.setPermutation(new ArrayList<Integer>(permutation));
		s.setTransformation(new Matrix4d(transformation));
		s.setAxisAngle(new AxisAngle4d(axisAngle));
		s.setFold(fold);
		s.setScores(scores);

		return s;
	}

	private void setupDistanceBox() {
		distanceThreshold = calcDistanceThreshold();
		box = new DistanceBox<Integer>(distanceThreshold);

		for (int i = 0; i < originalCoords.length; i++) {
			box.addPoint(originalCoords[i], i);
		}
	}

	private double calcDistanceThreshold() {
		double threshold = Double.MAX_VALUE;
		int n = subunits.getSubunitCount();
		List<Point3d> centers = subunits.getCenters();

		for (int i = 0; i < n - 1; i++) {
			Point3d pi = centers.get(i);
			for (int j = i + 1; j < n; j++) {
				Point3d pj = centers.get(j);
				threshold = Math.min(threshold, pi.distanceSquared(pj));
			}
		}
		double distanceThreshold = Math.sqrt(threshold);

		distanceThreshold = Math.max(distanceThreshold, parameters.getRmsdThreshold());

		return distanceThreshold;
	}

	/**
	 * Compare this.transformedCoords with the original coords. For each
	 * subunit, return the transformed subunit with the closest position.
	 * @return A list mapping each subunit to the closest transformed subunit
	 */
	private List<Integer> getPermutation() {
		List<Integer> permutation = new ArrayList<Integer>(transformedCoords.length);
		double sum = 0.0f;

		for (Point3d t: transformedCoords) {
			List<Integer> neighbors = box.getNeighborsWithCache(t);
			int closest = -1;
			double minDist = Double.MAX_VALUE;

			 for (int j : neighbors) {
				double dist = t.distanceSquared(originalCoords[j]);
				if (dist < minDist) {
					closest = j;
					minDist = dist;
				}
			}

			sum += minDist;
			if (closest == -1) {
				 break;
			}
			permutation.add(closest);
		}
		double rmsd = Math.sqrt(sum / transformedCoords.length);

		if (rmsd > distanceThreshold || permutation.size() != transformedCoords.length) {
			permutation.clear();
			return permutation;
		}

		// check uniqueness of indices
		Set<Integer> set = new HashSet<Integer>(permutation);

		// if size mismatch, clear permutation (its invalid)
		if (set.size() != originalCoords.length) {
			permutation.clear();
		}

		return permutation;
	}

	private void initialize() {
		// translation to centered coordinate system
		centroid = new Vector3d(subunits.getCentroid());
		// translation back to original coordinate system
		Vector3d reverse = new Vector3d(centroid);
		reverse.negate();
		centroidInverse.set(reverse);
		// Make sure matrix element m33 is 1.0. An old version vecmath did not set this element.
		centroidInverse.setElement(3, 3,  1.0);

		List<Point3d> centers = subunits.getCenters();
		int n = subunits.getSubunitCount();

		originalCoords = new Point3d[n];
		transformedCoords = new Point3d[n];

		for (int i = 0; i < n; i++) {
			originalCoords[i] = centers.get(i);
			transformedCoords[i] = new Point3d();
		}

		setupDistanceBox();
	}
}
