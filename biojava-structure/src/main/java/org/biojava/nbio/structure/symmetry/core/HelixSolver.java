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
import org.biojava.nbio.structure.geometry.SuperPositions;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import javax.vecmath.AxisAngle4d;
import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

import java.util.*;
import java.util.Map.Entry;

/**
 * 
 * 
 * @author Peter Rose
 *
 */
public class HelixSolver {

	private static final Logger logger = LoggerFactory
			.getLogger(HelixSolver.class);

	private QuatSymmetrySubunits subunits = null;
	private int fold = 1;
	private HelixLayers helixLayers = new HelixLayers();
	private QuatSymmetryParameters parameters = null;
	boolean modified = true;

	public HelixSolver(QuatSymmetrySubunits subunits, int fold,
			QuatSymmetryParameters parameters) {
		this.subunits = subunits;
		this.fold = fold;
		this.parameters = parameters;
	}

	public HelixLayers getSymmetryOperations() {
		if (modified) {
			solve();
			modified = false;
		}
		return helixLayers;
	}

	private void solve() {
		if (!preCheck()) {
			return;
		}

		HelicalRepeatUnit unit = new HelicalRepeatUnit(subunits);
		List<Point3d> repeatUnitCenters = unit.getRepeatUnitCenters();
		List<Point3d[]> repeatUnits = unit.getRepeatUnits();
		Set<List<Integer>> permutations = new HashSet<List<Integer>>();

		double minRise = parameters.getMinimumHelixRise() * fold; // for n-start
																	// helix,
																	// the rise
																	// must be
																	// steeper
		Map<Integer[], Integer> interactionMap = unit
				.getInteractingRepeatUnits();

		int maxLayerLineLength = 0;

		for (Entry<Integer[], Integer> entry : interactionMap.entrySet()) {
			Integer[] pair = entry.getKey();
			logger.debug("HelixSolver: pair: " + Arrays.toString(pair));
			
			int contacts = entry.getValue();
			Point3d[] h1 = CalcPoint.clonePoint3dArray(repeatUnits.get(pair[0]));
			Point3d[] h2 = CalcPoint.clonePoint3dArray(repeatUnits.get(pair[1]));

			// trial superposition of repeat unit pairs to get a seed
			// permutation
			Matrix4d transformation = SuperPositions.superposeAndTransform(h2, h1);

			double rmsd = CalcPoint.rmsd(h1, h2);
			double rise = getRise(transformation,
					repeatUnitCenters.get(pair[0]),
					repeatUnitCenters.get(pair[1]));
			double angle = getAngle(transformation);

			logger.debug(
					"Original rmsd: {}, Original rise {}, Original angle: {}",
					rmsd, rise, Math.toDegrees(angle));

			if (rmsd > parameters.getRmsdThreshold()) {
				continue;
			}

			if (Math.abs(rise) < minRise) {
				continue;
			}

			// determine which subunits are permuted by the transformation
			List<Integer> permutation = getPermutation(transformation);

			// check permutations for validity

			// don't save redundant permutations
			if (permutations.contains(permutation)) {
				continue;
			}
			permutations.add(permutation);
			logger.debug("Permutation: " + permutation);
			

			// keep track of which subunits are permuted
			Set<Integer> permSet = new HashSet<Integer>();
			int count = 0;
			boolean valid = true;
			for (int i = 0; i < permutation.size(); i++) {
				if (permutation.get(i) == i) {
					valid = false;
					break;
				}
				if (permutation.get(i) != -1) {
					permSet.add(permutation.get(i));
					permSet.add(i);
					count++;
				}

			}

			// a helix a repeat unit cannot map onto itself
			if (!valid) {
				logger.debug("Invalid mapping");
				continue;
			}

			// all subunits must be involved in a permutation
			if (permSet.size() != subunits.getSubunitCount()) {
				logger.debug("Not all subunits involved in permutation");
				continue;
			}

			// if all subunit permutation values are set, then it can't be
			// helical symmetry (must be cyclic symmetry)
			if (count == permutation.size()) {
				continue;
			}

			// superpose all permuted subunits
			List<Point3d> point1 = new ArrayList<Point3d>();
			List<Point3d> point2 = new ArrayList<Point3d>();
			List<Point3d> centers = subunits.getOriginalCenters();
			for (int j = 0; j < permutation.size(); j++) {
				if (permutation.get(j) != -1) {
					point1.add(new Point3d(centers.get(j)));
					point2.add(new Point3d(centers.get(permutation.get(j))));
				}
			}

			h1 = new Point3d[point1.size()];
			h2 = new Point3d[point2.size()];
			point1.toArray(h1);
			point2.toArray(h2);

			// calculate subunit rmsd if at least 3 subunits are available
			double subunitRmsd = 0;
			if (point1.size() > 2) {
				transformation = SuperPositions.superposeAndTransform(h2, h1);

				subunitRmsd = CalcPoint.rmsd(h1, h2);
				rise = getRise(transformation, repeatUnitCenters.get(pair[0]),
						repeatUnitCenters.get(pair[1]));
				angle = getAngle(transformation);

				logger.debug("Subunit rmsd: {}, Subunit rise: {}, Subunit angle: {}", subunitRmsd, rise, Math.toDegrees(angle));

				if (subunitRmsd > parameters.getRmsdThreshold()) {
					continue;
				}

				if (Math.abs(rise) < minRise) {
					continue;
				}

				if (subunitRmsd > parameters.getHelixRmsdToRiseRatio()
						* Math.abs(rise)) {
					continue;
				}
			}

			// superpose all C alpha traces
			point1.clear();
			point2.clear();
			List<Point3d[]> traces = subunits.getTraces();
			for (int j = 0; j < permutation.size(); j++) {
				if (permutation.get(j) == -1) {
					continue;
				}
				for (Point3d p : traces.get(j)) {
					point1.add(new Point3d(p));
				}
				for (Point3d p : traces.get(permutation.get(j))) {
					point2.add(new Point3d(p));
				}
			}

			h1 = new Point3d[point1.size()];
			h2 = new Point3d[point2.size()];
			point1.toArray(h1);
			point2.toArray(h2);
			Point3d[] h3 = CalcPoint.clonePoint3dArray(h1);
			transformation = SuperPositions.superposeAndTransform(h2, h1);

			Point3d xtrans = CalcPoint.centroid(h3);

			xtrans.negate();

			double traceRmsd = CalcPoint.rmsd(h1, h2);

			rise = getRise(transformation, repeatUnitCenters.get(pair[0]),
					repeatUnitCenters.get(pair[1]));
			angle = getAngle(transformation);

			logger.debug("Trace rmsd: " + traceRmsd);
			logger.debug("Trace rise: " + rise);
			logger.debug("Trace angle: " + Math.toDegrees(angle));
			logger.debug("Permutation: " + permutation);

			if (traceRmsd > parameters.getRmsdThreshold()) {
				continue;
			}

			if (Math.abs(rise) < minRise) {
				continue;
			}

			// This prevents translational repeats to be counted as helices
			if (angle < Math.toRadians(parameters.getMinimumHelixAngle())) {
				continue;
			}

			if (traceRmsd > parameters.getHelixRmsdToRiseRatio()
					* Math.abs(rise)) {
				continue;
			}

			AxisAngle4d a1 = new AxisAngle4d();
			a1.set(transformation);

			// save this helix rot-translation
			Helix helix = new Helix();
			helix.setTransformation(transformation);
			helix.setPermutation(permutation);
			helix.setRise(rise);
			// Old version of Vecmath on LINUX doesn't set element m33 to 1.
			// Here we make sure it's 1.
			transformation.setElement(3, 3, 1.0);
			transformation.invert();
			QuatSymmetryScores scores = QuatSuperpositionScorer.calcScores(
					subunits, transformation, permutation);
			scores.setRmsdCenters(subunitRmsd);
			helix.setScores(scores);
			helix.setFold(fold);
			helix.setContacts(contacts);
			helix.setRepeatUnits(unit.getRepeatUnitIndices());
			logger.debug("Layerlines: " + helix.getLayerLines());
			
			for (List<Integer> line : helix.getLayerLines()) {
				maxLayerLineLength = Math.max(maxLayerLineLength, line.size());
			}

			// TODO
			// checkSelfLimitingHelix(helix);

			helixLayers.addHelix(helix);

		}
		if (maxLayerLineLength < 3) {
			// System.out.println("maxLayerLineLength: " + maxLayerLineLength);
			helixLayers.clear();
		}

		return;
	}

	@SuppressWarnings("unused")
	private void checkSelfLimitingHelix(Helix helix) {
		HelixExtender he = new HelixExtender(subunits, helix);
		Point3d[] extendedHelix = he.extendHelix(1);

		int overlap1 = 0;
		for (Point3d[] trace : subunits.getTraces()) {
			for (Point3d pt : trace) {
				for (Point3d pe : extendedHelix) {
					if (pt.distance(pe) < 5.0) {
						overlap1++;
					}
				}
			}
		}

		extendedHelix = he.extendHelix(-1);

		int overlap2 = 0;
		for (Point3d[] trace : subunits.getTraces()) {
			for (Point3d pt : trace) {
				for (Point3d pe : extendedHelix) {
					if (pt.distance(pe) < 3.0) {
						overlap2++;
					}
				}
			}
		}
		System.out.println("SelfLimiting helix: " + overlap1 + ", " + overlap2);
	}

	private boolean preCheck() {
		if (subunits.getSubunitCount() < 3) {
			return false;
		}
		List<Integer> folds = this.subunits.getFolds();
		int maxFold = folds.get(folds.size() - 1);
		return maxFold > 1;
	}

	/**
	 * Returns a permutation of subunit indices for the given helix
	 * transformation. An index of -1 is used to indicate subunits that do not
	 * superpose onto any other subunit.
	 * 
	 * @param transformation
	 * @return
	 */
	private List<Integer> getPermutation(Matrix4d transformation) {
		double rmsdThresholdSq = Math
				.pow(this.parameters.getRmsdThreshold(), 2);

		List<Point3d> centers = subunits.getOriginalCenters();
		List<Integer> seqClusterId = subunits.getClusterIds();

		List<Integer> permutations = new ArrayList<Integer>(centers.size());
		double[] dSqs = new double[centers.size()];
		boolean[] used = new boolean[centers.size()];
		Arrays.fill(used, false);

		for (int i = 0; i < centers.size(); i++) {
			Point3d tCenter = new Point3d(centers.get(i));
			transformation.transform(tCenter);
			int permutation = -1;
			double minDistSq = Double.MAX_VALUE;
			for (int j = 0; j < centers.size(); j++) {
				if (seqClusterId.get(i) == seqClusterId.get(j)) {
					if (!used[j]) {
						double dSq = tCenter.distanceSquared(centers.get(j));
						if (dSq < minDistSq && dSq <= rmsdThresholdSq) {
							minDistSq = dSq;
							permutation = j;
							dSqs[j] = dSq;
						}
					}
				}
			}
			// can't map to itself
			if (permutations.size() == permutation) {
				permutation = -1;
			}

			if (permutation != -1) {
				used[permutation] = true;
			}

			permutations.add(permutation);
		}

		return permutations;
	}

	/**
	 * Returns the rise of a helix given the subunit centers of two adjacent
	 * subunits and the helix transformation
	 * 
	 * @param transformation
	 *            helix transformation
	 * @param p1
	 *            center of one subunit
	 * @param p2
	 *            center of an adjacent subunit
	 * @return
	 */
	private static double getRise(Matrix4d transformation, Point3d p1,
			Point3d p2) {
		AxisAngle4d axis = getAxisAngle(transformation);
		Vector3d h = new Vector3d(axis.x, axis.y, axis.z);
		Vector3d p = new Vector3d();
		p.sub(p1, p2);
		return p.dot(h);
	}

	/**
	 * Returns the pitch angle of the helix
	 * 
	 * @param transformation
	 *            helix transformation
	 * @return
	 */
	private static double getAngle(Matrix4d transformation) {
		return getAxisAngle(transformation).angle;
	}

	/**
	 * Returns the AxisAngle of the helix transformation
	 * 
	 * @param transformation
	 *            helix transformation
	 * @return
	 */
	private static AxisAngle4d getAxisAngle(Matrix4d transformation) {
		AxisAngle4d axis = new AxisAngle4d();
		axis.set(transformation);
		return axis;
	}
}
