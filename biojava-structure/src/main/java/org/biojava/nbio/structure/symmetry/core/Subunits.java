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

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Calc;
import org.biojava.nbio.structure.cluster.SubunitCluster;
import org.biojava.nbio.structure.symmetry.geometry.MomentsOfInertia;
import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

import java.util.*;

/**
 * A bean to represent information about the set of {@link Subunit} being
 * considered for symmetry detection. This class is a helper for the
 * {@link QuatSymmetryDetector} algorithm, since it calculates the
 * {@link MomentsOfInertia} and the centroids of each Subunit.
 * 
 * @author Peter Rose
 * @author Aleix Lafita
 * 
 */
public class Subunits {

	private List<Point3d[]> caCoords = new ArrayList<Point3d[]>();
	private List<Integer> sequenceClusterIds = new ArrayList<Integer>();

	private List<Integer> folds = new ArrayList<Integer>();
	private List<Point3d> originalCenters = new ArrayList<Point3d>();
	private List<Point3d> centers = new ArrayList<Point3d>();
	private List<Vector3d> unitVectors = new ArrayList<Vector3d>();

	private Point3d centroid;
	private MomentsOfInertia momentsOfInertia = new MomentsOfInertia();

	/**
	 * All input Lists should contain one element per subunit.
	 * 
	 * @param caCoords
	 *            CA coordinates of all subunits
	 * @param sequenceClusterIds
	 *            ID of the cluster that each subunit belongs to
	 * @param folds
	 *            Valid symmetry orders for this stoichiometry
	 */
	public Subunits(List<Point3d[]> caCoords, List<Integer> sequenceClusterIds,
			List<Integer> folds) {
		this.caCoords = caCoords;
		this.sequenceClusterIds = sequenceClusterIds;
		this.folds = folds;
	}

	/**
	 * Converts the List of {@link SubunitCluster} to a Subunit object.
	 * 
	 * @param clusters
	 *            List of SubunitCluster
	 */
	public Subunits(List<SubunitCluster> clusters) {

		// Loop through all subunits in the clusters and fill Lists
		for (int c = 0; c < clusters.size(); c++) {

			for (int s = 0; s < clusters.get(c).size(); s++) {
				sequenceClusterIds.add(c);

				Atom[] atoms = clusters.get(c).getAlignedAtomsSubunit(s);

				// Convert atoms to points
				Point3d[] points = new Point3d[atoms.length];
				for (int i = 0; i < atoms.length; i++)
					points[i] = new Point3d(atoms[i].getCoords());

				caCoords.add(points);
			}
		}

		// Fill in the folds with the function
		List<Integer> stoichiometry = new ArrayList<Integer>(clusters.size());
		for (int id = 0; id < clusters.size(); id++) {
			int size = clusters.get(id).size();
			stoichiometry.add(size);
		}
		folds = getValidFolds(stoichiometry);
	}

	public List<Point3d[]> getTraces() {
		return caCoords;
	}

	public int getSubunitCount() {
		run();
		if (centers == null) {
			return 0;
		}
		return centers.size();
	}

	public List<Integer> getSequenceClusterIds() {
		return sequenceClusterIds;
	}

	public List<Integer> getFolds() {
		return folds;
	}

	public String getStoichiometry() {

		// count number of members in each cluster
		Map<Integer, Integer> map = new TreeMap<Integer, Integer>();
		for (Integer id : sequenceClusterIds) {
			Integer value = map.get(id);
			if (value == null) {
				value = new Integer(1);
			} else {
				value++;
			}
			map.put(id, value);
		}

		List<Integer> stoichiometries = new ArrayList<Integer>(map.size());
		for (Integer key : map.keySet())
			stoichiometries.add(map.get(key));
		Collections.sort(stoichiometries);
		Collections.reverse(stoichiometries);

		// build formula string
		String alpha = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
		StringBuilder formula = new StringBuilder();
		for (int i = 0; i < stoichiometries.size(); i++) {
			String key = "?";
			if (i < alpha.length())
				key = alpha.substring(i, i + 1);

			formula.append(key);
			if (stoichiometries.get(i) > 1)
				formula.append(stoichiometries.get(i));
		}

		return formula.toString();
	}

	public int getCalphaCount() {
		int count = 0;
		for (Point3d[] trace : caCoords) {
			count += trace.length;
		}
		return count;
	}

	public int getLargestSubunit() {
		int index = -1;
		int maxLength = 0;
		for (int i = 0; i < caCoords.size(); i++) {
			int length = caCoords.get(i).length;
			if (length > maxLength) {
				index = i;
			}
		}
		return index;
	}

	public List<Point3d> getCenters() {
		run();
		return centers;
	}

	public List<Vector3d> getUnitVectors() {
		run();
		return unitVectors;
	}

	public List<Point3d> getOriginalCenters() {
		run();
		return originalCenters;
	}

	public Point3d getCentroid() {
		run();
		return centroid;
	}

	public MomentsOfInertia getMomentsOfInertia() {
		run();
		return momentsOfInertia;
	}

	private void run() {
		if (centers.size() > 0) {
			return;
		}
		calcOriginalCenters();
		calcCentroid();
		calcCenters();
		calcMomentsOfIntertia();
	}

	private void calcOriginalCenters() {
		for (Point3d[] trace : caCoords) {
			Point3d com = Calc.getCentroid(trace);
			originalCenters.add(com);
		}
	}

	private void calcCentroid() {
		Point3d[] orig = originalCenters.toArray(new Point3d[originalCenters
				.size()]);
		centroid = Calc.getCentroid(orig);
	}

	private void calcCenters() {
		for (Point3d p : originalCenters) {
			Point3d c = new Point3d(p);
			c.sub(centroid);
			centers.add(c);
			Vector3d v = new Vector3d(c);
			v.normalize();
			unitVectors.add(v);
		}
	}

	public Point3d getLowerBound() {
		Point3d lower = new Point3d();
		for (Point3d p : centers) {
			if (p.x < lower.x) {
				lower.x = p.x;
			}
			if (p.y < lower.y) {
				lower.y = p.y;
			}
			if (p.z < lower.z) {
				lower.z = p.z;
			}
		}
		return lower;
	}

	public Point3d getUpperBound() {
		Point3d upper = new Point3d();
		for (Point3d p : centers) {
			if (p.x > upper.x) {
				upper.x = p.x;
			}
			if (p.y > upper.y) {
				upper.y = p.y;
			}
			if (p.z > upper.z) {
				upper.z = p.z;
			}
		}
		return upper;
	}

	private void calcMomentsOfIntertia() {
		for (Point3d[] trace : caCoords) {
			for (Point3d p : trace) {
				momentsOfInertia.addPoint(p, 1.0f);
			}
		}
	}

	/**
	 * Find valid symmetry orders for a given stoichiometry. For instance, an
	 * A6B4 protein would give [1,2] because (A6B4)1 and (A3B2)2 are valid
	 * decompositions.
	 * 
	 * @param stoichiometry
	 *            List giving the number of copies in each chain cluster
	 * @return The common factors of the stoichiometry
	 */
	public static List<Integer> getValidFolds(List<Integer> stoichiometry) {

		List<Integer> denominators = new ArrayList<Integer>();

		if (stoichiometry.isEmpty())
			return denominators;

		int nChains = Collections.max(stoichiometry);

		// Remove duplicate stoichiometries
		Set<Integer> nominators = new TreeSet<Integer>(stoichiometry);

		// find common denominators
		for (int d = 1; d <= nChains; d++) {
			boolean isDivisable = true;
			for (Integer n : nominators) {
				if (n % d != 0) {
					isDivisable = false;
					break;
				}
			}
			if (isDivisable) {
				denominators.add(d);
			}
		}
		return denominators;
	}
}
