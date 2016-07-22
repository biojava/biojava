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
 * {@link QuatSymmetryDetector} algorithm, since it calculates and caches the
 * {@link MomentsOfInertia} and the centroids of each Subunit.
 * 
 * @author Peter Rose
 * @author Aleix Lafita
 * 
 */
public class QuatSymmetrySubunits {

	private List<Point3d[]> caCoords = new ArrayList<Point3d[]>();
	private List<Point3d> originalCenters = new ArrayList<Point3d>();
	private List<Point3d> centers = new ArrayList<Point3d>();
	private List<Vector3d> unitVectors = new ArrayList<Vector3d>();

	private List<Integer> folds = new ArrayList<Integer>();
	private List<Integer> clusterIds = new ArrayList<Integer>();
	private List<SubunitCluster> clusters;

	private Point3d centroid;
	private MomentsOfInertia momentsOfInertia = new MomentsOfInertia();

	/**
	 * Converts the List of {@link SubunitCluster} to a Subunit object.
	 * 
	 * @param clusters
	 *            List of SubunitCluster
	 */
	public QuatSymmetrySubunits(List<SubunitCluster> clusters) {

		this.clusters = clusters;

		// Loop through all subunits in the clusters and fill Lists
		for (int c = 0; c < clusters.size(); c++) {

			for (int s = 0; s < clusters.get(c).size(); s++) {

				clusterIds.add(c);
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

	public List<Integer> getClusterIds() {
		return clusterIds;
	}

	public List<SubunitCluster> getClusters() {
		return clusters;
	}

	public int getSubunitCount() {
		run();
		if (centers == null) {
			return 0;
		}
		return centers.size();
	}

	public List<Integer> getFolds() {
		return folds;
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

}
