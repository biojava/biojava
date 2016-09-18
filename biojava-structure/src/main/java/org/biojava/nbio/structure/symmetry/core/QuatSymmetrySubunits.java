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
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.cluster.SubunitCluster;
import org.biojava.nbio.structure.geometry.CalcPoint;
import org.biojava.nbio.structure.geometry.MomentsOfInertia;
import org.biojava.nbio.structure.symmetry.utils.SymmetryTools;

import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

import java.util.*;
import java.util.stream.Collectors;

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
					points[i] = atoms[i].getCoordsAsPoint3d();

				caCoords.add(points);
			}
		}

		// List number of members in each cluster
		List<Integer> stoichiometries = clusters.stream().map(c -> c.size())
				.collect(Collectors.toList());
		folds = SymmetryTools.getValidFolds(stoichiometries);
	}

	public List<Point3d[]> getTraces() {
		return caCoords;
	}

	public List<Integer> getClusterIds() {
		return clusterIds;
	}

	/**
	 * This method is provisional and should only be used for coloring Subunits.
	 * A new coloring schema has to be implemented to allow the coloring of
	 * Subunits, without implying one Subunit = one Chain.
	 * 
	 * @return A List of the Chain Ids of each Subunit
	 */
	public List<String> getChainIds() {
		
		List<String> chains = new ArrayList<String>(getSubunitCount());

		// Loop through all subunits in the clusters and fill Lists
		for (int c = 0; c < clusters.size(); c++) {
			for (int s = 0; s < clusters.get(c).size(); s++)
				chains.add(clusters.get(c).getSubunits().get(s).getName());
		}
		
		return chains;
	}

	/**
	 * This method is provisional and should only be used for coloring Subunits.
	 * A new coloring schema has to be implemented to allow the coloring of
	 * Subunits, without implying one Subunit = one Chain.
	 * 
	 * @return A List of the Model number of each Subunit
	 */
	public List<Integer> getModelNumbers() {
		
		List<Integer> models = new ArrayList<Integer>(getSubunitCount());

		// Loop through all subunits in the clusters and fill Lists
		for (int c = 0; c < clusters.size(); c++) {
			for (int s = 0; s < clusters.get(c).size(); s++) {

				Atom[] atoms = clusters.get(c).getAlignedAtomsSubunit(s);

				// TODO guess them chain and model (very ugly)
				Chain chain = atoms[0].getGroup().getChain();

				int model = 0;
				for (int m = 0; m < chain.getStructure().nrModels(); m++) {
					if (chain.getStructure().getModel(m).contains(chain)) {
						model = m;
						break;
					}
				}
				models.add(model);
			}
		}
		return models;
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
			Point3d com = CalcPoint.centroid(trace);
			originalCenters.add(com);
		}
	}

	private void calcCentroid() {
		Point3d[] orig = originalCenters.toArray(new Point3d[originalCenters
				.size()]);
		centroid = CalcPoint.centroid(orig);
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
