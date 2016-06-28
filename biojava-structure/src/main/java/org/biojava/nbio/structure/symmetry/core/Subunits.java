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
import org.biojava.nbio.structure.cluster.SubunitClustererMethod;
import org.biojava.nbio.structure.symmetry.geometry.MomentsOfInertia;
import org.biojava.nbio.structure.symmetry.geometry.SuperPosition;

import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

import java.util.*;

/**
 * A bean to represent information about the set of {@link Subunit} being
 * considered for symmetry detection. This class is a helper for the
 * {@link QuatSymmetryDetector} algorithm, since it calculates the
 * {@link MomentsOfInertia} and the centroids of each Subunit.
 * <p>
 * Delete public modifier when finished.
 * 
 * @author Peter Rose
 * @author Aleix Lafita
 * 
 */
public class Subunits {
	
	private List<Point3d[]> caCoords = new ArrayList<Point3d[]>();
	private List<Integer> sequenceClusterIds = new ArrayList<Integer>();

	@Deprecated
	// Return variable, guessed from clusters
	private List<Boolean> pseudoStoichiometry = new ArrayList<Boolean>();
	@Deprecated
	// Return variable, no longer valid
	private List<Double> minSequenceIdentity = new ArrayList<Double>();
	@Deprecated
	// Return variable, no longer valid
	private List<Double> maxSequenceIdentity = new ArrayList<Double>();

	@Deprecated
	// Subunits are not tied to a chain, used for coloring
	private List<String> chainIds = new ArrayList<String>();
	@Deprecated
	// Subunits are not tied to a model? used for coloring
	private List<Integer> modelNumbers = new ArrayList<Integer>();

	private List<Integer> folds = new ArrayList<Integer>();
	private List<Point3d> originalCenters = new ArrayList<Point3d>();
	private List<Point3d> centers = new ArrayList<Point3d>();
	private List<Vector3d> unitVectors = new ArrayList<Vector3d>();

	@Deprecated
	// Return variable, should be in the clusters
	private int nucleicAcidChainCount = 0;
	@Deprecated
	// This should be in QuatSymmetryResults
	private boolean pseudoSymmetric = false;

	private Point3d centroid;
	private MomentsOfInertia momentsOfInertia = new MomentsOfInertia();

	/**
	 * All input Lists should contain one element per subunit.
	 * 
	 * @param caCoords
	 *            CA coordinates of all subunits
	 * @param sequenceClusterIds
	 *            ID of the cluster that each subunit belongs to
	 * @param pseudoStoichiometry
	 *            Whether pseudosymmetry was used when clustering the subunit
	 * @param minSequenceIdentity
	 *            Minimum sequence identity to other cluster members
	 * @param maxSequenceIdentity
	 *            Maximum sequence identity to other cluster members
	 * @param folds
	 *            Valid symmetry orders for this stoichiometry
	 * @param chainIds
	 *            Chain ID for the subunit
	 * @param modelNumbers
	 *            Model number for the subunit
	 */
	@Deprecated
	public Subunits(List<Point3d[]> caCoords, List<Integer> sequenceClusterIds,
			List<Boolean> pseudoStoichiometry,
			List<Double> minSequenceIdentity, List<Double> maxSequenceIdentity,
			List<Integer> folds, List<String> chainIds,
			List<Integer> modelNumbers) {
		this.caCoords = caCoords;
		this.sequenceClusterIds = sequenceClusterIds;
		this.pseudoStoichiometry = pseudoStoichiometry;
		this.minSequenceIdentity = minSequenceIdentity;
		this.maxSequenceIdentity = maxSequenceIdentity;
		this.folds = folds;
		this.chainIds = chainIds;
		this.modelNumbers = modelNumbers;
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

			// TODO we should remove these variables
			minSequenceIdentity.add(0.0);
			maxSequenceIdentity.add(0.0);

			// Pseudostoichiometry means one structural cluster
			SubunitClustererMethod method = clusters.get(c)
					.getClustererMethod();
			boolean ps = (method == SubunitClustererMethod.STRUCTURE);
			ps = (ps || method == SubunitClustererMethod.INTERNAL_SYMMETRY);

			for (int s = 0; s < clusters.get(c).size(); s++) {
				sequenceClusterIds.add(c);
				pseudoStoichiometry.add(ps);

				Atom[] atoms = clusters.get(c).getAlignedAtomsSubunit(s);

				// Convert atoms to points
				Point3d[] points = new Point3d[atoms.length];
				for (int i = 0; i < atoms.length; i++)
					points[i] = new Point3d(atoms[i].getCoords());

				caCoords.add(points);

				// TODO guess them chain and model (very ugly)
				Chain chain = atoms[0].getGroup().getChain();
				String cid = chain.getId();
				chainIds.add(cid);
				
				int model = 0;
				for (int m = 0; m < chain.getStructure().nrModels(); m++){
					if (chain.getStructure().getModel(m).contains(chain)) {
						model = m;
						break;
					}
				}				
				modelNumbers.add(model);
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

	public boolean isPseudoStoichiometric() {
		for (Boolean b : pseudoStoichiometry) {
			if (b) {
				return true;
			}
		}
		return false;
	}

	public boolean isPseudoSymmetric() {
		return pseudoSymmetric;
	}

	public void setPseudoSymmetric(boolean pseudoSymmetric) {
		this.pseudoSymmetric = pseudoSymmetric;
	}

	public double getMinSequenceIdentity() {
		double minId = 1.0;
		for (double seqId : minSequenceIdentity) {
			minId = Math.min(seqId, minId);
		}
		return minId;
	}

	public double getMaxSequenceIdentity() {
		double maxId = 1.0;
		for (double seqId : maxSequenceIdentity) {
			maxId = Math.min(seqId, maxId);
		}
		return maxId;
	}

	public List<String> getChainIds() {
		return chainIds;
	}

	public List<Integer> getModelNumbers() {
		return modelNumbers;
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

		// build formula string
		String alpha = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
		StringBuilder formula = new StringBuilder();
		for (Integer stoich : stoichiometries) {
			String key = "?";
			if (stoich < alpha.length()) 
				key = alpha.substring(stoich, stoich + 1);
			
			formula.append(key);
			formula.append(stoich);
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

	/**
	 * @return the nucleicAcidChainCount
	 */
	public int getNucleicAcidChainCount() {
		run();
		return nucleicAcidChainCount;
	}

	/**
	 * @param nucleicAcidChainCount
	 *            the nucleicAcidChainCount to set
	 */
	public void setNucleicAcidChainCount(int nucleicAcidChainCount) {
		this.nucleicAcidChainCount = nucleicAcidChainCount;
	}

	public boolean overlaps(Subunits subunits) {
		Set<String> set1 = getSignatures(this);
		Set<String> set2 = getSignatures(subunits);
		set1.retainAll(set2);
		return set1.size() > 0;
	}

	public boolean contains(Subunits subunits) {
		Set<String> set1 = getSignatures(this);
		Set<String> set2 = getSignatures(subunits);
		return set1.containsAll(set2);
	}

	private static Set<String> getSignatures(Subunits subunits) {
		Set<String> set = new HashSet<String>(subunits.getSubunitCount());
		for (int i = 0; i < subunits.getSubunitCount(); i++) {
			set.add(subunits.getChainIds().get(i) + "_"
					+ subunits.getModelNumbers().get(i));
		}
		return set;
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
			Point3d com = SuperPosition.centroid(trace);
			originalCenters.add(com);
		}
	}

	private void calcCentroid() {
		Point3d[] orig = originalCenters.toArray(new Point3d[originalCenters
				.size()]);
		centroid = SuperPosition.centroid(orig);
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
