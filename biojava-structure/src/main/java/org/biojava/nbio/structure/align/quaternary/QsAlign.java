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
package org.biojava.nbio.structure.align.quaternary;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import javax.vecmath.Matrix4d;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Calc;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.multiple.Block;
import org.biojava.nbio.structure.align.multiple.BlockImpl;
import org.biojava.nbio.structure.align.multiple.BlockSet;
import org.biojava.nbio.structure.align.multiple.BlockSetImpl;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.align.multiple.MultipleAlignmentEnsembleImpl;
import org.biojava.nbio.structure.align.multiple.MultipleAlignmentImpl;
import org.biojava.nbio.structure.align.multiple.util.MultipleAlignmentScorer;
import org.biojava.nbio.structure.align.multiple.util.ReferenceSuperimposer;
import org.biojava.nbio.structure.cluster.Subunit;
import org.biojava.nbio.structure.cluster.SubunitCluster;
import org.biojava.nbio.structure.cluster.SubunitClusterer;
import org.biojava.nbio.structure.cluster.SubunitClustererParameters;
import org.biojava.nbio.structure.cluster.SubunitExtractor;
import org.biojava.nbio.structure.contact.Pair;
import org.biojava.nbio.structure.geometry.SuperPositions;
import org.biojava.nbio.structure.geometry.UnitQuaternions;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Quaternary Structure Alignment (QS-Align). The algorithm takes as input two
 * protein structures at the quaternary structure level (multiple chains or
 * subunits) and calculates the equivalent subunit matching and a residue-based
 * alignment, together with usual alignment quality scores.
 * 
 * @author Aleix Lafita
 * @since 5.0.0
 *
 */
public class QsAlign {

	private static final Logger logger = LoggerFactory.getLogger(QsAlign.class);

	public static QsAlignResult align(Structure s1, Structure s2,
			SubunitClustererParameters cParams, QsAlignParameters aParams)
			throws StructureException {
		return align(
				SubunitExtractor.extractSubunits(s1,
						cParams.getAbsoluteMinimumSequenceLength(),
						cParams.getMinimumSequenceLengthFraction(),
						cParams.getMinimumSequenceLength()),
				SubunitExtractor.extractSubunits(s2,
						cParams.getAbsoluteMinimumSequenceLength(),
						cParams.getMinimumSequenceLengthFraction(),
						cParams.getMinimumSequenceLength()), cParams, aParams);
	}

	public static QsAlignResult align(List<Subunit> s1, List<Subunit> s2,
			SubunitClustererParameters cParams, QsAlignParameters aParams)
			throws StructureException {

		QsAlignResult result = new QsAlignResult(s1, s2);

		// SETP 1: cluster each group of subunits O(N^2*L^2) - intra

		List<SubunitCluster> c1 = SubunitClusterer.cluster(s1, cParams).getClusters();
		List<SubunitCluster> c2 = SubunitClusterer.cluster(s2, cParams).getClusters();

		// STEP 2: match each subunit cluster between groups O(N^2*L^2) - inter
		Map<Integer, Integer> clusterMap = new HashMap<Integer, Integer>();
		for (int i = 0; i < c1.size(); i++) {
			for (int j = 0; j < c2.size(); j++) {

				if (clusterMap.keySet().contains(i))
					break;
				if (clusterMap.values().contains(j))
					continue;

				// Use structural alignment to match the subunit clusters
				if (c1.get(i).mergeStructure(c2.get(j),cParams)) {
					clusterMap.put(i, j);
				}
			}
		}

		logger.info("Cluster Map: " + clusterMap.toString());
		result.setClusters(c1);

		// STEP 3: Align the assemblies for each cluster match O(N^2*L)
		for (int globalKey : clusterMap.keySet()) {

			// Obtain the clusters
			SubunitCluster clust1 = c1.get(globalKey);
			SubunitCluster clust2 = c2.get(clusterMap.get(globalKey));

			// Take a cluster match as reference
			int index1 = 0;
			int index2 = clust1.size() - clust2.size();
			Map<Integer, Integer> subunitMap = new HashMap<Integer, Integer>();
			subunitMap.put(index1, index2);

			// Map cluster id to their subunit matching
			Map<Integer, Map<Integer, Integer>> clustSubunitMap = new HashMap<Integer, Map<Integer, Integer>>();
			clustSubunitMap.put(globalKey, subunitMap);

			// Change order of key set so that globalKey is first
			List<Integer> keySet = new ArrayList<Integer>(clusterMap.keySet());
			keySet.remove((Integer) globalKey);
			keySet.add(0, globalKey);

			for (int key : clusterMap.keySet()) {

				// Recover subunitMap if it is the reference, new one otherwise
				if (key == globalKey)
					subunitMap = clustSubunitMap.get(key);
				else
					subunitMap = new HashMap<Integer, Integer>();

				// Obtain the clusters of each subunit group
				clust1 = c1.get(key);
				clust2 = c2.get(clusterMap.get(key));

				// Get the initial subunit indices of each group
				index1 = 0;
				index2 = clust1.size() - clust2.size();

				for (int i = 0; i < index2; i++) {
					for (int j = index2; j < clust1.size(); j++) {

						if (subunitMap.keySet().contains(i))
							break;
						if (subunitMap.values().contains(j))
							continue;

						// Obtain cumulative transformation matrix
						Matrix4d transform = getTransformForClusterSubunitMap(
								c1, clustSubunitMap);

						// Obtain Atom arrays of the subunit pair to match
						Atom[] atoms1 = clust1.getAlignedAtomsSubunit(i);
						Atom[] atoms2 = clust1.getAlignedAtomsSubunit(j);

						// Obtain centroids and transform the second
						Atom centr1 = Calc.getCentroid(atoms1);
						Atom centr2 = Calc.getCentroid(atoms2);
						Calc.transform(centr2, transform);

						// 1- Check that the distance fulfills maximum
						double dCentroid = Calc.getDistance(centr1, centr2);
						if (dCentroid > aParams.getdCutoff()) {
							logger.debug(String.format("Subunit matching %d "
									+ "vs %d of cluster %d could not be "
									+ "matched, because centroid distance is "
									+ "%.2f", index1, index2, key, dCentroid));
							continue;
						}

						// Transform coordinates of second
						Atom[] atoms2c = StructureTools.cloneAtomArray(atoms2);
						Calc.transform(atoms2c, transform);

						// 2- Check the orientation metric condition
						double qOrient = UnitQuaternions.orientationAngle(
								Calc.atomsToPoints(atoms1),
								Calc.atomsToPoints(atoms2c), false);
						qOrient = Math.min(Math.abs(2*Math.PI - qOrient), qOrient);
						if (qOrient > aParams.getMaxOrientationAngle()) {
							logger.debug(String.format("Subunit matching %d "
									+ "vs %d of cluster %d could not be "
									+ "matched, because orientation metric is "
									+ "%.2f", i, j, key, qOrient));
							continue;
						}

						// 3- Check the RMSD condition
						double rmsd = Calc.rmsd(atoms1, atoms2c);
						if (rmsd > aParams.getMaxRmsd()) {
							logger.debug(String.format("Subunit matching %d "
									+ "vs %d of cluster %d could not be "
									+ "matched, because RMSD is %.2f", i,
									j, key, rmsd));
							continue;
						}

						logger.info(String.format("Subunit matching %d vs %d"
								+ " of cluster %d with centroid distance %.2f"
								+ ", orientation metric %.2f and RMSD %.2f",
								i, j, key, dCentroid, qOrient, rmsd));

						subunitMap.put(i, j);
					}
				}

				clustSubunitMap.put(key, subunitMap);
			}
			
			logger.info("Cluster Subunit Map: " + clustSubunitMap.toString());

			// Unfold the nested map into subunit map and alignment
			subunitMap = new HashMap<Integer, Integer>();
			List<Integer> alignRes1 = new ArrayList<Integer>();
			List<Integer> alignRes2 = new ArrayList<Integer>();
			List<Atom> atomArray1 = new ArrayList<Atom>();
			List<Atom> atomArray2 = new ArrayList<Atom>();

			for (int key : clustSubunitMap.keySet()) {

				// Obtain the cluster and the alignment in it
				SubunitCluster cluster = c1.get(key);
				List<List<Integer>> clusterEqrs = cluster
						.getMultipleAlignment().getBlock(0).getAlignRes();

				for (Entry<Integer, Integer> pair : clustSubunitMap.get(key)
						.entrySet()) {

					int i = pair.getKey();
					int j = pair.getValue();

					// Obtain the indices of the original Subunit Lists
					int orig1 = s1.indexOf(cluster.getSubunits().get(i));
					int orig2 = s2.indexOf(cluster.getSubunits().get(j));

					// Append rescaled aligned residue indices
					for (Integer eqr : clusterEqrs.get(i))
						alignRes1.add(eqr + atomArray1.size());
					for (Integer eqr : clusterEqrs.get(j))
						alignRes2.add(eqr + atomArray2.size());

					// Apend atoms to the arrays
					atomArray1.addAll(Arrays.asList(s1.get(orig1)
							.getRepresentativeAtoms()));
					atomArray2.addAll(Arrays.asList(s2.get(orig2)
							.getRepresentativeAtoms()));

					subunitMap.put(orig1, orig2);
				}
			}

			// Evaluate the goodness of the match with an alignment object
			MultipleAlignment msa = new MultipleAlignmentImpl();
			msa.setEnsemble(new MultipleAlignmentEnsembleImpl());
			msa.getEnsemble().setAtomArrays(
					Arrays.asList(new Atom[][] {
							atomArray1.toArray(new Atom[atomArray1.size()]),
							atomArray2.toArray(new Atom[atomArray2.size()]) }));

			// Fill in the alignment information
			BlockSet bs = new BlockSetImpl(msa);
			Block b = new BlockImpl(bs);
			List<List<Integer>> alignRes = new ArrayList<List<Integer>>(2);
			alignRes.add(alignRes1);
			alignRes.add(alignRes2);
			b.setAlignRes(alignRes);

			// Fill in the transformation matrices
			new ReferenceSuperimposer().superimpose(msa);

			// Calculate some scores
			MultipleAlignmentScorer.calculateScores(msa);

			// If it is the best match found so far store it
			if (subunitMap.size() > result.getSubunitMap().size()) {
				result.setSubunitMap(subunitMap);
				result.setAlignment(msa);
				logger.info("Better result found: " + result.toString());
			} else if (subunitMap.size() == result.getSubunitMap().size()) {
				if (result.getAlignment() == null) {
					result.setSubunitMap(subunitMap);
					result.setAlignment(msa);
				} else if (msa.getScore(MultipleAlignmentScorer.RMSD) < result
						.getRmsd()) {
					result.setSubunitMap(subunitMap);
					result.setAlignment(msa);
					logger.info("Better result found: " + result.toString());
				}
			}

		}

		return result;
	}

	/**
	 * Returns a pair of Atom arrays corresponding to the alignment of subunit
	 * matchings, in order of appearance. Superposition of the two Atom sets
	 * gives the transformation of the complex.
	 * <p>
	 * Utility method to cumulative calculate the alignment Atoms.
	 * 
	 * @param clusters
	 *            List of SubunitClusters
	 * @param clusterSubunitMap
	 *            map from cluster id to subunit matching
	 * @return pair of atom arrays to be superposed
	 */
	private static Pair<Atom[]> getAlignedAtomsForClusterSubunitMap(
			List<SubunitCluster> clusters,
			Map<Integer, Map<Integer, Integer>> clusterSubunitMap) {

		List<Atom> atomArray1 = new ArrayList<Atom>();
		List<Atom> atomArray2 = new ArrayList<Atom>();

		// For each cluster of subunits
		for (int key : clusterSubunitMap.keySet()) {

			// Obtain the cluster and the alignment in it
			SubunitCluster cluster = clusters.get(key);

			// For each subunit matching in the cluster
			for (Entry<Integer, Integer> pair : clusterSubunitMap.get(key)
					.entrySet()) {

				int i = pair.getKey();
				int j = pair.getValue();

				// Apend atoms to the arrays
				atomArray1.addAll(Arrays.asList(cluster
						.getAlignedAtomsSubunit(i)));
				atomArray2.addAll(Arrays.asList(cluster
						.getAlignedAtomsSubunit(j)));
			}

		}
		return new Pair<Atom[]>(
				atomArray1.toArray(new Atom[atomArray1.size()]),
				atomArray2.toArray(new Atom[atomArray2.size()]));
	}

	/**
	 * Returns the transformation matrix corresponding to the alignment of
	 * subunit matchings.
	 * <p>
	 * Utility method to cumulative calculate the alignment transformation.
	 * 
	 * @param clusters
	 *            List of SubunitClusters
	 * @param clusterSubunitMap
	 *            map from cluster id to subunit matching
	 * @return transformation matrix
	 * @throws StructureException
	 */
	private static Matrix4d getTransformForClusterSubunitMap(
			List<SubunitCluster> clusters,
			Map<Integer, Map<Integer, Integer>> clusterSubunitMap)
			throws StructureException {

		Pair<Atom[]> pair = getAlignedAtomsForClusterSubunitMap(clusters,
				clusterSubunitMap);

		return SuperPositions.superpose(Calc.atomsToPoints(pair.getFirst()), 
				Calc.atomsToPoints(pair.getSecond()));

	}
}
