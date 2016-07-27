package org.biojava.nbio.structure.align.quaternary;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import javax.vecmath.Matrix4d;

import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Calc;
import org.biojava.nbio.structure.SVDSuperimposer;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
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
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Quaternary Structure Alignment (QS-Align). The algorithm takes as input two
 * protein structures at the quaternary structure level (multiple interacting
 * chains) and calculates the equivalent cross chains and the optimal
 * superposition of the complexes, together with alignment quality scores.
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
		List<SubunitCluster> c1 = SubunitClusterer.cluster(s1, cParams);
		List<SubunitCluster> c2 = SubunitClusterer.cluster(s2, cParams);

		// STEP 2: match each subunit cluster between groups O(N^2*L^2) - inter
		Map<Integer, Integer> clusterMap = new HashMap<Integer, Integer>();
		for (int i = 0; i < c1.size(); i++) {
			for (int j = 0; j < c2.size(); j++) {

				if (clusterMap.keySet().contains(i))
					break;
				if (clusterMap.values().contains(j))
					continue;

				switch (cParams.getClustererMethod()) {

				case IDENTITY:
					if (c1.get(i).mergeIdentical(c2.get(j)))
						clusterMap.put(i, j);
					break;

				case SEQUENCE:
					try {
						if (c1.get(i).mergeSequence(c2.get(j),
								cParams.getSequenceIdentityThreshold(),
								cParams.getCoverageThreshold()))
							clusterMap.put(i, j);
					} catch (CompoundNotFoundException e) {
						logger.warn("Could compare by Sequence. {}",
								e.getMessage());
					}
					break;

				default: // case STRUCTURE:
					if (c1.get(i).mergeStructure(c2.get(j),
							cParams.getRmsdThreshold(),
							cParams.getCoverageThreshold()))
						clusterMap.put(i, j);
					break;
				}
			}
		}

		// STEP 3: Align the assemblies for each cluster match O(L^2*N+N^2*L)
		for (int globalKey : clusterMap.keySet()) {

			// Obtain the clusters
			SubunitCluster clust1 = c1.get(globalKey);
			SubunitCluster clust2 = c2.get(clusterMap.get(globalKey));

			// Take the cluster match as reference and obtain transformation
			int index1 = 0;
			int index2 = clust1.size() - clust2.size();

			Atom[] atoms1 = clust1.getAlignedAtomsSubunit(index1);
			Atom[] atoms2 = clust1.getAlignedAtomsSubunit(index2);

			SVDSuperimposer svd = new SVDSuperimposer(atoms1, atoms2);
			Matrix4d trans = svd.getTransformation();

			// Map the subunits of each cluster to their spatial equivalents
			Map<Integer, Map<Integer, Integer>> clustSubunitMap = new HashMap<Integer, Map<Integer, Integer>>();

			for (int key : clusterMap.keySet()) {

				// Obtain the clusters
				clust1 = c1.get(key);
				clust2 = c2.get(clusterMap.get(key));

				// Take the cluster match as reference and obtain transformation
				index1 = 0;
				index2 = clust1.size() - clust2.size();

				// Map the subunits of the cluster to their spatial equivalents
				Map<Integer, Integer> subunitMap = new HashMap<Integer, Integer>();

				for (int i = 0; i < index2; i++) {
					for (int j = index2; j < clust1.size(); j++) {

						if (subunitMap.keySet().contains(i))
							break;
						if (subunitMap.values().contains(j))
							continue;

						// Obtain centroids and transform the second
						Atom centr1 = Calc.getCentroid(clust1
								.getAlignedAtomsSubunit(i));
						Atom centr2 = Calc.getCentroid(clust1
								.getAlignedAtomsSubunit(j));
						Calc.transform(centr2, trans);

						if (Calc.getDistance(centr1, centr2) < aParams
								.getdCutoff())
							subunitMap.put(i, j);
					}
				}

				clustSubunitMap.put(key, subunitMap);
			}

			// Unfold the nested map into subunit map and alignment
			Map<Integer, Integer> subunitMap = new HashMap<Integer, Integer>();
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
			} else if (subunitMap.size() == result.getSubunitMap().size()) {
				if (result.getAlignment() == null) {
					result.setSubunitMap(subunitMap);
					result.setAlignment(msa);
				} else if (msa.getScore(MultipleAlignmentScorer.RMSD) < result
						.getRmsd()) {
					result.setSubunitMap(subunitMap);
					result.setAlignment(msa);
				}
			}

		}

		return result;
	}
}
