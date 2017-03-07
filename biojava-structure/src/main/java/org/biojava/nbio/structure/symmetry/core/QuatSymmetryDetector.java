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
 * Created on 2013-05-23
 *
 */
package org.biojava.nbio.structure.symmetry.core;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.cluster.Subunit;
import org.biojava.nbio.structure.cluster.SubunitCluster;
import org.biojava.nbio.structure.cluster.SubunitClusterer;
import org.biojava.nbio.structure.cluster.SubunitClustererParameters;
import org.biojava.nbio.structure.symmetry.utils.PowerSet;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.*;

/**
 * Detects the symmetry (global, pseudo, internal and local) of protein
 * structures.
 * <p>
 * The {@link SubunitClustererParameters} determine the subunit definition and
 * clustering, while the {@link QuatSymmetryParameters} determine the calculated
 * symmetry results (point group and axes).
 *
 * @author Peter Rose
 * @author Aleix Lafita
 *
 */
public class QuatSymmetryDetector {

	private static final Logger logger = LoggerFactory
			.getLogger(QuatSymmetryDetector.class);

	/** Prevent instantiation **/
	private QuatSymmetryDetector() {
	}

	/**
	 * Calculate GLOBAL symmetry results. This means that all {@link Subunit}
	 * are included in the symmetry.
	 *
	 * @param structure
	 *            protein chains will be extracted as {@link Subunit}
	 * @param symmParams
	 *            quaternary symmetry parameters
	 * @param clusterParams
	 *            subunit clustering parameters
	 * @return GLOBAL quaternary structure symmetry results
	 */
	public static QuatSymmetryResults calcGlobalSymmetry(Structure structure,
			QuatSymmetryParameters symmParams,
			SubunitClustererParameters clusterParams) {
		List<SubunitCluster> clusters = SubunitClusterer.cluster(structure,
				clusterParams);
		return calcGlobalSymmetry(clusters, symmParams);
	}

	/**
	 * Calculate GLOBAL symmetry results. This means that all {@link Subunit}
	 * are included in the symmetry.
	 *
	 * @param subunits
	 *            list of {@link Subunit}
	 * @param symmParams
	 *            quaternary symmetry parameters
	 * @param clusterParams
	 *            subunit clustering parameters
	 * @return GLOBAL quaternary structure symmetry results
	 */
	public static QuatSymmetryResults calcGlobalSymmetry(
			List<Subunit> subunits, QuatSymmetryParameters symmParams,
			SubunitClustererParameters clusterParams) {
		List<SubunitCluster> clusters = SubunitClusterer.cluster(subunits,
				clusterParams);
		return calcGlobalSymmetry(clusters, symmParams);
	}

	/**
	 * Calculate GLOBAL symmetry results. This means that all {@link Subunit}
	 * are included in the symmetry.
	 *
	 * @param clusters
	 *            list of {@link SubunitCluster}
	 * @param symmParams
	 *            quaternary symmetry parameters
	 * @return GLOBAL quaternary structure symmetry results
	 */
	public static QuatSymmetryResults calcGlobalSymmetry(
			List<SubunitCluster> clusters, QuatSymmetryParameters symmParams) {
		return calcQuatSymmetry(clusters, symmParams);
	}

	/**
	 * Returns a List of LOCAL symmetry results. This means that a subset of the
	 * {@link SubunitCluster} is left out of the symmetry calculation. Each
	 * element of the List is one possible LOCAL symmetry result.
	 * <p>
	 * Determine local symmetry if global structure is: (1) asymmetric, C1; (2)
	 * heteromeric (belongs to more than 1 subunit cluster); (3) more than 2
	 * subunits (heteromers with just 2 chains cannot have local symmetry)
	 *
	 * @param structure
	 *            protein chains will be extracted as {@link Subunit}
	 * @param symmParams
	 *            quaternary symmetry parameters
	 * @param clusterParams
	 *            subunit clustering parameters
	 * @return List of LOCAL quaternary structure symmetry results. Empty if
	 *         none.
	 */
	public static List<QuatSymmetryResults> calcLocalSymmetries(
			Structure structure, QuatSymmetryParameters symmParams,
			SubunitClustererParameters clusterParams) {
		List<SubunitCluster> clusters = SubunitClusterer.cluster(structure,
				clusterParams);
		return calcLocalSymmetries(clusters, symmParams);
	}

	/**
	 * Returns a List of LOCAL symmetry results. This means that a subset of the
	 * {@link SubunitCluster} is left out of the symmetry calculation. Each
	 * element of the List is one possible LOCAL symmetry result.
	 * <p>
	 * Determine local symmetry if global structure is: (1) asymmetric, C1; (2)
	 * heteromeric (belongs to more than 1 subunit cluster); (3) more than 2
	 * subunits (heteromers with just 2 chains cannot have local symmetry)
	 *
	 * @param subunits
	 *            list of {@link Subunit}
	 * @param symmParams
	 *            quaternary symmetry parameters
	 * @param clusterParams
	 *            subunit clustering parameters
	 * @return List of LOCAL quaternary structure symmetry results. Empty if
	 *         none.
	 */
	public static List<QuatSymmetryResults> calcLocalSymmetries(
			List<Subunit> subunits, QuatSymmetryParameters symmParams,
			SubunitClustererParameters clusterParams) {
		List<SubunitCluster> clusters = SubunitClusterer.cluster(subunits,
				clusterParams);
		return calcLocalSymmetries(clusters, symmParams);
	}

	/**
	 * Returns a List of LOCAL symmetry results. This means that a subset of the
	 * {@link SubunitCluster} is left out of the symmetry calculation. Each
	 * element of the List is one possible LOCAL symmetry result.
	 * <p>
	 * Determine local symmetry if global structure is: (1) asymmetric, C1; (2)
	 * heteromeric (belongs to more than 1 subunit cluster); (3) more than 2
	 * subunits (heteromers with just 2 chains cannot have local symmetry)
	 *
	 * @param clusters
	 *            list of {@link SubunitCluster}
	 * @param symmParams
	 *            quaternary symmetry parameters
	 * @return List of LOCAL quaternary structure symmetry results. Empty if
	 *         none.
	 */
	public static List<QuatSymmetryResults> calcLocalSymmetries(
			List<SubunitCluster> clusters, QuatSymmetryParameters symmParams) {

		List<QuatSymmetryResults> localSymmetries = new ArrayList<QuatSymmetryResults>();

		// If it is homomeric return empty
		if (clusters.size() < 2)
			return localSymmetries;

		// If there are less than 3 or more than maximum Subunits return empty
		QuatSymmetrySubunits subunits = new QuatSymmetrySubunits(clusters);
		if (subunits.getSubunitCount() < 3
				|| subunits.getSubunitCount() > symmParams
						.getMaximumLocalSubunits())
			return localSymmetries;

		// Start local symmetry calculations
		long start = System.nanoTime();

		// Calculate the power Set of the clusters
		Set<Set<SubunitCluster>> powerSet = new PowerSet<SubunitCluster>()
				.powerSet(new HashSet<SubunitCluster>(clusters));

		int combinations = 1;
		for (Set<SubunitCluster> cluster : powerSet) {

			// Break if time limit passed
			double time = (System.nanoTime() - start) / 1000000000;
			if (time > symmParams.getLocalTimeLimit()) {
				logger.warn("Exceeded time limit for local symmetry "
						+ "calculations. {} seconds elapsed. "
						+ "Local symmetry results may be incomplete.", time);
				break;
			}

			// Break if maximum number of results exceeded
			if (localSymmetries.size() > symmParams.getMaximumLocalResults()) {
				logger.warn("Exceeded maximum number of local symmetry "
						+ "results. {} results calculated. "
						+ "Local symmetry results may be incomplete.",
						localSymmetries.size());
				break;
			}

			// Break if maximum number of tried combinations exceeded
			if (combinations > symmParams.getMaximumLocalCombinations()) {
				logger.warn("Exceeded maximum number of local combinations."
						+ " {} combinations tried. "
						+ "Local symmetry results may be incomplete.",
						combinations);
				break;
			}

			// Do not use empty set or identity set
			if (cluster.size() == 0 || cluster.size() == clusters.size())
				continue;

			List<SubunitCluster> localClusters = new ArrayList<SubunitCluster>(
					cluster);
			QuatSymmetryResults localResult = calcGlobalSymmetry(localClusters,
					symmParams);

			if (!localResult.getSymmetry().equals("C1")) {
				localResult.setLocal(true);
				localSymmetries.add(localResult);
			}

			combinations++;
		}

		return localSymmetries;
	}

	private static QuatSymmetryResults calcQuatSymmetry(
			List<SubunitCluster> clusters, QuatSymmetryParameters parameters) {

		QuatSymmetrySubunits subunits = new QuatSymmetrySubunits(clusters);

		if (subunits.getSubunitCount() == 0)
			return null;

		RotationGroup rotationGroup = null;
		SymmetryPerceptionMethod method = null;
		if (subunits.getFolds().size() == 1) {
			// no symmetry possible, create empty ("C1") rotation group
			method = SymmetryPerceptionMethod.NO_ROTATION;
			rotationGroup = new RotationGroup();
			rotationGroup.setC1(subunits.getSubunitCount());
		} else if (subunits.getSubunitCount() == 2
				&& subunits.getFolds().contains(2)) {
			method = SymmetryPerceptionMethod.C2_ROTATION;
			QuatSymmetrySolver solver = new C2RotationSolver(subunits,
					parameters);
			rotationGroup = solver.getSymmetryOperations();
		} else {
			method = SymmetryPerceptionMethod.ROTATION;
			QuatSymmetrySolver solver = new RotationSolver(subunits, parameters);
			rotationGroup = solver.getSymmetryOperations();
		}

		QuatSymmetryResults results = new QuatSymmetryResults(clusters,
				rotationGroup, method);

		String symmetry = results.getSymmetry();

		// asymmetric structures cannot be pseudosymmetric
		if (symmetry.equals("C1"))
			results.setPseudosymmetric(false);

		// Check structures with Cn symmetry (n = 1, ...) for helical symmetry
		if (symmetry.startsWith("C")) {
			HelixSolver hc = new HelixSolver(subunits,
					rotationGroup.getOrder(), parameters);
			HelixLayers helixLayers = hc.getSymmetryOperations();

			if (helixLayers.size() > 0) {
				// if the complex has no symmetry (c1) or has Cn (n>=2) symmetry
				// and the
				// helical symmetry has a lower RMSD than the cyclic symmetry,
				// set helical symmetry
				// If the RMSD for helical and cyclic symmetry is similar, a
				// slight preference is
				// given to the helical symmetry by the helixRmsdThreshold
				// parameter.
				double cRmsd = rotationGroup.getScores().getRmsd();
				double hRmsd = helixLayers.getScores().getRmsd();
				// System.out.println("cRMSD: " + cRmsd + " hRMSD: " + hRmsd);
				double deltaRmsd = hRmsd - cRmsd;
				if (symmetry.equals("C1")
						|| (!symmetry.equals("C1") && deltaRmsd <= parameters
								.getHelixRmsdThreshold())) {
					method = SymmetryPerceptionMethod.ROTO_TRANSLATION;
					results = new QuatSymmetryResults(clusters, helixLayers,
							method);
				}
			}
		}

		return results;
	}
}
