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
package org.biojava.nbio.structure.cluster;

import java.util.ArrayList;
import java.util.List;

import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.symmetry.core.Stoichiometry;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * The SubunitClusterer takes as input a collection of {@link Subunit} and
 * returns a collection of {@link SubunitCluster}.
 * 
 * @author Aleix Lafita
 * @since 5.0.0
 * 
 */
public class SubunitClusterer {

	private static final Logger logger = LoggerFactory
			.getLogger(SubunitClusterer.class);

	/** Prevent instantiation **/
	private SubunitClusterer() {
	}

	public static Stoichiometry cluster(Structure structure,
			SubunitClustererParameters params) {
		List<Subunit> subunits = SubunitExtractor.extractSubunits(structure,
				params.getAbsoluteMinimumSequenceLength(),
				params.getMinimumSequenceLengthFraction(),
				params.getMinimumSequenceLength());
		return cluster(subunits, params);
	}

	public static Stoichiometry cluster(List<Subunit> subunits,
			SubunitClustererParameters params) {

		// The collection of clusters to return
		List<SubunitCluster> clusters = new ArrayList<SubunitCluster>();

		if (subunits.size() == 0)
			return new Stoichiometry(clusters);

		// First generate a new cluster for each Subunit
		for (Subunit s : subunits)
			clusters.add(new SubunitCluster(s));

		if (params.getClustererMethod() == SubunitClustererMethod.SEQUENCE ||
				params.getClustererMethod() == SubunitClustererMethod.SEQUENCE_STRUCTURE) {
			// Now merge clusters by SEQUENCE
			for (int c1 = 0; c1 < clusters.size(); c1++) {
				for (int c2 = clusters.size() - 1; c2 > c1; c2--) {
					try {
						if (clusters.get(c1).mergeSequence(clusters.get(c2), params)) {
							clusters.remove(c2);
						}

					} catch (CompoundNotFoundException e) {
						logger.warn("Could not merge by Sequence. {}",
								e.getMessage());
					}
				}
			}
		}

		if (params.getClustererMethod() == SubunitClustererMethod.STRUCTURE ||
				params.getClustererMethod() == SubunitClustererMethod.SEQUENCE_STRUCTURE) {
			// Now merge clusters by STRUCTURE
			for (int c1 = 0; c1 < clusters.size(); c1++) {
				for (int c2 = clusters.size() - 1; c2 > c1; c2--) {
					try {
						if (clusters.get(c1).mergeStructure(clusters.get(c2), params)) {
							clusters.remove(c2);
						}
					} catch (StructureException e) {
						logger.warn("Could not merge by Structure. {}", e.getMessage());
					}
				}
			}
		}

		if (params.isInternalSymmetry()) {
			// Now divide clusters by their INTERNAL SYMMETRY
			for (int c = 0; c < clusters.size(); c++) {
				try {
					clusters.get(c).divideInternally(params);
				} catch (StructureException e) {
					logger.warn("Error analyzing internal symmetry. {}",
							e.getMessage());
				}
			}

			// After internal symmetry merge again by structural similarity
			// Use case: C8 propeller with 3 chains with 3+3+2 repeats each
			for (int c1 = 0; c1 < clusters.size(); c1++) {
				for (int c2 = clusters.size() - 1; c2 > c1; c2--) {
					try {
						if (clusters.get(c1).mergeStructure(clusters.get(c2), params))
							clusters.remove(c2);
					} catch (StructureException e) {
						logger.warn("Could not merge by Structure. {}",
								e.getMessage());
					}
				}
			}
		}

		return new Stoichiometry(clusters);
	}
}
