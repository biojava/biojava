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
package org.biojava.nbio.phylo;

import org.forester.evoinference.distance.NeighborJoining;
import org.forester.evoinference.matrix.distance.BasicSymmetricalDistanceMatrix;
import org.forester.phylogeny.Phylogeny;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * The TreeConstructor uses the forester library to build different types of
 * phylogenetic trees.
 *
 * @author Scooter Willis
 * @author Aleix Lafita
 *
 */
public class TreeConstructor {

	private static final Logger logger = LoggerFactory
			.getLogger(TreeConstructor.class);

	/** Prevent instantiation */
	private TreeConstructor() {}

	public static Phylogeny distanceTree(BasicSymmetricalDistanceMatrix distM,
			TreeConstructorType constructor) {

		Phylogeny p = null;
		switch (constructor) {
		case NJ:
			NeighborJoining nj = NeighborJoining.createInstance();
			p = nj.execute(distM);
			p.setType(TreeType.DISTANCE.name);
			break;
		default:
			logger.warn("Only NJ Tree Constructor Supported!");
			break;
		}
		logger.info("Tree Completed");
		return p;
	}
}
