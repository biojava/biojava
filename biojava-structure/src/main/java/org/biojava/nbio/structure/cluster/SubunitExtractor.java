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

import org.biojava.nbio.structure.*;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * The SubunitExtractor extracts the information of each protein {@link Chain}
 * in a {@link Structure} and converts them into a List of {@link Subunit}.
 * 
 * @author Peter Rose
 * @author Aleix Lafita
 * @since 5.0.0
 * 
 */
public class SubunitExtractor {

	private static final Logger logger = LoggerFactory
			.getLogger(SubunitExtractor.class);

	/** Prevent instantiation **/
	private SubunitExtractor() {
	}

	/**
	 * Extract the information of each protein Chain in a Structure and converts
	 * them into a List of Subunit. The name of the Subunits is set to
	 * {@link Chain#getId()}.
	 * 
	 * 
	 * @param structure
	 *            Structure object with protein Chains
	 * @param absMinLen
	 *            {@link SubunitClustererParameters#getAbsoluteMinimumSequenceLength()}
	 * @param fraction
	 *            {@link SubunitClustererParameters#getMinimumSequenceLengthFraction()}
	 * @param minLen
	 *            {@link SubunitClustererParameters#getMinimumSequenceLength()}
	 * @return List of Subunits
	 */
	public static List<Subunit> extractSubunits(Structure structure,
			int absMinLen, double fraction, int minLen) {

		// The extracted subunit container
		List<Subunit> subunits = new ArrayList<Subunit>();

		for (Chain c : structure.getPolyChains()) {
			// Only take protein chains
			if (c.isProtein()) {
				Atom[] ca = StructureTools.getRepresentativeAtomArray(c);
				logger.debug("Chain " + c.getId() + "; CA Atoms: " + ca.length + "; SEQRES: " + c.getSeqResSequence());
				if (ca.length==0)
					continue;
				subunits.add(new Subunit(ca, c.getId(), null, structure));
			}
		}

		// Calculate the minimum length of a Subunit
		int adjustedMinLen = calcAdjustedMinimumSequenceLength(subunits,
				absMinLen, fraction, minLen);
		logger.debug("Adjusted minimum sequence length: " + adjustedMinLen);

		// Filter out short Subunits
		for (int s = subunits.size() - 1; s >= 0; s--) {
			if (subunits.get(s).size() < adjustedMinLen)
				subunits.remove(s);
		}

		return subunits;
	}

	/**
	 * Returns an adapted minimum sequence length. This method ensure that
	 * structure that only have short chains are not excluded by the
	 * minimumSequenceLength cutoff value.
	 * 
	 * @return adjustedMinimumSequenceLength
	 */
	private static int calcAdjustedMinimumSequenceLength(
			List<Subunit> subunits, int absMinLen, double fraction, int minLen) {

		int maxLength = Integer.MIN_VALUE;
		int minLength = Integer.MAX_VALUE;

		// Extract the length List, the min and the max
		List<Integer> lengths = new ArrayList<Integer>();
		for (int i = 0; i < subunits.size(); i++) {
			if (subunits.get(i).size() >= absMinLen) {
				maxLength = Math.max(subunits.get(i).size(), maxLength);
				minLength = Math.min(subunits.get(i).size(), minLength);
				lengths.add(subunits.get(i).size());

			}
		}

		int adjustedMinimumSequenceLength = minLen;

		if (lengths.size() < 2)
			return adjustedMinimumSequenceLength;

		// Calculate the median of the lengths
		double median = 0;
		Collections.sort(lengths);
		if (lengths.size() % 2 == 1) {
			int middle = (lengths.size() - 1) / 2;
			median = lengths.get(middle);
		} else {
			int middle2 = lengths.size() / 2;
			int middle1 = middle2 - 1;
			median = 0.5 * (lengths.get(middle1) + lengths.get(middle2));
		}

		// If the median * fraction is lower than the minLength
		if (minLength >= median * fraction) {
			adjustedMinimumSequenceLength = Math.min(minLength, minLen);
		}

		return adjustedMinimumSequenceLength;
	}
}
