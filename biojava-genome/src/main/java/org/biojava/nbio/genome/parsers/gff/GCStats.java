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
package org.biojava.nbio.genome.parsers.gff;

import org.biojava.nbio.core.sequence.DNASequence;

import java.util.Collection;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class GCStats {

	public static double getGCStats(Collection<DNASequence> sequences) {
		double gcCount = 0;
		double total = 0;

		for (DNASequence sequence : sequences) {
			char[] dna = sequence.toString().toCharArray();
			for (char d : dna) {
				if (d == 'G' || d == 'C' || d == 'g' || d == 'c') {
					gcCount++;
				}
				total++;
			}
		}

		return (gcCount / total) * 100.0;
	}

	public static double getGCStatsString(Collection<String> sequences) {
		double gcCount = 0;
		double total = 0;

		for (String sequence : sequences) {
			char[] dna = sequence.toCharArray();
			for (char d : dna) {
				if (d == 'G' || d == 'C' || d == 'g' || d == 'c') {
					gcCount++;
				}
				total++;
			}
		}

		return (gcCount / total) * 100.0;
	}
}
