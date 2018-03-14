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
package org.biojava.nbio.aaproperties;

import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * This class provides the protein properties at the level of individual amino acids.
 *
 * @author Yana Valasatava
 */
public class AminoAcidProperties {

	private static final Set<String> negChargedAAs = Stream.of("D", "E").collect(Collectors.toSet());
	private static final Set<String> posChargedAAs = Stream.of("K", "R", "H").collect(Collectors.toSet());
	private static final Set<String> polarAAs = Stream.of("D", "E", "K", "R", "H", "N", "Q", "S", "T", "Y")
			.collect(Collectors.toSet());
	
	/**
	 * At pH=7, two are negative charged: aspartic acid (Asp, D) and glutamic acid (Glu, E) (acidic side chains), 
	 * and three are positive charged: lysine (Lys, K), arginine (Arg, R) and histidine (His, H) (basic side chains).
	 * 
	 * @param aa The one-letter amino acid code 
	 * @return true if amino acid is charged
	 */
	public static final boolean isCharged(char aa) {
		if (negChargedAAs.contains(String.valueOf(aa))) {
			return true;
		}
		else if (posChargedAAs.contains(String.valueOf(aa))) {
			return true;
		}
		return false;
	}
	
	/**
	 * Returns the charge of amino acid. At pH=7, two are negative charged: aspartic acid (Asp, D) and glutamic acid (Glu, E) (acidic side chains), 
	 * and three are positive charged: lysine (Lys, K), arginine (Arg, R) and histidine (His, H) (basic side chains).
	 * 
	 * @param aa The one-letter amino acid code 
	 * @return the charge of amino acid (1 if positively charged, -1 if negatively charged, 0 if not charged)
	 */
	public static final int getChargeOfAminoAcid(char aa) {
		if (negChargedAAs.contains(String.valueOf(aa))) {
			return -1;
		}
		else if (posChargedAAs.contains(String.valueOf(aa))) {
			return 1;
		}
		return 0;
	}
	
	/**
	 * There are 10 amino acids: D, E, H, K, R, N, Q, S, T, Y, that are polar.
	 * 
	 * @param aa The one-letter amino acid code 
	 * @return true if amino acid is polar
	 */
	public static final boolean isPolar(char aa) {
		if (polarAAs.contains(String.valueOf(aa))) {
			return true;
		}
		return false;
	}
	
	/**
	 * There are 10 amino acids: D, E, H, K, R, N, Q, S, T, Y, that are polar.
	 * 
	 * @param aa The one-letter amino acid code
	 * @return the polarity of amino acid (1 if polar, 0 if not polar)
	 */
	public static final int getPolarityOfAminoAcid(char aa) {
		if (polarAAs.contains(String.valueOf(aa))) {
			return 1;
		}
		return 0;
	}
}
