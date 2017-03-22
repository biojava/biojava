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
	 * 
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
	 * 
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
	 * Returns the array of charges of each amino acid in a protein. At pH=7, two are negative charged: aspartic acid (Asp, D) and glutamic acid (Glu, E) (acidic side chains), 
	 * and three are positive charged: lysine (Lys, K), arginine (Arg, R) and histidine (His, H) (basic side chains).
	 * 
	 * @param aa The one-letter amino acid code
	 * 
	 * @return the array of charges of amino acids in the protein (1 if amino acid is positively charged, 
	 * -1 if negatively charged, 0 if not charged)
	 */
	public static final int[] getChargesOfAminoAcidsInProtein(String protein) {
		
		int[] charges = new int[protein.length()];
		for ( int i=0; i < protein.length(); i++ ) {
			char aa = protein.toCharArray()[i];
			charges[i] = getChargeOfAminoAcid(aa);
		}
		return charges;
	}
	
	/**
	 * There are 10 polar amino acids: D, E, H, K, R, N, Q, S, T, Y, that are polar.
	 * 
	 * @param aa The one-letter amino acid code
	 * 
	 * @return true if amino acid is polar
	 */
	public static final boolean isPolar(char aa) {
		if (polarAAs.contains(String.valueOf(aa))) {
			return true;
		}
		return false;
	}
	
	/**
	 * There are 10 polar amino acids: D, E, H, K, R, N, Q, S, T, Y, that are polar.
	 * 
	 * @param aa The one-letter amino acid code
	 * 
	 * @return the polarity of amino acid (1 if polar, 0 if not polar)
	 */
	public static final int getPolarityOfAminoAcid(char aa) {
		if (polarAAs.contains(String.valueOf(aa))) {
			return 1;
		}
		return 0;
	}
	
	/**
	 * Returns the array of polarity of each amino acid in a protein.
	 * 
	 * @param aa The one-letter amino acid code
	 * 
	 * @return the array of polarity of amino acids in the protein (1 if amino acid is polar, 0 if not)
	 */
	public static final int[] getPolarityOfAminoAcidsInProtein(String protein) {
		
		int[] polarity = new int[protein.length()];
		for ( int i=0; i < protein.length(); i++ ) {
			char aa = protein.toCharArray()[i];
			polarity[i] = getPolarityOfAminoAcid(aa);
		}
		return polarity;
	}
	
	public static void main(String[] args) {
		System.out.println(getChargeOfAminoAcid('D'));
		System.out.println(getChargeOfAminoAcid('K'));
		System.out.println(getChargeOfAminoAcid('A'));
	}
}
