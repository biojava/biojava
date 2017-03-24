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
package org.biojava.nbio.core.util;

import java.util.Set;
import java.util.stream.Collectors;

import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.RNASequence;
import org.biojava.nbio.core.sequence.compound.AmbiguityDNACompoundSet;
import org.biojava.nbio.core.sequence.compound.AmbiguityDNARNAHybridCompoundSet;
import org.biojava.nbio.core.sequence.compound.AmbiguityRNACompoundSet;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava.nbio.core.sequence.compound.DNACompoundSet;
import org.biojava.nbio.core.sequence.compound.RNACompoundSet;
import org.biojava.nbio.core.sequence.template.CompoundSet;
import org.biojava.nbio.core.sequence.template.Sequence;
import org.slf4j.LoggerFactory;

public class SequenceTools {
	private static final org.slf4j.Logger logger = LoggerFactory.getLogger(SequenceTools.class);

	protected static final String NUCLEOTIDE_LETTERS = "GCTAUXN";


	/**
	 * Cyclically permute the characters in {@code string} <em>forward</em> by {@code n} elements.
	 * @param string The string to permute
	 * @param n The number of characters to permute by; can be positive or negative; values greater than the length of the array are acceptable
	 */
	public static String permuteCyclic(String string, int n) {
		// single letters are char[]; full names are Character[]
		Character[] permuted = new Character[string.length()];
		char[] c = string.toCharArray();
		Character[] charArray = new Character[c.length];
		for (int i = 0; i < c.length; i++) {
			charArray[i] = c[i];
		}
		permuteCyclic(charArray, permuted, n);
		char[] p = new char[permuted.length];
		for (int i = 0; i < p.length; i++) {
			p[i] = permuted[i];
		}
		return String.valueOf(p);
	}

	/**
	 * Cyclically permute {@code array} <em>forward</em> by {@code n} elements.
	 * @param array The original result; will not be changed
	 * @param fill The permuted result will be filled into this array
	 * @param n The number of elements to permute by; can be positive or negative; values greater than the length of the array are acceptable
	 */
	public static <T> void permuteCyclic(T[] array, T[] fill, int n) {
		if (array.length != fill.length) throw new IllegalArgumentException("Lengths do not match");
		if (n < 0) n = array.length + n;
		while (n > array.length) {
			n -= array.length;
		}
		for (int i = 0; i < array.length; i++) {
			if (i + n < array.length) {
				fill[i] = array[i + n];
			} else {
				fill[i] = array[i - array.length + n];
			}
		}
	}

	/**
	 * Calculate the percentage of GCTAUXN in the sequence.
	 * Note that this omits some of the more esoteric ambiguity codes
	 * @param sequence
	 * @return
	 */
	public static int percentNucleotideSequence(String sequence)
	{
			if (sequence == null || sequence.length() == 0) return 0;

			int l = sequence.length();
			int n =0;

			for (int i = 0; i < l; i++)
			{
					if (NUCLEOTIDE_LETTERS.indexOf(sequence.charAt(i)) < 0 )
					{
							continue;
					}
					n++;
			}
			return (100 * n) / l;
	}

	/**
	 * Check if a sequence consists entirely of nucleotides (including all ambiguous codes)
	 * @param sequence
	 * @return
	 */
	public static boolean isNucleotideSequence(String sequence)
	{
			return isCompatibleSequence(sequence, AmbiguityDNARNAHybridCompoundSet.getDNARNAHybridCompoundSet());
	}
	/**
	 * Check if a sequence constists entirely of letters compatible with the specified compound set
	 * @param sequence
	 * @param set
	 * @return
	 */
	public static boolean isCompatibleSequence(String sequence, CompoundSet<?> set) {
		if (sequence == null || sequence.length() == 0) return false;

		Set<String> validLetters = set
				.getAllCompounds().stream()
				.map((nuc) -> nuc.getShortName())
				.collect(Collectors.toSet());
		return sequence.chars().allMatch((chr) -> validLetters.contains(chr) );
	}
	/**
	 * General method to create a sequence (of unknown type) from a string. Gaps are ignored.
	 * @param gappedSequenceString
	 * @return
	 */
	public static Sequence<?> getSequenceFromString(String gappedSequenceString){
		if (gappedSequenceString == null) return null;

		String sequenceString = gappedSequenceString.replace("-", "");

		try {
			if( isCompatibleSequence(sequenceString, DNACompoundSet.getDNACompoundSet()))
				return new DNASequence(sequenceString, DNACompoundSet.getDNACompoundSet());
			else if( isCompatibleSequence(sequenceString, RNACompoundSet.getRNACompoundSet()))
				return new RNASequence(sequenceString, RNACompoundSet.getRNACompoundSet());
			if( isCompatibleSequence(sequenceString, AmbiguityDNACompoundSet.getDNACompoundSet()))
				return new DNASequence(sequenceString, AmbiguityDNACompoundSet.getDNACompoundSet());
			else if( isCompatibleSequence(sequenceString, AmbiguityRNACompoundSet.getRNACompoundSet()))
				return new RNASequence(sequenceString, AmbiguityRNACompoundSet.getRNACompoundSet());
			else
				return new ProteinSequence(sequenceString, AminoAcidCompoundSet.getAminoAcidCompoundSet());
		} catch( CompoundNotFoundException e) {
			// Should never happen, since ProteinSequence accepts all 26 letters
			logger.error("Unrecognized compound");
			return null;
		}
	}
}
