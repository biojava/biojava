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

import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.template.Sequence;

public class SequenceTools {

	protected static final String NUCLEOTIDE_LETTERS = "GCTAUXN";

	/**
	 * Cyclically permute the characters in {@code string} <em>forward</em> by {@code n} elements.
	 * @param string The string to permute
	 * @param n The number of characters to permute by; can be positive or negative; values greater than the length of the array are acceptable
	 */
	public static String permuteCyclic(String string, int n) {
		String toMutate = string + string;
		n = n % string.length();
		if (n < 0) {
			n = string.length() + n;
		}
		return toMutate.substring(n, n + string.length());
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

	public static int percentNucleotideSequence(String sequence)
	{
			if (sequence == null || sequence.length() == 0) return 0;

			int l = sequence.length();
			int n =0;

			for (int i = 0; i < l; i++)
			{
					if (NUCLEOTIDE_LETTERS.indexOf(sequence.charAt(i)) < 0)
					{
							continue;
					}
					n++;
			}
			return (100 * n) / l;
	}

	public static boolean isNucleotideSequence(String sequence)
	{
			if (sequence == null || sequence.length() == 0) return false;

			int l = sequence.length();
			for (int i = 0; i < l; i++)
			{
					if (NUCLEOTIDE_LETTERS.indexOf(sequence.charAt(i)) < 0)
					{
							return false;
					}
			}
			return true;
	}

	/**
	 * Attempts to parse String as a DNA sequence first.<br/>
	 * If this fails it tries to  parse as a ProteinSequence.
	 * <br/>
	 * This method does not attempt to create an RNASequence.
	 * <p>
	 * Also, a sequence such as 'ATCGTA' which is both a
	 * peptide sequence and a DNA sequence, will always be returned 
	 * as a DNA sequence.
	 * </p>
	 * <p>
	 * An empty string argument returns a ProteinSequence of length 0.
	 * A null argument throws a {@link NullPointerException}
	 * @param sequence
	 * @return Either a DNASequence or a ProteinSequence
	 * @throws CompoundNotFoundException
	 */
	public Sequence<?> getSequenceFromString(String sequence) throws CompoundNotFoundException {


		if( isNucleotideSequence(sequence)) {
			return  new DNASequence(sequence);
		} else {
			return new ProteinSequence(sequence);

		}

	}

	/** A method to check whether an array of sequences contains at least two sequences having an equal length.
	 *
	 * @param sequences the array of {@link org.biojava.nbio.core.sequence.ProteinSequence} sequences
	 * @return true if any two sequences are of an equal length
	 */
	public static boolean equalLengthSequences(ProteinSequence[] sequences) {

		for (int i=0; i<sequences.length-1; i++) {
			if (sequences[i]==null)
				continue;
			for (int j=i+1; j<sequences.length; j++) {
				if (sequences[j]==null)
					continue;
				if (sequences[i].getLength() == sequences[j].getLength())
					return true;
			}
		}
		return false;
	}
}
