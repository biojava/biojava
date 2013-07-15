package org.biojava3.core.util;

import org.biojava3.core.sequence.DNASequence;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.template.Sequence;

public class SequenceTools {

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
    
    public Sequence<?> getSeqeunceFromString(String sequence){
    	
  
    	if( isNucleotideSequence(sequence)) {
    		return  new DNASequence(sequence);
    	} else {
    		return new ProteinSequence(sequence);
    	}
    	
    }
	
}
