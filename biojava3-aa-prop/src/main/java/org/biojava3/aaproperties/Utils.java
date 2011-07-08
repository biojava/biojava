package org.biojava3.aaproperties;

import java.util.Set;

public class Utils {

	/**
	 * Returns a value with the desired number of decimal places.
	 * 
	 * @param d
	 * 		value to round
	 * @param c
	 * 		number of decimal places desired. 
	 * 		Must be greater or equal to zero, otherwise, the given value d would be returned without any modification.
	 * @return
	 * 		a value with the given number of decimal places.
	 */
	public final static double roundToDecimals(double d, int c) {
		if(c < 0) return d;
		double p = Math.pow(10,c);
		d = d * p;
		double tmp = Math.round(d);
		return tmp/p;
	}
	
	/**
	 * Checks if given sequence contains invalid characters. Returns true if invalid characters are found, else return false.
	 * Note that any characters are deemed as valid only if it is found in cSet.
	 * 
	 * @param sequence
	 * 		protein sequence to be check.
	 * @param cSet
	 * 		the set of characters that are deemed valid.
	 * @return
	 * 		true if invalid characters are found, else return false.
	 */
	public final static boolean doesSequenceContainInvalidChar(String sequence, Set<Character> cSet){
		String seq = sequence.toUpperCase();
		for(char c:seq.toCharArray()){
			if(cSet.contains(c) == false) return true;
		}
		return false;
	}
	
	/**
	 * Checks if given sequence contains invalid characters. Returns true if invalid characters are found, else return false.
	 * Note that any character outside of the 20 standard protein amino acid codes are considered as invalid.
	 * 
	 * @param sequence
	 * 		protein sequence to be check.
	 * @return
	 * 		true if invalid characters are found, else return false.
	 */
	public final static boolean doesSequenceContainInvalidChar(String sequence){
		String seq = sequence.toUpperCase();
		for(char c:seq.toCharArray()){
			if(PeptideProperties.standardAASet.contains(c) == false) return true;
		}
		return false;
	}
	
	/**
	 * Return the number of invalid characters in sequence.
	 * 
	 * @param sequence
	 * 		protein sequence to count for invalid characters.
	 * @return
	 * 		the number of invalid characters in sequence.
	 */
	public final static int getNumberOfInvalidChar(String sequence){
		String seq = sequence.toUpperCase();
		int total = 0;
		for(char c:seq.toCharArray()){
			if(PeptideProperties.standardAASet.contains(c) == false) total++;
		}
		return total;
	}

	/**
	 * Returns a new sequence with all invalid characters being replaced by '-'.
	 * Note that any character outside of the 20 standard protein amino acid codes are considered as invalid.
	 * 
	 * @param sequence
	 * 		protein sequence to be clean
	 * @return
	 * 		a new sequence with all invalid characters being replaced by '-'.
	 */
	public final static String cleanSequence(String sequence){
		String seq = sequence.toUpperCase();
		String cleanSeq = "";
		for(char c:seq.toCharArray()){
			if(PeptideProperties.standardAASet.contains(c) == false){
				cleanSeq += "-";
			}else{
				cleanSeq += c;
			}
		}
		return cleanSeq;
	}
	
	/**
	 * Checks if the sequence contains invalid characters. 
	 * Note that any character outside of the 20 standard protein amino acid codes are considered as invalid.
	 * If yes, it will return a new sequence where invalid characters are replaced with '-'.
	 * If no, it will simply return the input sequence. 
	 * 
	 * @param sequence
	 * 		protein sequence to be check for invalid characters.
	 * @return
	 * 		a sequence with no invalid characters.
	 */
	public static final String checkSequence(String sequence){
		return checkSequence(sequence, null);
	}
	
	/**
	 * Checks if the sequence contains invalid characters. 
	 * Note that any character outside of the 20 standard protein amino acid codes are considered as invalid.
	 * If yes, it will return a new sequence where invalid characters are replaced with '-'.
	 * If no, it will simply return the input sequence. 
	 * 
	 * @param sequence
	 * 		protein sequence to be check for invalid characters.
	 * @param cSet
	 * 		character set which define the valid characters.
	 * @return
	 * 		a sequence with no invalid characters.
	 */
	public static final String checkSequence(String sequence, Set<Character> cSet){
		boolean containInvalid = false;
		if(cSet != null){
			containInvalid = sequence != null && doesSequenceContainInvalidChar(sequence, cSet);
		}else{
			containInvalid = sequence != null && doesSequenceContainInvalidChar(sequence);
		}
		if(containInvalid){
			String cSeq = cleanSequence(sequence);
			System.err.println("Warning: There exists invalid characters in the sequence. Computed results might not be precise.");
			System.err.println("To remove this warning: Please use org.biojava3.aaproperties.Utils.cleanSequence to clean sequence.");
			return cSeq;
		}
		else{
			return sequence;
		}
	}
	
	public static void main(String[] args){
		String seq = "MTADGPCRELLCQLRAAVRHRWWCx";
		System.out.println(doesSequenceContainInvalidChar(seq));
		System.out.println(cleanSequence(seq));
	}
}
