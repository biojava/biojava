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

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.HashSet;
import java.util.Set;

/**
 * This is a utility class that contains utility methods which will facilitates the coding of other methods
 *
 * @author kohchuanhock
 * @version 2011.08.22
 * @since 3.0.2
 */
public class Utils {

	private final static Logger logger = LoggerFactory.getLogger(Utils.class);

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
		for(char c:sequence.toCharArray()){
			if(!cSet.contains(c)) return true;
		}
		return false;
	}

	/**
	 * Return the number of invalid characters in sequence.
	 *
	 * @param sequence
	 * 		protein sequence to count for invalid characters.
	 * @param cSet
	 * 		the set of characters that are deemed valid.
	 * @param ignoreCase
	 * 		indicates if cases should be ignored
	 * @return
	 * 		the number of invalid characters in sequence.
	 */
	public final static int getNumberOfInvalidChar(String sequence, Set<Character> cSet, boolean ignoreCase){
		int total = 0;
		char[] cArray;
		if(ignoreCase) cArray = sequence.toUpperCase().toCharArray();
		else cArray = sequence.toCharArray();
		if(cSet == null) cSet = PeptideProperties.standardAASet;
		for(char c:cArray){
			if(!cSet.contains(c)) total++;
		}
		return total;
	}

	/**
	 * Returns a new sequence with all invalid characters being replaced by '-'.
	 * Note that any character outside of the 20 standard protein amino acid codes are considered as invalid.
	 *
	 * @param sequence
	 * 		protein sequence to be clean
	 * @param cSet
	 * 		user defined characters that are valid. Can be null. If null, then 20 standard protein amino acid codes will be considered as valid.
	 * @return
	 * 		a new sequence with all invalid characters being replaced by '-'.
	 */
	public final static String cleanSequence(String sequence, Set<Character> cSet){
		Set<Character> invalidCharSet = new HashSet<Character>();
		StringBuilder cleanSeq = new StringBuilder();
		if(cSet == null) cSet = PeptideProperties.standardAASet;
		for(char c:sequence.toCharArray()){
			if(!cSet.contains(c)){
				cleanSeq.append("-");
				invalidCharSet.add(c);
			}else{
				cleanSeq.append(c);
			}
		}

		// TODO: Should be StringJoiner once JDK8 used
		StringBuilder stringBuilder = new StringBuilder();
		for(char c: invalidCharSet){
			stringBuilder.append("\'" + c + "\'");
		}
		stringBuilder.deleteCharAt(stringBuilder.length()-1);
		stringBuilder.append(" are being replaced with '-'");
		logger.warn(stringBuilder.toString());

		return cleanSeq.toString();
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
			containInvalid = sequence != null && doesSequenceContainInvalidChar(sequence, PeptideProperties.standardAASet);
		}
		if(containInvalid){
			String cSeq = cleanSequence(sequence, cSet);
			logger.warn("There exists invalid characters in the sequence. Computed results might not be precise.");
			logger.warn("To remove this warning: Please use org.biojava.nbio.aaproperties.Utils.cleanSequence to clean sequence.");

			return cSeq;
		}
		else{
			return sequence;
		}
	}
}
