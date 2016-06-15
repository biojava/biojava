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
package org.biojava.nbio.aaproperties.profeat.convertor;

public class Convert2Charge extends Convertor{
	/**
	 * Class for the conversion of protein sequence into charge
	 */

	@Override
	public char convert(char c){
		switch(c){
		case 'K':case 'R':
			return group1;//Positive
		case 'A':case 'N':case 'C':case 'Q':case 'G':case 'H':case 'I':case 'L':case 'M':case 'F':case 'P':case 'S':case 'T':case 'W':case 'Y':case 'V':
			return group2;//Neutral
		case 'D':case 'E':
			return group3;//Negative
		default:
			return unknownGroup;//Non-standard AA
		}
	}

	private static String[] subCategory = {"Positive", "Neutral", "Negative"};
	@Override
	public String[] getGrouping(){return Convert2Charge.subCategory;}
	@Override
	public String getAttribute(){return "Charge";}
}
