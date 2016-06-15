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

public class Convert2Polarity extends Convertor{
	/**
	 * Class for the conversion of protein sequence into polarity
	 */

	@Override
	public char convert(char c){
		switch(c){
		case 'L':case 'I':case 'F':case 'W':case 'C':case 'M':case 'V':case 'Y':
			return group1;//Polarity value 4.9-6.2
		case 'P':case 'A':case 'T':case 'G':case 'S':
			return group2;//Polarity value 8.0-9.2
		case 'H':case 'Q':case 'R':case 'K':case 'N':case 'E':case 'D':
			return group3;//Polarity value 10.4-13.0
		default:
			return unknownGroup;//Non-standard AA
		}
	}
	private static String[] subCategory = {"Value_4.9-6.2", "Value_8.0-9.2", "Value_10.4-13.0"};
	@Override
	public String[] getGrouping(){return subCategory;}
	@Override
	public String getAttribute(){return "Polarity";}
}
