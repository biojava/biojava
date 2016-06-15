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

public class Convert2Polarizability extends Convertor{
	/**
	 * Class for the conversion of protein sequence into polarizability
	 */

	@Override
	public char convert(char c){
		switch(c){
		case 'G':case 'A':case 'S':case 'D':case 'T':
			return group1;//Polarizability value 0-0.08
		case 'C':case 'P':case 'N':case 'V':case 'E':case 'Q':case 'I':case 'L':
			return group2;//Polarizability value 0.128-0.186
		case 'K':case 'M':case 'H':case 'F':case 'R':case 'Y':case 'W':
			return group3;//Polarizability value 0.219-0.409
		default:
			return unknownGroup;//Non-standard AA
		}
	}

	private static String[] subCategory = {"Value_0-0.08", "Value_0.128-0.186", "Value_0.219-0.409"};
	@Override
	public String[] getGrouping(){return subCategory;}
	@Override
	public String getAttribute(){return "Polarizability";}
}
