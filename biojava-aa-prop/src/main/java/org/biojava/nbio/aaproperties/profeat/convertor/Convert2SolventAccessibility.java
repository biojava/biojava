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

public class Convert2SolventAccessibility extends Convertor{
	/**
	 * Class for the conversion of protein sequence into solvent accessibility
	 */
	@Override
	public char convert(char c){
		switch(c){
		case 'A':case 'L':case 'F':case 'C':case 'G':case 'I':case 'V':case 'W':
			return group1;//Buried
		case 'R':case 'K':case 'Q':case 'E':case 'N':case 'D':
			return group2;//Exposed
		case 'M':case 'P':case 'S':case 'T':case 'H':case 'Y':
			return group3;//Intermediate
		default:
			return unknownGroup;//Non-standard AA
		}
	}

	private static String[] subCategory = {"Buried", "Exposed", "Intermediate"};
	@Override
	public String[] getGrouping(){return subCategory;}
	@Override
	public String getAttribute(){return "Solvent Accessibility";}
}
