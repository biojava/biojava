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

public class Convert2SecondaryStructure extends Convertor{
	/**
	 * Class for the conversion of protein sequence into secondary structure
	 */
	@Override
	public char convert(char c){
		switch(c){
		case 'E':case 'A':case 'L':case 'M':case 'Q':case 'K':case 'R':case 'H':
			return group1;//Helix
		case 'V':case 'I':case 'Y':case 'C':case 'W':case 'F':case 'T':
			return group2;//Strand
		case 'G':case 'N':case 'P':case 'S':case 'D':
			return group3;//Coil
		default:
			return unknownGroup;//Non-standard AA
		}
	}
	private static String[] subCategory = {"Helix", "Strand", "Coil"};
	@Override
	public String[] getGrouping(){return subCategory;}
	@Override
	public String getAttribute(){return "Secondary Structure";}
}
