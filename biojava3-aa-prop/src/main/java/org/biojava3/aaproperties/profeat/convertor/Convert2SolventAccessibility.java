package org.biojava3.aaproperties.profeat.convertor;

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
