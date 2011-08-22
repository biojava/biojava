package org.biojava3.aaproperties.profeat.convertor;

public class Convert2Hydrophobicity extends Convertor {
	/**
	 * Class for the conversion of protein sequence into hydrophobicity 
	 */
	
	@Override
	public char convert(char c){
		switch(c){
		case 'R':case 'K':case 'E':case 'D':case 'Q':case 'N':
			return group1;//Polar
		case 'G':case 'A':case 'S':case 'T':case 'H':case 'P':case 'Y':
			return group2;//Neutral
		case 'C':case 'L':case 'V':case 'I':case 'M':case 'F':case 'W':
			return group3;//Hydrophobicity
		default:
			return unknownGroup;//Non-standard AA
		}
	}
	private static String[] subCategory = {"Polar", "Neutral", "Hydrophobicity"};
	@Override
	public String[] getGrouping(){return subCategory;}
	@Override
	public String getAttribute(){return "Hydrophobicity";}
}
