package org.biojava3.aaproperties.profeat.convertor;

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
