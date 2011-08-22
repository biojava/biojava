package org.biojava3.aaproperties.profeat.convertor;

public class Convert2NormalizedVanDerWaalsVolume extends Convertor{
	/**
	 * Class for the conversion of protein sequence into normalized van der waals volume
	 */
	
	@Override
	public char convert(char c){
		switch(c){
		case 'G':case 'A':case 'S':case 'T':case 'P':case 'D': case 'C':
			return group1;//Volume Range 0-2.78
		case 'N':case 'V':case 'E':case 'Q':case 'I':case 'L':
			return group2;//Volume Range 2.95-4.0
		case 'M':case 'H':case 'K':case 'F':case 'R':case 'Y':case 'W':
			return group3;//Volume Range 4.03-8.08
		default:
			return unknownGroup;//Non-standard AA
		}
	}
	private static String[] subCategory = {"Range_0-2.78", "Range_2.95-4.0", "Range_4.03-8.08"};
	@Override
	public String[] getGrouping(){return subCategory;}
	@Override
	public String getAttribute(){return "Normalized van der waals Volume";}
}
