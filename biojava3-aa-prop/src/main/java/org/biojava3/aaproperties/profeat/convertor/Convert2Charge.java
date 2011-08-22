package org.biojava3.aaproperties.profeat.convertor;

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
