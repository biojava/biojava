package org.biojava.bio.structure.domain.pdp;



public class CutValues {
	public double s_min;
	public int site2;
	public boolean first_cut;

	public double AD;

	public CutValues(){
		s_min = 100;
		site2 = 0;
		first_cut = true;
	}
	
	@Override
	public String toString() {
		return "CutValues [s_min=" + s_min + ", site2=" + site2 +
		", AD=" + AD
		+ "]";
	}



}
