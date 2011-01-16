package org.biojava.bio.structure.domain.pdp;

import java.util.Arrays;


public class CutSites {
	int ncuts;
	int[] cut_sites;
	
	public CutSites(){
		ncuts = 0;
		
		cut_sites = new int[PDPParameters.MAX_CUTS];
	}
	
	@Override
	public String toString() {
		return "CutSites [ncuts=" + ncuts + ", cut_sites="
				+ Arrays.toString(cut_sites) + "]";
	}
	public int getNcuts() {
		return ncuts;
	}
	public void setNcuts(int ncuts) {
		this.ncuts = ncuts;
	}
	public int[] getCut_sites() {
		return cut_sites;
	}
	public void setCut_sites(int[] cut_sites) {
		this.cut_sites = cut_sites;
	}
	
	
	
}
