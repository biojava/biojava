package org.biojava.bio.structure.domain.pdp;

public class PDPDistanceMatrix {
	int[][] dist;
	int nclose;
	int[] iclose ;
	int[] jclose ;
	
	public PDPDistanceMatrix(){
		
	}

	public int[][] getDist() {
		return dist;
	}

	public void setDist(int[][] dist) {
		this.dist = dist;
	}

	public int getNclose() {
		return nclose;
	}

	public void setNclose(int nclose) {
		this.nclose = nclose;
	}

	public int[] getIclose() {
		return iclose;
	}

	public void setIclose(int[] iclose) {
		this.iclose = iclose;
	}

	public int[] getJclose() {
		return jclose;
	}

	public void setJclose(int[] jclose) {
		this.jclose = jclose;
	}
	
	
		

}
