package org.biojava.nbio.structure.secstruc;

/**
 * A secondary structure element (SSE) is an object representing a block
 * of sequential residues that share the same secondary structure type.
 * 
 * @author Aleix Lafita
 * @since 4.1.1
 *
 */
public class SecStrucElement {

	private SecStrucType type;
	private int start;
	private int end;

	/**
	 * Create a new SSE object. The start will be the smallest index residue
	 * and the end will be the biggest. The two indices cannot be the same, 
	 * since the SSE will have length 0.
	 * 
	 * @param type
	 * @param res1
	 * @param res2
	 */
	public SecStrucElement(SecStrucType type, int res1, int res2){
		
		this.type = type;
		if (res1 < res2){
			this.start = res1;
			this.end = res2;
		} else if (res2 < res1){
			this.start = res2;
			this.end = res1;
		} else throw new IllegalArgumentException(
				"start and end residues cannot be equal in SecStrucElement");
	}
	
	public SecStrucType getType(){
		return type;
	}
	
	public int getLength(){
		return end-start;
	}
}
