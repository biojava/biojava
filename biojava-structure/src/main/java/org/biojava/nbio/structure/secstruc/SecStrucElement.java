package org.biojava.nbio.structure.secstruc;

import org.biojava.nbio.structure.ResidueNumber;
import org.biojava.nbio.structure.ResidueRange;

/**
 * A secondary structure element (SSE) is an object representing a block
 * of sequential residues that share the same secondary structure type
 * inside a protein chain.
 * 
 * @author Aleix Lafita
 * @since 4.1.1
 *
 */
public class SecStrucElement {

	private SecStrucType type;
	private ResidueRange range;
	private int index;
	private int length;

	/**
	 * Create a new SSE object. 
	 * The start and end residue numbers cannot be the same.
	 * 
	 * @param type object describing the type of SS
	 * @param start first residue of the SSE
	 * @param end final residue of the SSE
	 * @param length number of residues included in the SSE
	 * @param chainID the chain ID
	 */
	public SecStrucElement(SecStrucType type, ResidueNumber start, 
			ResidueNumber end, int length, String chainID){

		this.type = type;
		this.length = length;
		range = new ResidueRange(chainID, start, end);
		if (start.equals(end)) throw new IllegalArgumentException(
				"start and end residues cannot be equal in a SecStrucElement");
	}

	/** 
	 * Returns the {@link SecStrucType} of this element.
	 * @return
	 */ 
	public SecStrucType getType(){
		return type;
	}

	/**
	 * Return the number of residues in the SSE.
	 * @return
	 */
	public int getLength(){
		return length;
	}

	/** 
	 * Returns the ID of this element. 
	 * The ID is the concatenation of the type letter and the 
	 * numerical element identifier (e.g. H1, S1, ...).
	 * @return
	 */
	public String getId() {
		return type.type+index+"";
	}

	/** 
	 * Returns the residue range of this SSE.
	 * @return 
	 */ 
	public ResidueRange getRange() {
		return range;
	}

	public String toString() {
		return getId() + ": " + range.getStart() + " - " + range.getEnd();
	}
	
}
