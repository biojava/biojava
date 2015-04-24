package org.biojava.nbio.structure.align.model;

import java.util.List;

/**
 * A Block is a Data Structure that stores aligned positions of a multiple alignment that fulfill the following conditions:
 * 
 * 		1- Residues are in a sequential order (increasing)
 * 		2- At least two structures have a residue in every column (no empty columns or columns with one residue = no alignment).
 *
 * It is part of a {@link BlockSet} instance, named as parent, which is in turn part of a {@link MultipleAlignment} instance.
 * 
 * @author Aleix Lafita
 * 
 */
public interface Block extends Cloneable{
	
	/**
	 * Creates and returns an identical copy of this object.
	 * @return Block identical copy of this object.
	 */
	public Object clone();
	
	/** 
     * Set the back-reference to its parent BlockSet.
     * @param parent the parent BlockSet.
     * @see #getBlockSet()
     */
	public void setBlockSet(BlockSet parent);

	/** 
     * Returns the parent BlockSet of the Block.
     * Returns null if there is no referenced object. 
     * @return BlockSet the parent BlockSet of the Block, or null.
     * @see #setBlockSet(BlockSet)
     */
	public BlockSet getBlockSet();

	/**
	 * Returns the double List containing the aligned residues for each structure.
	 * alignRes.get(structure).get(residue) = alignRes.get(size).get(length)
	 * @return List a double List of aligned residues for each structure.
	 * @see #setAlignRes()
	 */
	public List<List<Integer>> getAlignRes();

	/**
	 * Set the double List containing the aligned residues for each structure.
	 * @param alignRes a double List of Integers with the aligned residues.
	 * @see #getAlignRes()
	 */
	public void setAlignRes(List<List<Integer>> alignRes);

	/**
	 * Returns the number of aligned positions (columns) in the Block.
	 * @return int number of aligned residues
	 * @see #size()
	 * @see #getAlignRes()
	 */
	public int length();
	
	/**
	 * Returns the number of aligned structures (rows) in the Block.
	 * @return int number of aligned structures
	 * @see #length()
	 * @see #getAlignRes()
	 */
	public int size();
}
