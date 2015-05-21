package org.biojava.nbio.structure.align.model;

import java.util.List;

/**
 * A Block is a Data Structure that stores aligned positions of a multiple alignment that fulfill the following conditions:
 * 		1- Residues are in a sequential order (increasing or decreasing)
 * 		2- At least two structures have a residue in every column (no empty columns or columns with one residue = no alignment).
 * A collection of Blocks, named {@link BlockSet}, allows the description of circular permutations (CP) and non-sequential
 * alignments.
 * Every Block object is part of a {@link BlockSet} instance, its parent, which has in turn a {@link MultipleAlignment} instance as parent.
 * 
 * @author Aleix Lafita
 * 
 */
public interface Block {
	
	/**
	 * Creates and returns an identical copy of this block.
	 * @return Block identical copy of this object.
	 */
	public Block clone();
	
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
	 * alignRes.get(structure).get(residue) = alignRes.get(size).get(length).
	 * Initializes the variable if it is null.
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
	 * Returns the total number of aligned positions (columns) in the Block.
	 * @return int number of aligned residues.
	 * @see #getCoreLength();
	 * @see #size()
	 */
	public int length();
	
	/**
	 * Returns the number of aligned structures (rows) in the Block.
	 * @return int number of aligned structures.
	 * @see #length()
	 * @see #getCoreLength()
	 */
	public int size();
	
	/**
	 * Returns the number of aligned positions (columns) without gaps in the Block.
	 * @return int number of aligned residues.
	 * @see #updateCoreLength()
	 * @see #length()
	 * @see #size()
	 */
	public int getCoreLength();
	
	/**
	 * Calculates and sets the number of aligned positions without gaps in the Block.
	 * @see #getCoreLength()
	 * @see #length()
	 * @see #size()
	 */
	public int updateCoreLength();
	
}
