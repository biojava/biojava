package org.biojava.nbio.structure.align.model;

import java.util.List;

/**
 * A BlockSet is a Data Structure to store the aligned positions of a multiple alignment as a collection of {@link Block}.
 * It allows flexible alignments, non-sequential alignments and circular permutations, thanks to the multiple Block format.
 * It is part of a {@link MultipleAlignment} instance, named as parent.
 *
 * @author Aleix Lafita
 * 
 */
public interface BlockSet extends Cloneable{
	
	/**
	 * Creates and returns an identical copy of this object.
	 * @return BlockSet identical copy of this object.
	 */
	public Object clone();
	
	/** 
     * Returns the parent MultipleAlignment of the BlockSet.
     * Returns null if there is no referenced object. 
     * @return MultipleAlignment the parent MultipleAlignment of the BlockSet, or null.
     * @see #setMultipleAlignment(MultipleAlignment)
     */
	public MultipleAlignment getMultipleAlignment();
	
	/** 
     * Set the back-reference to its parent MultipleAlignment.
     * @param parent the parent MultipleAlignment.
     * @see #getMultipleAlignment()
     */
	public void setMultipleAlignment(MultipleAlignment parent);
	
	/**
	 * Returns the List of alignment Blocks of the BlockSet.
	 * @return List of alignment Blocks.
	 * @see #setBlocks(List)
	 */
	public List<Block> getBlocks();
	
	/**
	 * Set the List of alignment Blocks of the BlockSet.
	 * @param blocks List of alignment Blocks.
	 * @see #getBlocks()
	 */
	public void setBlocks(List<Block> blocks);
	
	/**
	 * Returns the 3D transformation information of the alignment as a Pose object.
	 * @return Pose the 3D transformation information.
	 * @see #setPose(Pose)
	 */
	public Pose getPose();
	
	/**
	 * Set the 3D transformation information of the alignment as a Pose object.
	 * @param pose the Pose instance containing the 3D transformation information.
	 */
	public void setPose(Pose pose);
	
	/**
	 * Returns the total number of aligned residues (columns) in the alignment: the sum of all Block lengths.
	 * @return int the total number of aligned residues.
	 * @see #size()
	 * @see #getBlockNum()
	 */
	public int length();

	/**
	 * Returns the number of aligned structures in the BlockSet.
	 * @return int number of aligned structures
	 * @see #length()
	 * @see #getBlockNum()
	 */
	public int size();
	
	/**
	 * Returns the number of alignment Blocks in the BlockSet.
	 * @return int number of Blocks
	 * @see #length()
	 * @see #size()
	 */
	public int getBlockNum();
	
}
