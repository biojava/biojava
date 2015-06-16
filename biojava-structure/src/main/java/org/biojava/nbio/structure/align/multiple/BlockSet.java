package org.biojava.nbio.structure.align.multiple;

import java.util.List;

import javax.vecmath.Matrix4d;

/**
 * A BlockSet is a Data Structure to store a flexible alignment part of a 
 * multiple alignment as a collection of {@link Block}s.
 * <p>
 * It allows non-sequential alignments and circular permutations, 
 * thanks to the multiple {@link Block} format.
 * <p>
 * Every BlockSet has a unique transformation 4D matrix, which describes 
 * the 3D superimposition of the structures in this particular region of the alignment.
 * <p>
 * A collection of BlockSets, in a {@link MultipleAlignment}, allows the description 
 * of multiple alignments with various flexible parts.
 * Every BlockSet object is part of a {@link MultipleAlignment} instance, its parent.
 *
 * @author Aleix Lafita
 * 
 */
public interface BlockSet extends ScoresCache {
	
	/**
	 * Creates and returns an identical copy of this blockset, including a deep
	 * clone of all constituent {@link Block}s.
	 * @return BlockSet identical copy of this object.
	 */
	public BlockSet clone();
	
	/** 
     * Returns the parent MultipleAlignment of the BlockSet.
     * Returns null if there is no referenced object.
     * @return MultipleAlignment the parent MultipleAlignment of the BlockSet, or null.
     * @see #setMultipleAlignment(MultipleAlignment)
     */
	public MultipleAlignment getMultipleAlignment();
	
	/** 
     * Set the back-reference to its parent MultipleAlignment.
     * <p>
     * Neither removes this BlockSet from its previous alignment, if any, nor
     * adds it to the new parent. Calling code should assure that links to
     * and from the ensemble are consistent and free of memory leaks.
     * @param parent the parent MultipleAlignment.
     * @see #getMultipleAlignment()
     */
	public void setMultipleAlignment(MultipleAlignment parent);
	
	/**
	 * Returns the List of alignment Blocks of the BlockSet.
	 * It initializes a new List of Blocks if it is null.
	 * @return List of alignment Blocks.
	 * @see #setBlocks(List)
	 */
	public List<Block> getBlocks();
	
	/**
	 * Set the List of alignment Blocks of the BlockSet.
	 * <p>
	 * Also calls {@link Block#setBlockSet(BlockSet)} for each argument
	 * @param blocks List of alignment Blocks.
	 * @see #getBlocks()
	 */
	public void setBlocks(List<Block> blocks);
	
	/**
	 * Returns a transformation matrix for each structure giving the
	 * 3D superimposition information of the multiple structure alignment.
	 * @return the 3D superimposition information of the alignment
	 */
	public List<Matrix4d> getTransformations();
	
	/**
	 * Set a new superposition for the structures.
	 * <p>
	 * This may trigger other properties to update which depend on the superposition.
	 * @param matrices
	 */
	public void setTransformations(List<Matrix4d> transformations);

	/**
	 * Returns the total number of aligned residues (columns) in the alignment: the sum of all Block lengths.
	 * @return int the total number of aligned residues.
	 * @see #updateLength()
	 * @see #getCoreLength()
	 * @see #size()
	 * @see #getBlockNum()
	 */
	public int length();
	
	/**
	 * Returns the number of aligned residues (columns) without gaps in the alignment: the sum of all Block core lengths.
	 * @return int the total number of aligned residues.
	 * @see #updateCoreLength()
	 * @see #length()
	 * @see #size()
	 * @see #getBlockNum()
	 */
	public int getCoreLength();

	/**
	 * Returns the number of aligned structures in the BlockSet.
	 * @return int number of aligned structures
	 * @see #length()
	 * @see #getBlockNum()
	 */
	public int size();
	
	/**
	 * Clear scores and other properties which depend on the specific alignment.
	 * <p>
	 * This can free memory and ensures consistency for cached variables.
	 */
	public void clear();
}
