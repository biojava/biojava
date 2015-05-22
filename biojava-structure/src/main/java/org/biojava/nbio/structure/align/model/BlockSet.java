package org.biojava.nbio.structure.align.model;

import java.util.List;

import javax.vecmath.Matrix4d;

/**
 * A BlockSet is a Data Structure to store aligned positions of a multiple alignment as a collection of {@link Block}s.
 * It allows non-sequential alignments and circular permutations, thanks to the multiple {@link Block} format.
 * Every BlockSet has a {@link Pose} object associated, which describes the 3D superimposition of the structures and some values
 * associated with it (RMSD, TMscore, background distances, etc.).
 * A collection of BlockSets, in a {@link MultipleAlignment}, allows the description of flexible multiple alignments.
 * Every BlockSet object is part of a {@link MultipleAlignment} instance, its parent.
 *
 * @author Aleix Lafita
 * 
 */
public interface BlockSet extends ScoresCache {
	
	/**
	 * Creates and returns an identical copy of this blockset, including a deep
	 * clone of all constituent blocks.
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
	 * @param blocks List of alignment Blocks.
	 * @see #getBlocks()
	 */
	public void setBlocks(List<Block> blocks);
	
	/**
	 * Returns a transformation matrix for each structure giving the
	 * 3D superimposition information of the multiple structure alignment.
	 * @return the 3D superimposition information of the alignment
	 * @throws StructureAlignmentException 
	 */
	public List<Matrix4d> getTransformations() throws StructureAlignmentException;
	
	/**
	 * Set a new superposition for the structures.
	 * 
	 * This may trigger other properties to update which depend on the superposition.
	 * @param matrices
	 */
	public void setTransformations(List<Matrix4d> transformations) throws StructureAlignmentException;

	/**
	 * Returns the total number of aligned residues (columns) in the alignment: the sum of all Block lengths.
	 * @return int the total number of aligned residues.
	 * @throws StructureAlignmentException if there are no Blocks.
	 * @see #updateLength()
	 * @see #getCoreLength()
	 * @see #size()
	 * @see #getBlockNum()
	 */
	public int length() throws StructureAlignmentException;
	
	/**
	 * Returns the number of aligned residues (columns) without gaps in the alignment: the sum of all Block core lengths.
	 * @return int the total number of aligned residues.
	 * @throws StructureAlignmentException if there are no Blocks.
	 * @see #updateCoreLength()
	 * @see #length()
	 * @see #size()
	 * @see #getBlockNum()
	 */
	public int getCoreLength() throws StructureAlignmentException;

	/**
	 * Returns the number of aligned structures in the BlockSet.
	 * @return int number of aligned structures
	 * @throws StructureAlignmentException if the BlockSet is empty.
	 * @see #length()
	 * @see #getBlockNum()
	 */
	public int size() throws StructureAlignmentException;
}
