package org.biojava.nbio.structure.align.model;

import java.util.List;

import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.model.Pose.PoseMethod;

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
	 * Returns the 3D transformation information of the alignment as a Pose object.
	 * Initializes a new Pose if it is null.
	 * @return Pose the 3D superimposition information.
	 * @throws StructureAlignmentException 
	 * @see #updatePose(PoseMethod)
	 */
	public Pose getPose() throws StructureAlignmentException;
	
	/**
	 * Calculates and sets the new 3D superimposition information in the Pose of the BlockSet part.
	 * Methods: REFERENCE (align everything to the first structure, the master), 
	 * 			MEDIAN (take the closest structure to all others in average as the master and align everything to it),
	 * 			CONSENSUS (build a consensus structure and align everything to it)
	 * @param method PoseMethod indicating one of the methods listed above, to be used in the superimposition.
	 * @throws StructureException
	 * @throws StructureAlignmentException
	 * @see #getPose()
	 */
	public void updatePose(PoseMethod method) throws StructureException, StructureAlignmentException;
	
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
	 * Calculates and sets the total number of aligned residues (columns) in the alignment: the sum of all Block lengths.
	 * @throws StructureAlignmentException
	 * @see #length()
	 * @see #updateCoreLength()
	 * @see #getCoreLength()
	 */
	public void updateLength() throws StructureAlignmentException;
	
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
	 * Calculates and sets the number of aligned residues (columns) without gaps in the alignment: the sum of all Block core lengths.
	 * @throws StructureAlignmentException
	 * @see #getCoreLength()
	 * @see #length()
	 * @see #updateLength()
	 */
	public void updateCoreLength() throws StructureAlignmentException;
	
	/**
	 * Returns the number of aligned structures in the BlockSet.
	 * @return int number of aligned structures
	 * @throws StructureAlignmentException if the BlockSet is empty.
	 * @see #length()
	 * @see #getBlockNum()
	 */
	public int size() throws StructureAlignmentException;
	
	/**
	 * Returns the number of alignment Blocks in the BlockSet.
	 * @return int number of Blocks
	 * @throws StructureAlignmentException if the BlockSet is empty.
	 * @see #length()
	 * @see #size()
	 */
	public int getBlockNum() throws StructureAlignmentException;
	
	/**
	 * Calls all the update methods for the Cache variables.
	 * @throws StructureAlignmentException
	 * @param method PoseMethod indicating one of the methods listed above, to be used in the superimposition.
	 * @throws StructureException 
	 * @see #updateLength()
	 * @see #updateCoreLength()
	 * @see #updatePose(PoseMethod)
	 */
	public void updateCache(PoseMethod method) throws StructureAlignmentException, StructureException;
	
}
