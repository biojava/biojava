package org.biojava.nbio.structure.align.model;

import java.util.List;

import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.model.Pose.PoseMethod;

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
	 * Returns the global similarity measure among all the structures aligned.
	 * @return double similarity measure.
	 * @see #updateSimilarity()
	 */
	public double getSimilarity();
	
	/**
	 * Calculates and sets the new value of the similarity measure of the alignment.
	 * @see #getSimilarity()
	 */
	public void updateSimilarity();
	
	/**
	 * Returns the coverage of the alignment.
	 * @return double coverage as a value in [0,1].
	 * @see #updateCoverage()
	 */
	public double getCoverage();
	
	/**
	 * Calculates and sets the new value of coverage of the alignment.
	 * @see #getCoverage()
	 */
	public void updateCoverage();
	
	/**
	 * Returns the 3D transformation information of the alignment as a Pose object.
	 * @return Pose the 3D superimposition information.
	 * @see #setPose(Pose)
	 * @see #updatePose(PoseMethod)
	 */
	public Pose getPose();
	
	/**
	 * Set the 3D transformation information of the alignment as a Pose object.
	 * @param pose the Pose instance containing the 3D superimposition information.
	 * @see #updatePose(PoseMethod)
	 */
	public void setPose(Pose pose);
	
	/**
	 * Calculates and sets the new 3D superimposition information in the Pose.
	 * Methods: REFERENCE (align everything to the first structure, the master), 
	 * 			MEDIAN (take the closest structure to all others in average as the master and align everything to it),
	 * 			CONSENSUS (build a consensus structure and align everything to it)
	 * @param method PoseMethod indicating one of the methods listed above, to be used in the superposition.
	 * @throws StructureException 
	 * @see #setPose(Pose)
	 * @see #getPose()
	 */
	public void updatePose(PoseMethod method) throws StructureException;
	
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
