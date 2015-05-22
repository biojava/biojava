package org.biojava.nbio.structure.align.model;

import java.util.List;

import javax.vecmath.Matrix4d;

/**
 * A MultipleAlignment is a Data Structure to store the core information of a multiple structure alignment, as a return type.
 * Each alignment is described as a collection of {@link BlockSet}s that define the aligned positions, 
 * a collection of structure identifiers (i,e. Atom arrays), information about the 3D superimposition in {@link Pose},
 * and creation properties (algorithm, version, etc).
 * A collection of MultipleAlignments that share the Atom arrays and creation properties form a {@link MultipleAlignmentEnsemble}.
 * Every MultipleAlignment has a {@link MultipleAlignmentEnsemble} as its parent.
 *
 * @author Aleix Lafita
 * 
 */
public interface MultipleAlignment extends ScoresCache {
	
	/**
	 * Creates and returns an identical copy of this alignment, including a deep
	 * clone of all constituent blocks sets.
	 * @return MultipleAlignment identical copy of this object.
	 */
	public MultipleAlignment clone();
	
	/** 
     * Returns the parent Ensemble of the MultipleAlignment.
     * Returns null if there is no referenced object.
     * @return MultipleAlignmentEnsemble the parent MultipleAlignment of the BlockSet, or null.
     * @see #setParent(MultipleAlignmentEnsemble)
     */
	public MultipleAlignmentEnsemble getEnsemble();
	
	/** 
     * Set the back-reference to its parent Ensemble.
     * @param parent the parent MultipleAlignmentEnsemble.
     * @see #getParent()
     */
	public void setEnsemble(MultipleAlignmentEnsemble parent);

	/**
	 * Returns the BlockSet List of the multiple structure alignment.
	 * Initializes the variable if it is null.
	 * @return List of BlockSets that describe the aligned residues of all the structures.
	 * @see #getBlockSetNum()
	 * @see #setBlockSets(List)
	 */
	public List<BlockSet> getBlockSets();

	/**
	 * Sets the List of BlockSet List of the specified alignment. The optimal alignment is always stored at position 0.
	 * @param blockSets the List of BlockSets that describe the aligned residues.
	 * @see #getBlockSets()
	 */
	public void setBlockSets(List<BlockSet> blockSets);

	/**
	 * Convenience method to get a list of all blocks from all blocksets
	 * @return
	 * @throws StructureAlignmentException
	 */
	public List<Block> getBlocks();

	/**
	 * Returns the List of Strings that represent the multiple sequence alignment of all the structures.
	 * @return List of Strings multiple sequence alignment
	 * @see #updateAlnSequences()
	 */
	public List<String> getAlnSequences();


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
	 * Returns the number of aligned structures in the MultipleAlignment.
	 * @return int number of aligned structures
	 * @throws StructureAlignmentException 
	 * @see #length()
	 * @see #getCoreLength()
	 * @see #getBlockSetNum()
	 */
	public int size() throws StructureAlignmentException;

	/**
	 * Returns the total number of aligned residues (columns) in the multiple alignment: the sum of all BlockSet lengths.
	 * @return int the total number of aligned residues in the alignment.
	 * @throws StructureAlignmentException if there are no BlockSets.
	 * @see #updateLength()
	 * @see #getCoreLength()
	 * @see #size()
	 * @see #getBlockSetNum()
	 */
	public int length() throws StructureAlignmentException;

	/**
	 * Returns the number of aligned residues (columns) without gaps in the alignment: the sum of all BlockSet core lengths.
	 * @return int the total number of aligned residues.
	 * @throws StructureAlignmentException if there are no BlockSets.
	 * @see #updateCoreLength()
	 * @see #length()
	 * @see #size()
	 * @see #getBlockNum()
	 */
	public int getCoreLength() throws StructureAlignmentException;
}
